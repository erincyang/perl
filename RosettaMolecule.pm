# RosettaMolecule package
# act on a Molecule using Rosetta
#
# 2003-2007, John Karanicolas, Baker group, UW
#
# derived from Molecule

package RosettaMolecule;

require 5.004;

use strict;

use FileHandle;
use IPC::Open2;

use Molecule;

use GenUtil;
use UWutil;
use Sequence;

use vars qw ( @ISA $rosettaExec $rosettaData $nnmakeExec $Vall_name $PD_Exec $PD_ParamExec $PD_MS_Exec $PD_Trans_Exec );

@ISA = ("Molecule");

## data: chain[] -> { id atom[] res[] resinx[]
## data:              xcoor[] ycoor[] zcoor[] }
## data:   atom[] -> { atominx atomname resname resnum chain
## data:               xcoor ycoor zcoor hyd bb seg aux1 aux2
## data:               aux3 }
## data:   res[]  -> { name num start end seg valid chain score des loop }
## molecule information with substructures containing atom and residue
## information as well as a residue lookup table and coordinate cache
## arrays.

## data: chainlookup -> { chainid ... }
## lookup hash table for multiple chains

## data: defchain
## currently selected default chain

## data: segmentlist[] -> { name first last }
## list of segment IDs

## data: ssbond[] -> { chain1 resnum1 chain2 resnum2 }
## disulfide bonds

## data: par -> { protein, series, no_seq, currFile,
## data:          pathfile, outdir, etc.... }
## Rosetta parameters

## data: score -> { score, env, pair, vdw, hs, ss, sheet,
## data:            r-sigma, cb, rg, co, rama, hb_srbb,
## data:            hb_lrbb, etc.... }
## Rosetta score terms

BEGIN {

  if ( ( exists $ENV{'ROSETTAEXEC'} ) && ( $ENV{'ROSETTAEXEC'} ne "" ) ) {
    $rosettaExec = $ENV{'ROSETTAEXEC'};
  }
  else {
    undef $rosettaExec;
  }

  if ( ( exists $ENV{'ROSETTADATA'} ) && ( $ENV{'ROSETTADATA'} ne "" ) ) {
    $rosettaData = $ENV{'ROSETTADATA'};
  }
  else {
    undef $rosettaData;
  }

  if ( ( exists $ENV{'ROSETTAMAKEFRAGS'} )
    && ( $ENV{'ROSETTAMAKEFRAGS'} ne "" ) )
  {
    $nnmakeExec = $ENV{'ROSETTAMAKEFRAGS'};
  }
  else {
    undef $nnmakeExec;
  }

  if ( ( exists $ENV{'VALLNAME'} ) && ( $ENV{'VALLNAME'} ne "" ) ) {
    $Vall_name = $ENV{'VALLNAME'};
  }
  else {
    undef $Vall_name;
  }

  $PD_Exec = "~johnk/bin/PatchDock/patch_dock.Linux";
  $PD_ParamExec = "~johnk/bin/PatchDock/buildParams.pl";
  $PD_MS_Exec = "~johnk/bin/PatchDock/buildMS.pl";
  $PD_Trans_Exec = "~johnk/bin/PatchDock/transOutput.pl";

}

## constructor: new([PDBfile],[parameters])
## creates a new RosettaMolecule object and reads a PDB structures
## if a file name is given

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $farg=shift;
  my @pArr=@_;

  my $self = $class->SUPER::new();

  $self->resetPars(@pArr);

  $self->readPDB($farg)
    if ( defined $farg );

  $self->{par}->{protein} = &_randID();

  return $self;
}

# Destructor
DESTROY {
  my $self = shift;

  &GenUtil::remove( $self->{par}->{protein} . ".pdb" )
    if ( &GenUtil::exists( $self->{par}->{protein} . ".pdb" ) );
}

## method: resetPars( [pars] )
## setup or reset all parameters to their default values. Pars also
## may be passed in, which will be applied after resetting.

sub resetPars {
  my $self = shift;
  my @pArr=@_;

  my %parhash = (

    # General Rosetta parameters
    series => "XX",    # series code

    # Fragment insertion parameters
    no_seq => undef,    # no sequence available yet

    # Current file
    currFile => undef,    # copy of the current molecule if available

    # Currently running job parameters
    pathfile    => "paths.txt",      # name of the Rosetta path file
    outdir      => "./output/",      # output directory
    stdout      => "rosetta.out",    # write Rosetta stdout here
    stderr      => "rosetta.err",    # write Rosetta stderr here
    tempInfiles => "infiles.tmp",    # a list of temporary input files
    outputPDBs  => "output.expd",    # a list of expected output PDBs
    cmdLog      => 1,                # print the Rosetta command-line to the log
    nstruct     => 1,                # number of structures to build
    decoystats  => 1,                # include decoystats flag
    monotone_line_search  => 1,      # include monotone line search flag
    gu_exempt   => undef,            # a file containing a list of group_uns to ignore
    interface_ds  => 0,              # include interface_ds flag
    tight_core  => 0,                # do tight core analysis

    # Loop mode parameters
    loop_fa_input  => 0,             # fa_input flag for loop mode
    loop_fa_output => 1,             # fa_output flag for loop mode
    loop_fa_refine => 0,             # fa_refine flag for loop mode
    loop_fix_natsc => 1,             # fix_natsc flag for loop mode
    cen_rlx        => 0,             # centroid relax
    FA_rlx         => 0,             # full-atom relax

    # Interface mode parameters
    mutfile => "rev.inp",            # mutfile for WT reversion
    revfile => "rev.out",            # output file for WT reversion

    # Design mode parameters
    resfile    => "resfile", # resfile for design mode
    repackDist => 6.0,       # repack within this distance of a designed residue
    interfaceDist  => 5.0,   # design within this distance of an interface
    desock         => 0,     # do iterative design and dockmin
    desrub         => 0,     # do iterative design and backrub (interface only)
    des_int        => 0,     # use interface design mode
    fix_target_seq => 0,     # design to fixed target in interface mode
    water_on_fixed_chain => 0,     # put explicit water on the fixed target chain, in interface mode
    ex_rot         => 0,     # use additional rotamers in repack/design
    favor_native   => undef, # energy by which to favor the native residue
    favor_polar    => undef, # energy by which to favor polar residues
    favor_nonpolar => undef, # energy by which to favor non-polar residues
    favor_charged  => undef, # energy by which to favor charged residues
    favor_aromatic => undef, # energy by which to favor aromatic residues
    favored_PDB_seq => undef, # the name of a PDB file containing the sequence to favor

    sidechain_entropy_penalty => undef, # penalize for each chi angle, to get small sidechains

    gb             => 0,     # use gen_born flag
    soft_rep       => 0,     # use soft_rep flag
    mcmin          => 0,     # use mcmin_trials flag
    design_trials  => 0,     # use design_trials flag
    design_hydrophobic_only => 0,  # allow only hydrophobics in design
    alaScan        => 0,     # use fast_ala_scan flag
    desSwapLoops   => 0,     # swap a specific region with preexisting loops
    ARinfo         => undef, # AR-specific info

    # Docking mode parameters
    dockMin       => 0,      # refine via minimization only
    dockMcm       => 0,      # refine via monte carlo minimization
    dockPert      => undef,  # string of 3 numbers defining perturbation
    dockGenerateDiversity => 0,  # variation of dockpert to generate diverse starting points
    dockRand      => 0,      # randomize, for starting structures
    dockRTmin     => 0,      # RT minimize
    dockConsFname => undef,  # filename containing data for site constraint file

    # Relax mode parameters
    minimize      => 0,      # minimize flag

    # General flags
    multi_chain => 0,        # force multi_chain system

    # Explicit water
    exp_water     => 0,      # explicit water
    use_input_sc  => 1,      # use_input_sc flag
    use_het_water => 0,      # read_hetero_h2o flag

    # jk secret rosetta flags!
    jk_interface => 0,       # interface design weights
    no_trp       => 0,       # disallow trp
    no_new_CG    => 0,       # disallow Cys, Gly unless present in the native
    build_stacks => 0,       # build aromatic stacking
    incorporate_hotspot => 0, # build a Tyr hotspot (new method)
    build_tyr_hotspot => 0,   # build a Tyr hotspot
    refine_hotspot => 0,      # refine a Tyr hotspot
    hotspot_restore_sidechains => "",   # restore the sidechains from the provided "native"
                                        # at the beginning of the hotspot code

    # Charmm min parameters
    chamin_int => 0,         # minimize only interface atoms

    # For reversion to WT
    revertScoreCut => 0.5,    # revert if difference is less than this

    # Skip the cleaning step
    doClean => 0,             # clean the output files (fit to template)

    # Renumber across chains
    doRenumber => 1,          # renumbering better if calling Rosetta

    junk => undef             # junk
  );

  $self->{par} = \%parhash;

  $self->_getpar(@pArr);

  return;
}

## method: readPDB(PDBfile)
## reads a protein structures from a PDB file.
## undef $self if the file is empty

sub readPDB {
  my $self=shift;
  my $PDBname=shift;

  if ( !&GenUtil::exists($PDBname) ) {
    undef $self;
    printf STDERR "Warning: PDBfile %s does not exist\n", $PDBname;
    return;

  }
  elsif ( !&GenUtil::zexists($PDBname) ) {
    undef $self;
    printf STDERR "Warning: PDBfile %s is empty\n", $PDBname;
    return;

  }
  else {

    # Ensure there's a "TER" or an "END"
    # (if this is not present, a Rosetta job probably wrote a partial PDB then stopped)
    my $inpdb = &GenUtil::getInputFile($PDBname);
    my $done=0;
    while (<$inpdb>) {
      my $inline = $_;
      if ( !$done ) {
        chomp($inline);
        $done = 1 if ( ( $inline =~ /^TER/ ) || ( $inline =~ /^END/ ) );
      }
    }
    undef $inpdb;
    if ( !$done ) {
      undef $self;
      return;
    }

    $self->SUPER::readPDB($PDBname,ignoressbond=>1,ignorehet=>1);
#    $self->SUPER::readPDB($PDBname,ignoressbond=>1,ignorehet=>0);

#    $self->removeHetero();
    $self->findSSBonds();

    $self->renumberAcrossChains() if ( $self->{par}->{doRenumber} );

    $self->resetResidueNameAll( "HSD", "HIS" );
    $self->resetResidueNameAll( "HSE", "HIS" );

    $self->updateScore($PDBname);

    my $c = $self->getChain();
    if ( ( !defined $c->{id} ) || ( $c->{id} eq "-" ) ) {
      $self->setChain("A");
    }

    $self->{par}->{currFile} = $PDBname;

    # How to do design and loop modes (ie. don't)
    foreach my $c ( @{ $self->activeChains() } ) {
      if ( $#{ $c->{res} } >= 0 ) {
        for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
          $c->{res}->[$ir]->{des} = "NATRO"
            if ( !defined $c->{res}->[$ir]->{des} );
          $c->{res}->[$ir]->{loop} = 0
            if ( !defined $c->{res}->[$ir]->{loop} );
        }
      }
    }
  }

  return;
}

## method: clone()
## clone the current molecule (ie. copy constructor)

sub clone {
  my $self = shift;
  my $valid = shift;

  # Clone the Molecule info
  my $n = $self->SUPER::clone($valid);

  # Clone the score hash
  $n->{score}={};
  %{$n->{score}}=%{$self->{score}};

  # Clone the par hash
  $n->{par}={};
  %{$n->{par}}=%{$self->{par}};

  # Undefine currFile if only part of the system is cloned
  undef $self->{par}->{currFile} if ($valid);

  return $n;
}

## method: renumberAcrossChains()
## polymorphism for renumberAcrossChains method in Molecule class

sub renumberAcrossChains {
  my $self = shift;
  my @args = @_;

  return if ( !defined $self->activeChains() );

  undef $self->{par}->{currFile};

  return $self->SUPER::renumberAcrossChains(@args);
}

## method: renumber()
## polymorphism for renumber method in Molecule class

sub renumber {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};

  return $self->SUPER::renumber(@args);
}

## method: numberReset()
## polymorphism for renumber method in Molecule class

sub numberReset {
  my $self = shift;
  my @args = @_;

  if ( $self->{par}->{doRenumber} ) {
      &GenUtil::log("Could not apply number reset when doRenumber is true");
      die "Could not apply number reset when doRenumber is true";
  } else {
      undef $self->{par}->{currFile};
      return $self->SUPER::numberReset(@args);
  }

  return;
}

## method: shiftResNumber()
## polymorphism for shiftResNumber method in Molecule class

sub shiftResNumber {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};

  if ( $self->{par}->{doRenumber} ) {
      $self->renumberAcrossChains();
  } else {
      return $self->SUPER::shiftResNumber(@args);
  }

  return;
}

## method: updateScore(PDBfile, [getByRes], [getVsPDB])
## read scores from PDBfile
## note: B-factors will not be updated

sub updateScore {
  my $self=shift;
  my $PDBname=shift;
  my $getByRes=shift;
  my $getVsPDB=shift;

  $getByRes = 0 if ( !defined $getByRes );
  $getVsPDB = 0 if ( !defined $getVsPDB );

  if ( !&GenUtil::exists($PDBname) ) {
    printf STDERR "Warning: PDBfile %s for updateScore does not exist\n",
      $PDBname;
    return;

  }
  elsif ( !&GenUtil::zexists($PDBname) ) {
    printf STDERR "Warning: PDBfile %s for updateScore is empty\n", $PDBname;
    return;

  }
  else {

    # Get possible score info
    my $inpdb = &GenUtil::getInputFile($PDBname);
    undef $self->{score};
    $self->{score} = {};
    $self->{uns} = ();
    $self->{gu} = ();
    my $gotScoresByRes = 0;
    my $gotScoresVsPDB = 0;
    while (<$inpdb>) {

      if (/:/) {
        chomp();
        my @p = split( /[:|,| ]+/, $_ );
        shift(@p) if ( $p[0] eq "" );
        while ( $#p >= 0 ) {
          my $k = shift(@p);
          if ( ( defined $p[0] ) && ( $p[0] !~ /[:|,]/ ) ) {
            $self->{score}->{$k} = ( shift(@p) + 0.0 );
          }
          else {
            $self->{score}->{$k} = undef;
          }
        }

      }
      elsif (/UNS /) {
        chomp();
        my @p = split( ' ', $_ );
        my $newds = {};
        $newds->{tag} = $p[0];
        $newds->{type} = $p[1];
        $newds->{restype} = $p[4];
        $newds->{resnum} = $p[5];
        $newds->{atomname} = $p[10];
        $newds->{SASA1} = $p[11];
        $newds->{SASA2} = $p[12];
        $newds->{SASA3} = $p[13];
        push( @{ $self->{uns} }, $newds );

      }
      elsif (/GU /) {
	  chomp();
	  my @p = split( ' ', $_ );
	  my $newgu = {};

	  my $offset=0;
	  $offset=1 if ($p[0] ne "GU");
	  $newgu->{tag} = undef;
	  $newgu->{tag} = $p[0] if ($offset > 0);
	  $newgu->{restype} = $p[$offset+1];
	  $newgu->{type} = $p[$offset+2];
	  $newgu->{resnum} = $p[$offset+3];
	  push( @{ $self->{gu} }, $newgu );

      }
      elsif ( (/^res aa /) && ( !$gotScoresByRes ) && ($getByRes) ) {

        my $done = 0;
        while ( !$done ) {
          my $inline = <$inpdb>;
          if ( $inline =~ /tot/ ) {
            $done = 1;
          }
          else {
            chomp($inline);
            my @p = split( ' ', $inline );
            my $resnum = $p[0];
            foreach my $c ( @{ $self->activeChains() } ) {
              if ( $#{ $c->{res} } >= 0 ) {
                my $r = $self->getResidueInChain( $resnum, $c );
                if ( defined $r ) {
                  if ( $#p >= 0 ) {
                    $r->{score}->{fa_atr} = $p[2];
                    $r->{score}->{fa_rep} = $p[3];
                    $r->{score}->{LJ} = $p[2] + $p[3];
                    $r->{score}->{fa_sol} = $p[4];
                    $r->{score}->{dun} = $p[7];
                    $r->{score}->{HB} = $p[9];
                    $r->{score}->{score} = $p[$#p];
                  }
                }
              }
            }
          }
        }
        $gotScoresByRes = 1;

      }
      elsif ( (/^ *energies-average/) && ( !$gotScoresVsPDB ) && ($getVsPDB) ) {
        if ( <$inpdb> =~ /^ *res  aa / ) {
          my $done = 0;
          while ( !$done ) {
            my $inline = <$inpdb>;
            if ( $inline =~ /tot/ ) {
              chomp($inline);
              my @p = split( ' ', $inline );
              $self->{VsPDB} = {};
              $self->{VsPDB}->{LJatr} = $p[1];
              $self->{VsPDB}->{LJrep} = $p[2];
              $self->{VsPDB}->{sol} = $p[3];
              $self->{VsPDB}->{HB} = $p[7];
              $self->{VsPDB}->{LJ} = $p[9];
              $self->{VsPDB}->{tot} = $p[10];
              $self->{VsPDB}->{SASApack} = $p[11];
              $self->{VsPDB}->{SASAprob} = $p[12];
              $done = 1;
            }
            else {
              chomp($inline);
              my @p = split( ' ', $inline );
              my $resnum = $p[0];
              foreach my $c ( @{ $self->activeChains() } ) {
                if ( $#{ $c->{res} } >= 0 ) {
                  my $r = $self->getResidueInChain( $resnum, $c );
                  if ( ( defined $r ) && ( $#p >= 0 ) ) {
                    $r->{VsPDB}->{LJatr} = $p[3];
                    $r->{VsPDB}->{LJrep} = $p[4];
                    $r->{VsPDB}->{sol} = $p[5];
                    $r->{VsPDB}->{dun} = $p[7];
                    $r->{VsPDB}->{HB} = $p[9];
                    $r->{VsPDB}->{LJ} = $p[11];
                    $r->{VsPDB}->{tot} = $p[12];
                    $r->{VsPDB}->{SASApack} = $p[13];
                    $r->{VsPDB}->{SASAprob} = $p[15];
                  }
                }
              }
            }
          }
        }
        $gotScoresVsPDB = 1;

      }
    }
    undef $inpdb;

  }

  return;
}

## method: updateBfactors(PDBfile)
## update B-factors from a PDB file.
## Note: expects residue numbering to match

sub updateBfactors {
  my $self = shift;
  my $PDBname = shift;

  my $Bfac = {};

  my $mol = RosettaMolecule->new($PDBname);
  foreach my $rc ( @{ $mol->activeChains() } ) {
    foreach my $a ( @{ $rc->{atom} } ) {
      my $key = $a->{resnum} . $a->{atomname};
      $Bfac->{$key} = $a->{aux2};
    }
  }
  undef $mol;

  foreach my $c ( @{ $self->activeChains() } ) {
    foreach my $a ( @{ $c->{atom} } ) {
      my $key = $a->{resnum} . $a->{atomname};
      if ( defined $key ) {
        $a->{aux2} = $Bfac->{$key};
      }
      else {
        printf STDERR "Key %s not found\n", $key;
      }
    }
  }

  return;
}

## method: writePDB(PDBfile,[full | coor])
## write PDB file
## "full" mode includes score info, default "coor" mode doesn't.
## PDBfile is simply copied from existing PDBfile if possible,
## but specification of mode forces writing from current

sub writePDB {
  my $self = shift;
  my $PDBname = shift;
  my $mode = shift;
  my %wpar = @_;

  my $compress = $wpar{compress};
  $compress = 1 if ( !defined $compress );

  $PDBname = $self->{par}->{protein} . ".pdb" if ( !defined $PDBname );

  if ( (! defined $mode ) && ( defined $self->{par}->{currFile} ) &&
       ( -e $self->{par}->{currFile} ) ) {
    &GenUtil::copy( $self->{par}->{currFile}, $PDBname )
      if ( $PDBname ne $self->{par}->{currFile} );
  }
  else {

    $mode = "coor" if ( !defined $mode );

    foreach my $c ( @{ $self->activeChains() } ) {
      if ( $#{ $c->{res} } >= 0 ) {
        for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
          my $r = $c->{res}->[$ir];
          $self->resetResidueName( "HIS", $r->{num}, $c->{id} )
            if ( $r->{name} =~ /HSD|HSE|HSP/ );
        }
      }
    }

    $self->SUPER::writePDB($PDBname,"translate"=>"GENERIC","compress"=>$compress,"TER"=>1);

    if ( lc($mode) eq "full" ) {
      my $junkfile = "junk.pdb";
      &GenUtil::rename( $PDBname, $junkfile );
      my $cleanPDB = &GenUtil::getInputFile($junkfile);
      my $outPDB = &GenUtil::getOutputFile($PDBname);
      while (<$cleanPDB>) {
        if ( (/^ATOM/) || (/^TER/) || (/^SSBOND/) ) {
          chomp;
          printf $outPDB "%s\n", $_;
        }
      }
      undef $cleanPDB;
      &GenUtil::remove($junkfile);

      printf $outPDB "TER\n";
      foreach my $k ( keys %{ $self->{score} } ) {
        printf $outPDB "%s:  %f\n", $k, $self->{score}->{$k};
      }
      printf $outPDB "END\n";
      undef $outPDB;

      &GenUtil::compress($PDBname) if ($compress);

    }
  }

  return;
}

## function: setRosettaExec
## Set the Rosetta executable to be used
## (overrides the environment variable ROSETTAEXEC)

sub setRosettaExec {
  $rosettaExec = shift;
  return;
}

## function: reportRosettaExec
## Report the current Rosetta executable

sub reportRosettaExec {
  return $rosettaExec;
}

## function: useSvnRosettaExec
## Use the Rosetta executable from env var SVN_ROSETTAEXEC
## (overrides the environment variable ROSETTAEXEC)

sub useSvnRosettaExec {
  if ( ( exists $ENV{'SVN_ROSETTAEXEC'} ) && ( $ENV{'SVN_ROSETTAEXEC'} ne "" ) )
  {
    $rosettaExec = $ENV{'SVN_ROSETTAEXEC'};
  }
  else {
    die "SVN_ROSETTAEXEC environment variable has not been set";
  }
  return;
}

## function: setRosettaData
## Set the Rosetta data directory to be used
## (overrides the environment variable ROSETTADATA)

sub setRosettaData {
  $rosettaData = shift;
  return;
}

## method: setBfactorsResScore()
## set Bfactors to a residue score field, if available

sub setBfactorsResScore {
  my $self = shift;
  my $field = shift;

  $field = "score" if ( !defined $field );

  die "No chains defined for setBfactorsResScore"
    if ( !defined $self->activeChains() );

  my $VsPDB = 0;
  if ( uc($field) =~ /VSPDB/ ) {
    $VsPDB = 1;
    ( $field = $field ) =~ s/VSPDB//gi;
  }

  my $EresList = ();
  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      my $a = $c->{atom};
      for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
        my $r = $c->{res}->[$ir];
        if ( ( !$VsPDB ) && ( defined $r->{score}->{$field} ) ) {
          my $Eres = $r->{score}->{$field};
          for ( my $ai = $r->{start} ; $ai <= $r->{end} ; $ai++ ) {
            my $atom = $a->[$ai];
            $atom->{aux2} = $Eres;
            if ( $atom->{atomname} eq "CA" ) {
              my $f = {};
              $f->{resnum} = $r->{num};
              $f->{val} = $Eres;
              push( @{$EresList}, $f );
            }
          }

        }
        elsif ( ($VsPDB) && ( defined $r->{VsPDB}->{$field} ) ) {
          my $Eres = $r->{VsPDB}->{$field};
          for ( my $ai = $r->{start} ; $ai <= $r->{end} ; $ai++ ) {
            my $atom = $a->[$ai];
            $atom->{aux2} = $Eres;
            if ( $atom->{atomname} eq "CA" ) {
              my $f = {};
              $f->{resnum} = $r->{num};
              $f->{val} = $Eres;
              push( @{$EresList}, $f );
            }
          }

        }
        else {
          die "Did not recognize field " . $field;
        }
      }
    }
  }

  undef $self->{par}->{currFile};
  return $EresList;
}

## method: setBfactorsRadii()
## set Bfactors to the Rosetta atomic radii
## note: radii taken from etable.h

sub setBfactorsRadii {
  my $self = shift;

  die "No chains defined for setBfactorsRadii"
    if ( !defined $self->activeChains() );

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{atom} } >= 0 ) {
      foreach my $a ( @{ $c->{atom} } ) {
        my $rad = 1.0;
        my $type = $a->{atomname};

        if ( $type =~ /^C/ ) {
          $rad = 2.0;
        }
        elsif ( $type =~ /^N/ ) {
          $rad = 1.75;
        }
        elsif ( $type =~ /^O/ ) {
          $rad = 1.55;
        }
        elsif ( $type =~ /^S/ ) {
          $rad = 1.9;
        }
        elsif ( $type eq "H" ) {
          $rad = 1.0;
        }
        elsif ( $type =~ /^.W/ ) {
          $rad = 1.4;
        }
        else {

          my $resname = $a->{resname};
          if ( $resname eq "ARG" ) {
            if ( ( $type eq "HE" ) || ( $type =~ /HH/ ) ) {
              $rad = 1.0;
            }
            else {
              $rad = 1.2;
            }
          }
          elsif ( $resname eq "LYS" ) {
            if ( $type =~ /HZ$/ ) {
              $rad = 1.0;
            }
            else {
              $rad = 1.2;
            }
          }
          elsif ( $resname eq "ASN" ) {
            if ( $type =~ /HD2$/ ) {
              $rad = 1.0;
            }
            else {
              $rad = 1.2;
            }
          }
          elsif ( $resname eq "GLN" ) {
            if ( $type =~ /HE2$/ ) {
              $rad = 1.0;
            }
            else {
              $rad = 1.2;
            }
          }
          elsif ( ( $resname eq "SER" ) || ( $resname eq "THR" ) ) {
            if ( $type =~ /^HG/ ) {
              $rad = 1.0;
            }
            else {
              $rad = 1.2;
            }
          }
          elsif ( $resname eq "TYR" ) {
            if ( $type eq "HH" ) {
              $rad = 1.0;
            }
            else {
              $rad = 1.2;
            }
          }
          elsif ( $resname eq "HIS" ) {
            if ( ( $type eq "HD1" ) || ( $type eq "HE2" ) ) {
              $rad = 1.0;
            }
            else {
              $rad = 1.2;
            }
          }
          elsif ( $resname eq "TRP" ) {
            if ( $type eq "HE1" ) {
              $rad = 1.0;
            }
            else {
              $rad = 1.2;
            }
          }
          else {

            # Non-polar and aromatic H's
            $rad = 1.2;
          }
        }

        $a->{aux2} = $rad;
      }
    }
  }

  undef $self->{par}->{currFile};
  return;
}

## method: setBfactorsFromAux3()
## polymorphism

sub setBfactorsFromAux3 {
  my $self = shift;

  $self->SUPER::setBfactorsFromAux3();
  undef $self->{par}->{currFile};

  return;
}

## method: operateOnAllAtoms(operator)
## polymorphism

sub operateOnAllAtoms {
  my $self = shift;
  my $operator = shift;

  $self->SUPER::operateOnAllAtoms($operator);
  undef $self->{par}->{currFile};

  return;
}

## method: operateOnAllAtomsVsRef(refMol, operator)
## polymorphism

sub operateOnAllAtomsVsRef {
  my $self = shift;
  my $refMol = shift;
  my $operator = shift;

  $self->SUPER::operateOnAllAtomsVsRef($refMol,$operator);
  undef $self->{par}->{currFile};

  return;
}

## method: residueSASA(res,scOnly)
## return the total SASA for a residue
## Note: assumes the B-factor field contains the atomic SASA

sub residueSASA {
  my $self = shift;
  my $res = shift;
  my $scOnly = shift;

  my $cid = $res->{chain};
  my $start = $res->{start};
  my $end = $res->{end};

  my $SASA = 0.0;
  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $c->{id} eq $cid ) {
      for ( my $ai = $start ; $ai <= $end ; $ai++ ) {
        my $atom = $c->{atom}->[$ai];
        if ( ( !$scOnly ) || ( !$atom->{bb} ) ) {
	  # Don't include negative SASA values
          $SASA += $atom->{aux2} if ($atom->{aux2} >= 0.0);
        }
      }
    }
  }

  return $SASA;
}

## method: numDesRes()
## return a hash of the total number of residues of each design type

sub numDesRes {
  my $self = shift;

  my $nr = {};

  return $nr
    if ( !defined $self->activeChains() );

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      my $cid = $c->{id};
      for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
        my $type = $c->{res}->[$ir]->{des};
        $nr->{$type} = 0 if ( !defined $nr->{$type} );
        $nr->{$type}++;
      }
    }
  }

  return $nr;
}

## method: getScore()
## return the current score if defined

sub getScore {
  my $self = shift;

  return $self->{score}->{score}
    if ( ( defined $self->{score} )
    && ( defined $self->{score}->{score} ) );
  return 9999.9;
}

## method: getScoreHash()
## return the current score hash

sub getScoreHash {
  my $self = shift;

  return $self->{score} if ( defined $self->{score} );
  my $null = {};
  $null->{score} = 999.9;
  return $null;
}

## method: getRegionScore(startRes,stopRes)
## return the total score of a certain region

sub getRegionScore {
  my $self = shift;
  my $startRes = shift;
  my $stopRes = shift;

  return undef
    if ( !defined $self->activeChains() );

  $startRes = 1 if ( !defined $startRes );
  $stopRes = $self->numres() if ( !defined $stopRes );

  my $E = 0.0;
  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
        $E += $c->{res}->[$ir]->{score}
          if ( ( $c->{res}->[$ir]->{num} >= $startRes )
          && ( $c->{res}->[$ir]->{num} <= $stopRes ) );
      }
    }
  }

  return $E;
}

## method: getChainScore(cid)
## return the total score of a certain chain

sub getChainScore {
  my $self = shift;
  my $cid = shift;

  my $c = $self->getChain($cid);
  die "Cannot find chain %s", $cid if ( !defined $c->{res} );

  my $E = 0.0;
  if ( $#{ $c->{res} } >= 0 ) {
    for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
      $E += $c->{res}->[$ir]->{score};
    }
  }

  return $E;
}

## method: need_Rosetta(earr)
## Given a list of expected filenames,
## determine whether to call Rosetta

sub need_Rosetta {
  my $self = shift;
  my $earr = shift;

  my $unpacked_files = $self->unpack_concat_output();

  # Look for any expected output files which are not present
  my $need_file = 0;
  foreach my $efile ( @{$earr} ) {
    $need_file = 1;
    if ( ( defined $unpacked_files->{$efile} ) ||
	 ( defined $unpacked_files->{$efile . ".gz"} ) ){
      $need_file = 0;
    } elsif ( &GenUtil::exists($efile) ) {
      $need_file = 0;
    } else {
      my $prefix=&PDButil::prefixAndPath($efile);
      my $restored_fname = $prefix.".restored.pdb";
      if ( &GenUtil::exists($restored_fname) ) {
	$need_file = 0;
      }
    }
    return 1 if ( $need_file );
  }

  return 0;
}

## function: unpack_concat_file()
## Given a list of expected filenames,
## determine whether to call Rosetta

sub unpack_concat_file {
  my $concat_fname = shift;

  return if ( ! &GenUtil::exists($concat_fname) );

  &GenUtil::log("unpacking ".$concat_fname);

  my $cfile = &GenUtil::getInputFile($concat_fname);
  my $pdbfile;
  my $new_output_file=0;
  my $pdb_fname="-";
  while (<$cfile>) {
    my $inline = $_;
    chomp($inline);
    if ( $inline =~ /^HEADER  CURRENT FILENAME/) {
      undef $pdbfile;
      &GenUtil::compress( $pdb_fname ) if ( $pdb_fname ne "-");
      my @p = split(' ',$inline);
      $pdb_fname=$p[3];
      ( $pdb_fname = $pdb_fname ) =~ s/.gz$//gi;
      $new_output_file=1;
    } else {
      if ( $new_output_file ) {
	$pdbfile = &GenUtil::getOutputFile( $pdb_fname );
	$new_output_file=0;
      }
      printf $pdbfile "%s\n", $inline;
    }
  }
  undef $pdbfile;
  &GenUtil::compress( $pdb_fname ) if ( $pdb_fname ne "-");
  undef $cfile;

  &GenUtil::remove($concat_fname);

  return;

}

## method: unpack_concat_output()
## Given a list of expected filenames,
## determine whether to call Rosetta

sub unpack_concat_output {
  my $self = shift;

  my $outdir = $self->{par}->{outdir};
  my $unpacked_files = {};
  &unpack_concat_file($outdir . "outpdbs");
  my $num_processes = &UWutil::reportClusterCpus();
  $num_processes = 100 if ( $num_processes == 1 );
  for ( my $i = 0 ; $i < $num_processes ; $i++ ) {
    &unpack_concat_file($outdir . $i . "/outpdbs");
    if ( &GenUtil::exists($outdir . $i . "/outlist.txt") ) {
      my $listfile = &GenUtil::getInputFile($outdir . $i . "/outlist.txt");
      while (<$listfile>) {
	my $inline = $_;
	chomp($inline);
	( $inline = $inline ) =~ s/^\.\///gi;
	# Note: mark all generated files as unpacked
	$unpacked_files->{$inline} = 1;
      }
      undef $listfile;
    }
  }

  return $unpacked_files;

}

## method: null_call([$startArr],[parameters])
## Setup input files as if they're output files
## Note: no call to Rosetta

sub null_call {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $startArr = $self->_setupFromStartArr($startArr);

  my $tempInfiles = &GenUtil::getOutputFile( $self->{par}->{tempInfiles} );

  my $filelist;
  my $earr;
  my $num_processes;

  ($filelist,$earr,$num_processes) =
      $self->_setupOutfiles($startArr,0,1);
  printf $tempInfiles "%s\n", $filelist;

  $self->_rosettaSetup($tempInfiles);

  undef $tempInfiles;

  # Link output files from input files
  &GenUtil::log("linking files in null_call");
  my $ename = $self->{par}->{outdir} . $self->{par}->{outputPDBs};
  my $elist = &GenUtil::getInputFile($ename);
  while (<$elist>) {
    my $inline = $_;
    chomp($inline);
    my $f = {};
    ( $f->{filename}, $f->{Rosout}, $f->{origfile} ) =
      split( ' ', $inline );
    if ( ( $f->{origfile} =~ /gz$/ ) && ( $f->{Rosout} !~ /gz$/ ) ) {
      $f->{Rosout} .= ".gz";
    }
    &GenUtil::link( $f->{origfile}, $f->{Rosout} );
  }
  undef $elist;

  return $self->_processRosettaOutput();
}

## method: computeScore([$startArr],[parameters])
## Run Rosetta score mode

sub computeScore {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $startArr = $self->_setupFromStartArr($startArr);

  my $tempInfiles = &GenUtil::getOutputFile( $self->{par}->{tempInfiles} );

  my $cmd = $rosettaExec;
  $cmd .= " " . $self->{par}->{series} . " " . $self->{par}->{protein} . " A";

  $cmd .= " -score";
  $cmd .= " -fa_input";
  $cmd .= " -fa_output";

  my $filelist;
  my $earr;
  my $num_processes;

  ($filelist,$earr,$num_processes) =
      $self->_setupOutfiles($startArr,0);
  $cmd .= " -l " . $filelist;
  printf $tempInfiles "%s\n", $filelist;

  $cmd .= $self->_rosettaSetup($tempInfiles);

  undef $tempInfiles;

  if ( !$self->need_Rosetta($earr) ) {
    &GenUtil::log("no need to call Rosetta in score mode");
    return $self->_processRosettaOutput();
  }

  &GenUtil::log("calling Rosetta in score mode");
  if ( $self->{par}->{cmdLog} ) {
    &GenUtil::log( "Rosetta cmd: " . $cmd );
  }

  if ( &UWutil::reportClusterCpus() > 0 ) {
    &UWutil::runCmdCluster($cmd, $num_processes, $self->{par}->{outdir}, $self->{par}->{stdout}, $self->{par}->{stderr}, 1);
    sleep(10);
    &UWutil::waitForClusterCompletion();
  }
  else {
    &UWutil::runCmd($cmd, $self->{par}->{stdout}, $self->{par}->{stderr});
  }

  return $self->_processRosettaOutput();
}

## method: idealize([$startArr],[parameters])
## Run Rosetta idealize mode
## Idealize is for only the current structure,
## unless a pointer to an array of
## RosettaMolecules is passed. In this case,
## the current is used as a template and designs
## are carried out for the incoming RosettaMolecules.

sub idealize {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $startArr = $self->_setupFromStartArr($startArr);

  # jk doesn't work for idealize due to unusual energy landscape
  $self->{par}->{monotone_line_search} = 0;

  my $tempInfiles = &GenUtil::getOutputFile( $self->{par}->{tempInfiles} );

  my $cmd = $rosettaExec;
  $cmd .= " " . $self->{par}->{series} . " " . $self->{par}->{protein} . " A";

  $cmd .= " -idealize";
  $cmd .= " -constant_seed"; # not sure why this is needed...
  $cmd .= " -fa_input -fa_output";

  my $filelist;
  my $earr;
  my $num_processes;

  ($filelist,$earr,$num_processes) =
      $self->_setupOutfiles($startArr,0);
  $cmd .= " -l " . $filelist;
  printf $tempInfiles "%s\n", $filelist;

  $cmd .= $self->_rosettaSetup($tempInfiles);

  undef $tempInfiles;

  if ( !$self->need_Rosetta($earr) ) {
    &GenUtil::log("no need to call Rosetta in idealize mode");
    return $self->_processRosettaOutput();
  }

  &GenUtil::log("calling Rosetta in idealize mode");
  if ( $self->{par}->{cmdLog} ) {
    &GenUtil::log( "Rosetta cmd: " . $cmd );
  }

  if ( &UWutil::reportClusterCpus() > 0 ) {
    &UWutil::runCmdCluster($cmd, $num_processes, $self->{par}->{outdir}, $self->{par}->{stdout}, $self->{par}->{stderr}, 1);
    sleep(10);
    &UWutil::waitForClusterCompletion();
  }
  else {
    &UWutil::runCmd($cmd, $self->{par}->{stdout}, $self->{par}->{stderr});
  }

  return $self->_processRosettaOutput();
}

## method: design([$startArr],[parameters])
## Run Rosetta design mode
## Design is for only the current structure,
## unless a pointer to an array of
## RosettaMolecules is passed. In this case,
## the current is used as a template and designs
## are carried out for the incoming RosettaMolecules.

sub design {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $startArr = $self->_setupFromStartArr($startArr);

  my $tempInfiles = &GenUtil::getOutputFile( $self->{par}->{tempInfiles} );

  my $cmd = $rosettaExec;
  $cmd .= " " . $self->{par}->{series} . " " . $self->{par}->{protein} . " A";

  $cmd .= " -design";
  $cmd .= " -output_coord";

  $cmd .= " -no_trp"        if ( $self->{par}->{no_trp} );
  $cmd .= " -mcmin_trials " if ( $self->{par}->{mcmin} );
  $cmd .= " -design_trials " if ( $self->{par}->{design_trials} );

  # disallow Cys,Gly at positions which were not Cys,Gly previously
  $cmd .= " -no_new_CG" if ( $self->{par}->{no_new_CG} );

  # flag for building aromatic stacking
  $cmd .= " -build_stacks" if ( $self->{par}->{build_stacks} );

  # flag for building a Tyr hotspot
  $cmd .= " -build_tyr_hotspot" if ( $self->{par}->{build_tyr_hotspot} );

  # flag for refining a Tyr hotspot
  $cmd .= " -refine_hotspot -hotspot_ene_scaling 1.2"
      if ( $self->{par}->{refine_hotspot} );

  $cmd .= " -design_hydrophobic_only"
      if ( $self->{par}->{design_hydrophobic_only} );

  # AR specific information
  $cmd .= " -AR_conservation_file " . $self->{par}->{ARinfo}
      if ( defined $self->{par}->{ARinfo} );

  if ( ! $self->{par}->{incorporate_hotspot} ) {
    $cmd.=" -favor_native_residue ".$self->{par}->{favor_native}
      if ( defined $self->{par}->{favor_native} );
    $cmd .= " -favor_polar ".$self->{par}->{favor_polar}
      if ( defined $self->{par}->{favor_polar} );
    $cmd .= " -favor_nonpolar ".$self->{par}->{favor_nonpolar}
      if ( defined $self->{par}->{favor_nonpolar} );
    $cmd .= " -favor_charged ".$self->{par}->{favor_charged}
      if ( defined $self->{par}->{favor_charged} );
    $cmd .= " -favor_aromatic ".$self->{par}->{favor_aromatic}
      if ( defined $self->{par}->{favor_aromatic} );
    $cmd .= " -favor_PDB_sequence ".$self->{par}->{favored_PDB_seq}
      if ( defined $self->{par}->{favored_PDB_seq} );
  }

  if ( $self->{par}->{fixInterfaceUNS} ) {
    $cmd .= " -fix_interface_UNS";
    $self->setParameter(identifyInterfaceUNS=>1);
  }

  if ( $self->{par}->{incorporate_hotspot} ) {

    $cmd .= " -incorporate_hotspot";
    $cmd .= " -hotspot_on_AR";
    $cmd .= " -set_interface_cutoff 7.0";
#    $cmd .= " -HS_WRITE_ALL";
    $cmd .= " -HS_NO_PRESERVE_IG";
    $cmd .= " -fa_input";
    $cmd .= " -fa_output";

    if ( $self->{par}->{dockPert} ne "") {
      $cmd .= " -dock_pert ".$self->{par}->{dockPert};
      $cmd .= " -use_score12";
      $cmd .= " -no_filters";
      $cmd .= " -norepack1 -norepack2";
      if ( defined $self->{par}->{dockConsFname} ) {
	my $cstFname = $self->{par}->{protein} . ".cst";
	$self->setupDockingCst($cstFname);
	printf $tempInfiles "%s\n", $cstFname;
      }
      $cmd .= " -atom_vdw_set highres";
      #  $cmd.=" -atom_vdw_set hybrid";
    }

    $cmd .= " -hotspot_restore_sidechains ".$self->{par}->{hotspot_restore_sidechains}
      if ( $self->{par}->{hotspot_restore_sidechains} ne "");

#    $cmd .= " -favor_aromatic -0.2";
#    $cmd .= " -use_hs_bbind_conformer hs_conf.txt";


# jk TESTING
#    my @AR_loop_flist = glob("~johnk/new_interface/AR_loop_library/*.aligned.pdb");
#    if ( $#AR_loop_flist >= 0 ) {
#      $cmd .= " -AR_module_list AR_loop_list.txt";
#      my $AR_list_file = &GenUtil::getOutputFile( "AR_loop_list.txt" );
#      while ( $#AR_loop_flist >= 0 ) {
#	printf $AR_list_file "%s\n", pop(@AR_loop_flist);
#      }
#      undef $AR_list_file;
#    }

# jk TESTING
#    $cmd .= " -tight_memory_restrictions -MB_limit_for_rpes 500";
#    $cmd .= " -packer_precompute_only";

    $self->setParameter(multi_chain=>1);
    $self->setParameter(use_input_sc=>1);
    $self->setParameter(decoystats=>1);
    $self->setParameter(interface_ds=>1);
    $self->setParameter(jk_interface=>1);


  } elsif ( $self->{par}->{identifyInterfaceUNS} ) {
    $cmd .= " -identify_interface_UNS";
    $self->setParameter(multi_chain=>1);
    $self->setParameter(use_input_sc=>1);
    $self->setParameter(exp_water=>1);
    $self->setParameter(decoystats=>1);
    $self->setParameter(interface_ds=>1);
    $self->setParameter(jk_interface=>1);

  } elsif ( $self->{par}->{alaScan} ) {
    $cmd .= " -fast_ala_scan scan.out";

  } elsif ( $self->{par}->{desock} ) {
    $cmd .= " -desock";

  } elsif ( $self->{par}->{desrub} ) {
    $cmd .= " -desrub";

  } elsif ( $self->{par}->{des_int} ) {
    $cmd .= " -design_inter";
    printf $tempInfiles "interface_residues\n";

  } else {
    $cmd .= " -fixbb";
  }

#  if ( ! $self->{par}->{incorporate_hotspot} ) {
#    $cmd .= " -ndruns " . $self->{par}->{nstruct};
#  }

  if ( defined $self->{par}->{sidechain_entropy_penalty} ) {
    $cmd .= " -sidechain_entropy_penalty " . $self->{par}->{sidechain_entropy_penalty};
  }

  printf $tempInfiles "redesign\n";

  if ( $self->{par}->{fix_target_seq} > 0 ) {
    $cmd .= " -fix_target_seq " . $self->{par}->{fix_target_seq};
    $cmd .= " -water_on_fixed_chain " if ( $self->{par}->{water_on_fixed_chain} );
  }

#  $cmd .= " -rot_opt";

  my $filelist;
  my $earr;
  my $num_processes;

  ($filelist,$earr,$num_processes) =
      $self->_setupOutfiles($startArr,1);
  $cmd .= " -l " . $filelist;
  printf $tempInfiles "%s\n", $filelist;

  $cmd .= $self->_rosettaSetup($tempInfiles);

  undef $tempInfiles;

  if ( !$self->need_Rosetta($earr) ) {
    &GenUtil::log("no need to call Rosetta in design mode");
    return $self->_processRosettaOutput();
  }

  &GenUtil::log("calling Rosetta in design mode");
  if ( $self->{par}->{cmdLog} ) {
    &GenUtil::log( "Rosetta cmd: " . $cmd );
  }

  my $high_mem_requirements = 0;
  $high_mem_requirements = 1 if ( $self->{par}->{incorporate_hotspot} );

  if ( &UWutil::reportClusterCpus() > 0 ) {
# Note: line below is to force running on whole cluster...
#    undef $num_processes if ( $self->{par}->{incorporate_hotspot} );
    &UWutil::runCmdCluster($cmd, $num_processes, $self->{par}->{outdir}, $self->{par}->{stdout},
			   $self->{par}->{stderr}, 1, $high_mem_requirements);
    sleep(10);
    &UWutil::waitForClusterCompletion();
  } else {
    &UWutil::runCmd($cmd, $self->{par}->{stdout}, $self->{par}->{stderr});
  }

  return $self->_processRosettaOutput();
}

## method: alaScan([$startArr],[parameters])
## Determine stability difference of an Ala at each position

sub alaScan {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $self->{par}->{alaScan} = 1;

  return $self->design($startArr);
}

## method: identifyUNS([$startArr],[parameters])
## Call Rosetta to identify unsatisfied buried Hbond donors/acceptors

sub identifyUNS {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $self->{par}->{fixInterfaceUNS}=0;
  $self->{par}->{identifyInterfaceUNS}=0;
  $self->{par}->{incorporate_hotspot}=0;
  $self->{par}->{alaScan} = 0;
  $self->{par}->{desock} = 0;
  $self->{par}->{desrub} = 0;
  $self->{par}->{des_int} = 0;

  $self->{par}->{decoystats}=1;
  $self->{par}->{exp_water}=1;

  $self->setAllDesType(undef,"NATRO");

  return $self->design($startArr);
}

## method: fixInterfaceUNS([$startArr],[parameters])
## Call Rosetta to fix interface unsatisfied buried Hbond donors/acceptors

sub fixInterfaceUNS {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $self->{par}->{fixInterfaceUNS} = 1;

  return $self->design($startArr);
}

## method: identifyInterfaceUNS([$startArr],[parameters])
## Call Rosetta to identify interface unsatisfied buried Hbond donors/acceptors

sub identifyInterfaceUNS {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $self->{par}->{identifyInterfaceUNS} = 1;

  return $self->design($startArr);
}

## method: incorporateHotspot([$startArr],[parameters])
## Design interfaces which contain hotspots

sub incorporateHotspot {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $self->{par}->{incorporate_hotspot} = 1;

  return $self->design($startArr);
}


## method: patchdock([parameters])
## Run Patchdock
## Applies only to the current structure,
## CANNOT take a pointer to an array of RosettaMolecules.
## Returns a pointer to an array of RosettaMolecules.

sub patchdock {
  my $self = shift;
  $self->_getpar(@_);

  my $tempInfiles = &GenUtil::getOutputFile( $self->{par}->{tempInfiles} );
  printf $tempInfiles "pd_cmd.out\n";
  printf $tempInfiles "pd_cmd.err\n";

  &GenUtil::log("Writing patchdock input PDBs");

  # Store chains in separate files
  my @chainlist=();
  foreach my $c ( @{ $self->activeChains() } ) {
    push(@chainlist,$c->{id});
  }

  die "Too many chains for patchdock" if ( $#chainlist > 1 );
  die "Too few chains for patchdock" if ( $#chainlist < 1 );

  my $param_cmd = $PD_ParamExec;
  my $ms_cmd = $PD_MS_Exec;
  my $chain_num=1;
  foreach my $cid ( @chainlist ) {
    $self->selectChain($cid);
    my $Smol=$self->clone(1);
    my $fname = $cid.".pd_start.pdb";
    $Smol->writePDB($fname, "coor", compress => 0);
    $param_cmd .= " " . $fname;
    $ms_cmd .= " " . $fname;
    printf $tempInfiles "%s\n", $fname;
    printf $tempInfiles "%s.ms\n", $fname;
    # Apply site constraints for this chain
    if ( defined $self->{par}->{dockConsFname} ) {
      my $siteFname = "site".$chain_num.".txt";
      $self->setupPatchdockCst($siteFname,$cid);
      printf $tempInfiles "%s\n", $siteFname;
    }
    $chain_num++;
  }
  $self->selectChain();
  printf $tempInfiles "params.txt\n";
  &GenUtil::log("Writing patchdock param file");
  &UWutil::runCmd($param_cmd, "pd_cmd.out", "pd_cmd.err");
  &GenUtil::log("Writing patchdock ms files");
  &UWutil::runCmd($ms_cmd, "pd_cmd.out", "pd_cmd.err");

  printf $tempInfiles "patch_dock.out\n";
  printf $tempInfiles "patch_dock.log\n";
  my $pd_cmd = $PD_Exec . " params.txt patch_dock.out";
  &GenUtil::log("Calling patchdock");
  &UWutil::runCmd($pd_cmd, "pd_cmd.out", "pd_cmd.err");

  # Build all the output PDBs
  &GenUtil::log("Building patchdock output PDBs");
  # Note: no error message is produced if we ask for more output than are available
  my $num_output_pdbs = 999999;
  $num_output_pdbs = 500;
  my $trans_cmd = $PD_Trans_Exec . " patch_dock.out 1 " . $num_output_pdbs;
  &UWutil::runCmd($trans_cmd, "pd_cmd.out", "pd_cmd.err");

  # Move the output PDBs to the output directory, list them so they'll be found upon processing output
  my $ename = $self->{par}->{outputPDBs};
  printf $tempInfiles "%s\n",$ename;
  my $elist = &GenUtil::getOutputFile($ename);
  &GenUtil::log("Renaming patchdock output PDBs");
  my $pdbnum=1;
  my $exists=1;
  while ( $exists ) {
    my $origname = "patch_dock.out.".$pdbnum.".pdb";
    if ( &GenUtil::exists($origname) ) {
      my $newname=$self->{par}->{outdir} . $pdbnum . ".pdb";
      &GenUtil::rename($origname,$newname);
      printf $elist "%s %s %s\n", $newname, $newname, "start.pdb";
      # Run patchdock output through Molecule, to remove END (and generally clean up)
      my $mol=Molecule->new();
      $mol->readPDB($newname);
      $mol->writePDB($newname);
    } else {
      $exists=0;
    }
    $pdbnum++;
  }

  undef $tempInfiles;

  return $self->_processRosettaOutput();
}

## method: dock([$startArr],[parameters])
## Run Rosetta centroid-mode docking mode
## Applies only to the current structure,
## unless a pointer to an array of
## RosettaMolecules is passed. In this case,
## the current is used as a template and designs
## are carried out for the incoming RosettaMolecules.

sub dock {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $startArr = $self->_setupFromStartArr($startArr);

  my $tempInfiles = &GenUtil::getOutputFile( $self->{par}->{tempInfiles} );

  my $cmd = $rosettaExec;
  $cmd .= " " . $self->{par}->{series} . " " . $self->{par}->{protein} . " A";

  $cmd .= " -dock";
  $cmd .= " -use_score12";

  my $known_type=0;
  if ( $self->{par}->{dockMin} ) {
    $cmd .= " -dock_min";
    $known_type=1;
  }
  if ( $self->{par}->{dockMcm} ) {
    $cmd .= " -dock_mcm";
    $known_type=1;
  }
  if ( defined $self->{par}->{dockPert} ) {
    $cmd .= " -dock_pert " . $self->{par}->{dockPert};
    $known_type=1;
#    $cmd .= " -dock_generate_diversity" if ( $self->{par}->{dockGenerateDiversity} );
  }
  if ( $self->{par}->{dockRand} ) {
    $cmd .= " -randomize1 -randomize2 -I_sc_filter 999";
    $known_type=1;
  }
  if ( ! $known_type ) {
    die "Cannot determine what type of docking to perform";
  }

  if ( ( $self->{par}->{dockMin} ) || ( $self->{par}->{dockMcm} ) ) {
    $cmd .= " -dockFA -scorefxn 12";
    if ( $self->{par}->{dockRTmin} ) {
      $cmd .= " -dock_rtmin";
    }
  }
  if ( ( ! $self->{par}->{dockRand} ) && ( ! $self->{par}->{dockPert} ) ) {
    $cmd .= " -docking_local_refine";
  }
  if ( ! $self->{par}->{dockRand} ) {
    $cmd .= " -no_filters";
  }
  if ( ( ! $self->{par}->{dockMin} ) && ( ! $self->{par}->{dockMcm} ) ) {
    $cmd .= " -norepack1 -norepack2";
  }

  if ( defined $self->{par}->{dockConsFname} ) {
    my $cstFname = $self->{par}->{protein} . ".cst";
    $self->setupDockingCst($cstFname);
    printf $tempInfiles "%s\n", $cstFname;
  }

  $cmd .= " -fa_input";
  $cmd .= " -fa_output";

  $cmd .= " -atom_vdw_set highres";
  #  $cmd.=" -atom_vdw_set hybrid";

  my $filelist;
  my $earr;
  my $num_processes;

  ($filelist,$earr,$num_processes) =
      $self->_setupOutfiles($startArr,0);
  $cmd .= " -l " . $filelist;
  printf $tempInfiles "%s\n", $filelist;

  $cmd .= $self->_rosettaSetup($tempInfiles);

  undef $tempInfiles;

  if ( !$self->need_Rosetta($earr) ) {
    &GenUtil::log("no need to call Rosetta in docking mode");
    return $self->_processRosettaOutput();
  }

  &GenUtil::log("calling Rosetta in docking mode");
  if ( $self->{par}->{cmdLog} ) {
    &GenUtil::log( "Rosetta cmd: " . $cmd );
  }

  my $high_mem_requirements = 0;
  $high_mem_requirements = 1 if ( $self->{par}->{dockMin} || $self->{par}->{dockMcm} );

  if ( &UWutil::reportClusterCpus() > 0 ) {
    &UWutil::runCmdCluster($cmd, $num_processes, $self->{par}->{outdir}, $self->{par}->{stdout},
			   $self->{par}->{stderr}, 1, $high_mem_requirements);
    sleep(10);
    &UWutil::waitForClusterCompletion();
  }
  else {
    &UWutil::runCmd($cmd, $self->{par}->{stdout}, $self->{par}->{stderr});
  }

  return $self->_processRosettaOutput();
}

## method: condenseDockingOutput()
## Bring filenames closer together, so that Rosetta won't skip those it believes are done
## (because Rosetta checks for the existance of the last output file for a given input file)

sub condenseDockingOutput {
  my $self = shift;

  &GenUtil::log("Condensing docking output");

  # Remove all "in_progress" PDBs
  my @pdbs_in_progress = glob($self->{par}->{outdir}."*pdb.in_progress");
  foreach my $progress_fname ( @pdbs_in_progress ) {
    &GenUtil::remove($progress_fname);
  }
  my $dirnum=0;
  while ( -d $dirnum ) {
    my @pdbs_in_progress = glob($self->{par}->{outdir}.$dirnum."/*pdb.in_progress");
    foreach my $progress_fname ( @pdbs_in_progress ) {
      &GenUtil::remove($progress_fname);
    }
  }

  # Collect a list of the expected output PDBs
  my $ename = $self->{par}->{outdir}.$self->{par}->{outputPDBs};
  return 0 if ( ! &GenUtil::exists($ename) );
  my $elist = &GenUtil::getInputFile($ename);
  my @output_names = ();
  while (<$elist>) {
    my $inline = $_;
    chomp($inline);
    my $f = {};
    my @p = split( ' ', $inline );
    push( @output_names, $p[1] );
  }
  undef $elist;
  &GenUtil::log("Expecting ".(($#output_names)+1)." output PDBs");

  # Adjust the numbering to be sequential starting from 1
  my $newnum=0;
  for (my $orignum=0; $orignum <= $#output_names; $orignum++) {
    my $origname=$output_names[$orignum];
    if ( &GenUtil::exists($origname) ) {
      if ( $orignum != $newnum ) {
	my $newname=$output_names[$newnum];
	&GenUtil::rename($origname,$newname);
      }
      $newnum++;
    }
  }
  &GenUtil::log("Condensed names for ".$newnum." output PDBs");

  return $newnum;
}

## method: dockGenerateFunnel([$startArr],[parameters])
## Generate data for docking funnels

sub dockGenerateFunnel {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $self->resetDockFlags();
  $self->{par}->{dockMcm} = 1;
  $self->{par}->{dockPert} = "3 8 8";
  $self->{par}->{dockGenerateDiversity} = 0;
  $self->{par}->{ex_rot} = 2;
  $self->{par}->{dockRTmin} = 1;
  $self->{par}->{interface_ds} = 1;

  return $self->dock($startArr);
}

## method: resetDockFlags()
## Reset docking flags, between calls to "dock" subroutine

sub resetDockFlags {
  my $self = shift;

  $self->{par}->{dockMin} = 0;
  $self->{par}->{dockMcm} = 0;
  $self->{par}->{dockPert} = undef;
  $self->{par}->{dockGenerateDiversity} = 0;
  $self->{par}->{dockRand} = 0;
  $self->{par}->{dockRTmin} = 0;
  $self->{par}->{dockConsFname} = undef;

  return;
}

## method: minimize([$startArr],[parameters])
## Determine stability difference of an Ala at each position

sub minimize {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $self->{par}->{minimize} = 1;

  return $self->relax($startArr);
}

## method: relax([$startArr],[parameters])
## Run Rosetta relax mode
## Relax is for only the current structure,
## unless a pointer to an array of
## RosettaMolecules is passed. In this case,
## the current is used as a template and designs
## are carried out for the incoming RosettaMolecules.

sub relax {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $startArr = $self->_setupFromStartArr($startArr);

  my $tempInfiles = &GenUtil::getOutputFile( $self->{par}->{tempInfiles} );

  my $cmd = $rosettaExec;
  $cmd .= " " . $self->{par}->{series} . " " . $self->{par}->{protein} . " A";

  $cmd .= " -relax -farlx";
  $cmd .= " -minimize" if ($self->{par}->{minimize});
  $cmd .= " -fa_input -fa_output";

  my $filelist;
  my $earr;
  my $num_processes;

  ($filelist,$earr,$num_processes) =
      $self->_setupOutfiles($startArr,0);
  $cmd .= " -l " . $filelist;
  printf $tempInfiles "%s\n", $filelist;

  $cmd .= $self->_rosettaSetup($tempInfiles);

  undef $tempInfiles;

  if ( !$self->need_Rosetta($earr) ) {
    &GenUtil::log("no need to call Rosetta in relax mode");
    return $self->_processRosettaOutput();
  }

  &GenUtil::log("calling Rosetta in relax mode");
  if ( $self->{par}->{cmdLog} ) {
    &GenUtil::log( "Rosetta cmd: " . $cmd );
  }

  if ( &UWutil::reportClusterCpus() > 0 ) {
    &UWutil::runCmdCluster($cmd, $num_processes, $self->{par}->{outdir}, $self->{par}->{stdout}, $self->{par}->{stderr}, 1);
    sleep(10);
    &UWutil::waitForClusterCompletion();
  }
  else {
    &UWutil::runCmd($cmd, $self->{par}->{stdout}, $self->{par}->{stderr});
  }

  return $self->_processRosettaOutput();
}

## function: clean_infile_list($inFname)
## Go through infiles.list, remove any ".gz"

sub clean_infile_list {
  my $inFname = shift;

  my $inF = &GenUtil::getInputFile($inFname);
  my $outF = &GenUtil::getOutputFile($inFname.".tmp");

  while (<$inF>) {
    my $inline=$_;
    chomp($inline);
    ( $inline = $inline ) =~ s/.gz$//gi;
    printf $outF "%s\n",$inline;
  }
  undef $inF;
  undef $outF;

  &GenUtil::rename($inFname.".tmp",$inFname);

  return;

}

## method: abinitFrag([$startArr],[parameters])
## Call Rosetta ab initio mode for fragment insertions
## Acts on only the current structure,
## unless a pointer to an array of
## RosettaMolecules is passed. In this case,
## the current is used as a template and modeling
## is carried out for the incoming RosettaMolecules

sub abinitFrag {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $startArr = $self->_setupFromStartArr($startArr);

  my $tempInfiles = &GenUtil::getOutputFile( $self->{par}->{tempInfiles} );

  my $cmd = $rosettaExec;
  $cmd .= " " . $self->{par}->{series} . " " . $self->{par}->{protein} . " A";

  $cmd .= " -regions";
  if ( ( !defined $self->{par}->{no_seq} ) || ( $self->{par}->{no_seq} == 1 ) )
  {
    $cmd .= " -no_seq";
    $cmd .= " -scorefxn 70";
  }

  my $filelist = "infiles.list";
  $cmd .= " -l " . $filelist;
  my $flist = &GenUtil::getOutputFile($filelist);
  printf $tempInfiles "%s\n", $filelist;
  my $ename = $self->{par}->{outdir} . $self->{par}->{outputPDBs};
  my $elist = &GenUtil::getOutputFile($ename);
  my $earr = ();
  my $fnum = 1;
  my $numOutfiles = 0;

  foreach my $startfile ( @{$startArr} ) {
    my $outdir = ($fnum % 100) . "/";
    my $fulldir = $self->{par}->{outdir}.$outdir;
    &GenUtil::makeDir($fulldir) if ( !-d $fulldir );
    my $prefix = "temp" . $fnum;
    my $fname = $prefix . ".pdb";
    printf $flist "%s %s\n", $startfile->{filename}, $outdir.$fname;
    my $efile;
    for ( my $i = 1 ; $i <= $self->{par}->{nstruct} ; $i++ ) {

      # Note: abinitio mode output is different from loop mode!
      $efile = $fulldir . $self->{par}->{series};
      $efile .= $self->{par}->{protein} . &GenUtil::zPad( $i, 4 ) . ".pdb";
      my $cfile = $fulldir . $self->{par}->{series};
      $cfile .=
        $self->{par}->{protein} . &GenUtil::zPad( $i, 4 ) . ".clean.pdb";
      printf $elist "%s %s %s\n", $cfile, $efile, $startfile->{filename};
      $numOutfiles++;
    }
    push( @{$earr}, $efile );
    my $regionFname = $prefix . ".regions";
    $self->writeRegionFile($regionFname);
    printf $tempInfiles "%s\n", $regionFname;
    $fnum++;
  }
  undef $flist;
  undef $elist;

  my $aa3 = "aa" . $self->{par}->{protein} . "A03_05.200_v1_3";
  my $aa9 = "aa" . $self->{par}->{protein} . "A09_05.200_v1_3";
  $self->_makeFragLib( $self->{par}->{outdir} . "makeFrag.out" )
    if ( ( !-e $aa3 ) || ( !-e $aa9 ) );

  &clean_infile_list($filelist);

  $cmd .= $self->_rosettaSetup($tempInfiles);

  undef $tempInfiles;

  if ( !$self->need_Rosetta($earr) ) {
    &GenUtil::log("no need to call Rosetta in ab initio mode");
    return $self->_processRosettaOutput();
  }

  &GenUtil::log("calling Rosetta in ab initio mode");
  if ( $self->{par}->{cmdLog} ) {
    &GenUtil::log( "Rosetta cmd: " . $cmd );
  }

  if ( &UWutil::reportClusterCpus() > 0 ) {
    &UWutil::runCmdCluster($cmd, $numOutfiles, $self->{par}->{outdir}, $self->{par}->{stdout}, $self->{par}->{stderr}, 1);
    sleep(10);
    &UWutil::waitForClusterCompletion();
  }
  else {
    &UWutil::runCmd($cmd, $self->{par}->{stdout}, $self->{par}->{stderr});
  }

  return $self->_processRosettaOutput();
}

## method: loopCenRlx([$startArr],[parameters])
## Call Rosetta (centroid) relax mode

sub loopCenRlx {
  my $self = shift;

  $self->setParameter( loop_fix_natsc => 0 );
  $self->setParameter( loop_fa_input => 0 );
  $self->setParameter( loop_fa_output => 0 );
  $self->setParameter( loop_fa_refine => 0 );
  $self->setParameter( cen_rlx => 1 );
  $self->setParameter( FA_rlx => 0 );
  $self->setParameter( no_seq => 1 );

  return $self->loops(@_);

}

## method: loopFaRlx([$startArr],[parameters])
## Call Rosetta (full-atom) relax mode

sub loopFaRlx {
  my $self = shift;

  $self->setParameter( loop_fix_natsc => 1 );
  $self->setParameter( loop_fa_input => 1 );
  $self->setParameter( loop_fa_output => 1 );
  $self->setParameter( loop_fa_refine => 0 );
  $self->setParameter( cen_rlx => 0 );
  $self->setParameter( FA_rlx => 1 );
  $self->setParameter( no_seq => 0 );

  return $self->loops(@_);

}

## method: loopFrag([$startArr],[parameters])
## Call Rosetta loop mode for fragment insertions

sub loopFrag {
  my $self = shift;

  $self->setParameter( loop_fix_natsc => 0 );
  $self->setParameter( loop_fa_input => 0 );
  $self->setParameter( loop_fa_output => 0 );
  $self->setParameter( loop_fa_refine => 0 );
  $self->setParameter( cen_rlx => 0 );
  $self->setParameter( FA_rlx => 0 );
  $self->setParameter( no_seq => 0 );

  return $self->loops(@_);

}

## method: loopRefine([$startArr],[parameters])
## Call Rosetta loop mode for all-atom refinement

sub loopRefine {
  my $self = shift;

  $self->setParameter( loop_fix_natsc => 1 );
  $self->setParameter( loop_fa_input => 1 );
  $self->setParameter( loop_fa_output => 1 );
  $self->setParameter( loop_fa_refine => 1 );
  $self->setParameter( cen_rlx => 0 );
  $self->setParameter( FA_rlx => 0 );
  $self->setParameter( no_seq => 0 );

  return $self->loops(@_);

}

## method: loops([$startArr],[parameters])
## Run Rosetta loop mode
## Loop mode is for only the current structure,
## unless a pointer to an array of
## RosettaMolecules is passed. In this case,
## the current is used as a template and loop mode
## is carried out for the incoming RosettaMolecules

sub loops {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $startArr = $self->_setupFromStartArr($startArr);

  my $tempInfiles = &GenUtil::getOutputFile( $self->{par}->{tempInfiles} );

  my $cmd = $rosettaExec;
  $cmd .= " " . $self->{par}->{series} . " " . $self->{par}->{protein} . " A";

  if ( ( $self->{par}->{cen_rlx} ) || ( $self->{par}->{FA_rlx} ) ) {
    $cmd .= " -relax";
    $cmd .= " -cenrlx" if ( $self->{par}->{cen_rlx} );
    $cmd .= " -farlx" if ( $self->{par}->{FA_rlx} );
  }
  else {
    $cmd .= " -fa_refine" if ( $self->{par}->{loop_fa_refine} );
  }

  $cmd .= " -loops";

  if ( ( $self->{par}->{FA_rlx} ) || ( $self->{par}->{loop_fa_refine} ) ) {
    $cmd .= " -minimize";
    $cmd .= " -scorefxn 12";
  }

  $cmd .= " -no_seq -scorefxn 70"
    if ( ( $self->{par}->{cen_rlx} ) && ( $self->{par}->{no_seq} ) );

  my $filelist;
  my $earr;
  my $num_processes;

  ($filelist,$earr,$num_processes) =
      $self->_setupOutfiles($startArr,0);
  $cmd .= " -l " . $filelist;
  printf $tempInfiles "%s\n", $filelist;

  my $aa3 = "aa" . $self->{par}->{protein} . "A03_05.200_v1_3";
  my $aa9 = "aa" . $self->{par}->{protein} . "A09_05.200_v1_3";
  $self->_makeFragLib( $self->{par}->{outdir} . "makeFrag.out" )
    if ( ( !-e $aa3 ) || ( !-e $aa9 ) );
  $cmd .= " -fa_input"  if ( $self->{par}->{loop_fa_input} );
  $cmd .= " -fa_output" if ( $self->{par}->{loop_fa_output} );
  $cmd .= " -fix_natsc" if ( $self->{par}->{loop_fix_natsc} );

  &clean_infile_list($filelist);

  $cmd .= $self->_rosettaSetup($tempInfiles);

  undef $tempInfiles;

  if ( !$self->need_Rosetta($earr) ) {
    &GenUtil::log("no need to call Rosetta in loop mode");
    return $self->_processRosettaOutput();
  }

  &GenUtil::log("calling Rosetta in loop mode");
  if ( $self->{par}->{cmdLog} ) {
    &GenUtil::log( "Rosetta cmd: " . $cmd );
  }

  if ( &UWutil::reportClusterCpus() > 0 ) {
    &UWutil::runCmdCluster($cmd, $num_processes, $self->{par}->{outdir}, $self->{par}->{stdout}, $self->{par}->{stderr}, 1);
    sleep(10);
    &UWutil::waitForClusterCompletion();
  }
  else {
    &UWutil::runCmd($cmd, $self->{par}->{stdout}, $self->{par}->{stderr});
  }

  return $self->_processRosettaOutput();
}

## method: charmm_mini([$startArr],[parameters])
## Call charmm for minimization
## Operates in a manner such that input and output are
## as similar as possible to those of the Rosetta
## subroutines above

sub charmm_mini {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  my $chaMiniScript = &GenUtil::findExecutable("charmm_mini.pl");
  die "cannot find the charmm minimization script"
    if ( !defined $chaMiniScript );

  $startArr = $self->_setupFromStartArr($startArr);

  &GenUtil::makeDir( $self->{par}->{outdir} ) if ( !-d $self->{par}->{outdir} );

  my $tempInfiles = &GenUtil::getOutputFile( $self->{par}->{tempInfiles} );

  my $filelist = "infiles.list";
  printf $tempInfiles "%s\n", $filelist;

  my $flist = &GenUtil::getOutputFile($filelist);
  my $ename = $self->{par}->{outdir} . $self->{par}->{outputPDBs};
  my $elist = &GenUtil::getOutputFile($ename);
  my $earr = ();
  my $fnum = 1;
  my $numOutfiles = 0;

  foreach my $startfile ( @{$startArr} ) {
    my $outdir = ($fnum % 100) . "/";
    my $fulldir = $self->{par}->{outdir}.$outdir;
    &GenUtil::makeDir($fulldir) if ( !-d $fulldir );
    my $efile = $fulldir . $self->{par}->{series} . "temp" . $fnum . ".mini.pdb";
    printf $flist "%s %s\n", $startfile->{filename}, $outdir.$efile;
    my $cfile = $fulldir . $self->{par}->{series} . "temp";
    $cfile .= $fnum . ".mini.clean.pdb";
    printf $elist "%s %s %s\n", $cfile, $efile, $startfile->{filename};
    push( @{$earr}, $efile );
    $numOutfiles++;
    $fnum++;
  }

  undef $flist;
  undef $elist;
  my $cmd = $chaMiniScript;
  $cmd .= " -interface" if ( $self->{par}->{chamin_int} );
  $cmd .= " " . $filelist;
  printf $tempInfiles "%s\n", &UWutil::reportJobFilename()
    if ( &UWutil::reportClusterCpus() > 0 );
  undef $tempInfiles;

  if ( !$self->need_Rosetta($earr) ) {
    &GenUtil::log("no need to call Charmm for minimization");
    return $self->_processRosettaOutput();
  }

  &GenUtil::log("calling Charmm for minimization");

  if ( &UWutil::reportClusterCpus() > 0 ) {
    &UWutil::runCmdCluster($cmd, $numOutfiles, $self->{par}->{outdir}, $self->{par}->{stdout}, $self->{par}->{stderr}, 1);
    sleep(10);
    &UWutil::waitForClusterCompletion();
  }
  else {
    &UWutil::runCmd($cmd, $self->{par}->{stdout}, $self->{par}->{stderr});
  }

  return $self->_processRosettaOutput();
}

## method: writePathFile([parameters])
## generates a pathfile for rosetta and creates output directories if needed

sub writePathFile {
  my $self = shift;

  $self->_getpar(@_);

  my $outdir = $self->{par}->{outdir};
  if ( ( $outdir !~ /^\./ ) && ( $outdir !~ /^\// ) ) {
    $outdir = "./" . $outdir;
  }

  my $file = &GenUtil::getOutputFile( $self->{par}->{pathfile} );
  printf $file "Rosetta Input/Output Paths (order essential)\n";
  printf $file
    "path is first '/', './',or  '../' to next whitespace, must end with '/'\n";
  printf $file "INPUT PATHS:\n";
  printf $file "pdb1	                        ./\n";
  printf $file "pdb2	                        ./\n";
  printf $file "alternate data files            /work/johnk/rosetta_database/\n";
  printf $file "fragments                       ./\n";
  printf $file "structure dssp,ssa (dat,jones)  ./\n";
  printf $file "sequence fasta,dat,jones        ./\n";
  printf $file "constraints                     ./\n";
  printf $file "starting structure              ./\n";
  printf $file "data files                      %s\n", $rosettaData;
  printf $file "OUTPUT PATHS:\n";
  printf $file "movie                           %s\n", $outdir;
  printf $file "pdb path                        %s\n", $outdir;
  printf $file "score                           %s\n", $outdir;
  printf $file "status                          %s\n", $outdir;
  printf $file "user                            %s\n", $outdir;
  printf $file "FRAGMENTS: (use '*****' in place of pdb name and chain)\n";
  printf $file
    "2                                      number of valid fragment files\n";
  printf $file "3                                      frag file 1 size\n";
  printf $file "aa*****03_05.200_v1_3                               name\n";
  printf $file "9                                      frag file 2 size\n";
  printf $file "aa*****09_05.200_v1_3                               name\n";
  printf $file
"-------------------------------------------------------------------------\n";
  printf $file "CVS information:\n";
  printf $file "\$Revision: 1.10 \$\n";
  printf $file "\$Date: 2002/07/16 16:47:31 \$\n";
  printf $file "\$Author: rohl \$\n";
  printf $file
"-------------------------------------------------------------------------\n";
  undef $file;

  &GenUtil::makeDir($outdir) if ( !-d $outdir );
  return;
}

## method: writeResfile([resfileName])
## Write a resfile

sub writeResfile {
  my $self = shift;
  my $resfileName = shift;

  die "No chains defined for writeResfile"
    if ( !defined $self->activeChains() );

  $resfileName = $self->{par}->{resfile} if ( !defined $resfileName );

  my $resout = &GenUtil::getOutputFile($resfileName);
  printf $resout " start\n";

  my $totres = 1;
  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      my $cid = $c->{id};
      if ( $cid ne "+" ) {
	$cid = " " if ( $cid eq "" );
	for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
	  my $ri = $c->{res}->[$ir];
	  printf $resout " %s %4d %4d %5s\n", $cid, $totres, $ri->{num},
	    $ri->{des};
	  $totres++;
	}
      }
    }
  }
  undef $resout;

  return;
}

## method: readResfile([resfileName])
## Read a resfile, set the {des} fields accordingly

sub readResfile {
  my $self = shift;
  my $resfileName = shift;

  $resfileName = $self->{par}->{resfile} if ( !defined $resfileName );

  my $resin = &GenUtil::getInputFile($resfileName);

  my $got_start=0;
  while (<$resin>) {
    my $inline=$_;
    chomp($inline);
    if ( $got_start) {
      my @p = split( ' ', $inline );
      my $c = $self->getChain( $p[0] );
      die "Cannot find chain ".$p[0]." in readResfile" if ( !defined $c->{res} );
      my $type = $p[3];
      $type .= "  " . $p[4] if ( $#p >= 4 );
      my $r = $self->getResidueInChain( $p[2], $c );
      $r->{des} = $type if ( defined $r );
    } elsif ($inline =~ /^ start/) {
      $got_start=1;
    }
  }
  undef $resin;

  return;
}

## method: restrictSeqToNative()
## Set all {des} fields which are not NATRO to
## NOTAA restricting to native sequence

sub restrictSeqToNative {
  my $self = shift;

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
        my $ri = $c->{res}->[$ir];
	if ( $ri->{des} ne "NATRO" ) {
	    $ri->{des}="NOTAA  ".
		&Sequence::notRes($Sequence::_seqabbrev{$ri->{name}});
	}
      }
    }
  }

  return;
}

## method: setupDockingCst()
## Read either a .cst or .cons file,
## write a .cst (docking site constraint file)
## Note: at the moment, only supports ONE site per chain

sub setupDockingCst {
  my $self = shift;
  my $outFname = shift;

  my $inFname = $self->{par}->{dockConsFname};
  return if ( !defined $inFname );

  $outFname = $self->{par}->{protein} . ".cst"
    if ( !defined $outFname );

  # Read the whole file at once
  my $infile = &GenUtil::getInputFile($inFname);
  my @lines = <$infile>;

  # Figure out whether it's a .cst or a .cons file
  # (it's a .cst if it has chain ID's (including "_"))
  my $cst = 0;
  foreach my $inline (@lines) {
    $cst = 1 if ( $inline =~ /[A-Z]/ );
    $cst = 1 if ( $inline =~ /_/ );
  }

  if ($cst) {

    # If we already have a .cst file, just link to it
    &GenUtil::link( $inFname, $outFname );
  }
  else {

    # Filter val is 1 greater than the sum of all the -ve terms
    # (ie. output decoys must use none of the +ve sites, all of the -ve ones)
    my $filter_val = 1.0;

    my $chash = {};
    if ( $#lines > -1 ) {
      my $numSites = $#lines + 1;
      my $outF = &GenUtil::getOutputFile($outFname);
      printf $outF "%d\n", $numSites;

      foreach my $inline (@lines) {
	chomp($inline);
	my @p = split( ' ', $inline );
	my $reslist = $p[0];
	my $bonus = -100;
	$bonus = $p[1] if ( $#p > 0 );
        printf $outF "%s\n", $bonus;
	$filter_val += $bonus if ($bonus < 0);

	my @resarr = split(/,/, $reslist);
	my $cid;
	foreach my $resnum ( @resarr ) {
	  foreach my $c ( @{ $self->activeChains() } ) {
	    my $r = $self->getResidueInChain( $resnum, $c );
	    $cid = $c->{id} if ( defined $r );
	  }

	  die "Could not find a chain with residue " . $resnum
	    if ( !defined $cid );
	  printf $outF "%d %s\n", $resnum, $cid;
	}

	printf $outF "0 %s\n", $cid;
      }

      printf $outF "filter %f\n", $filter_val;
      undef $outF;

    }
  }

  return;
}

## method: setupPatchdockCst()
## Read a .cons file, write a site?.txt (Patchdock active site file)

sub setupPatchdockCst {
  my $self = shift;
  my $outFname = shift;
  my $desired_cid = shift;

  die if ( !defined $outFname );
  die if ( !defined $desired_cid );

  my $inFname = $self->{par}->{dockConsFname};
  return if ( !defined $inFname );

  # Read the whole file at once
  my $infile = &GenUtil::getInputFile($inFname);
  my @lines = <$infile>;
  undef $infile;

  # Figure out whether it's a .cst or a .cons file
  # (it's a .cst if it has chain ID's (including "_"))
  my $cst = 0;
  foreach my $inline (@lines) {
    $cst = 1 if ( $inline =~ /[A-Z]/ );
    $cst = 1 if ( $inline =~ /_/ );
  }

  die "cst files are not supported" if ($cst);

  my $chash = {};
  if ( $#lines > -1 ) {
    my $numSites = $#lines + 1;
    my $outF = &GenUtil::getOutputFile($outFname);

    foreach my $inline (@lines) {
      chomp($inline);
      my @p = split( ' ', $inline );
      my $reslist = $p[0];
      my $bonus = -1;
      $bonus = $p[1] if ( $#p > 0 );

      if ( $bonus < 0 ) {
	my @resarr = split(/,/, $reslist);
	my $found_cid;
	foreach my $resnum ( @resarr ) {
	  foreach my $c ( @{ $self->activeChains() } ) {
	    my $r = $self->getResidueInChain( $resnum, $c );
	    $found_cid = $c->{id} if ( defined $r );
	  }

	  if ( $found_cid eq $desired_cid ) {
	    printf $outF "%d %s\n", $resnum, $found_cid;
	  }
	}
      }
    }

    undef $outF;

  }

  return;
}

## method: setRepackNearby([cid])
## Set residues near designed residues to repack,
## return an array of residue numbers to repack.
## If a chain ID is supplied, only residues on
## this chain are repacked

sub setRepackNearby {
  my $self = shift;
  my $cid = shift;

  my $repackDist = $self->{par}->{repackDist};
  my $repackList = ();

  if ( defined $self->activeChains() ) {
    foreach my $c1 ( @{ $self->activeChains() } ) {
      if ( $#{ $c1->{res} } >= 0 ) {
        for ( my $ir = 0 ; $ir <= $#{ $c1->{res} } ; $ir++ ) {
          my $ri = $c1->{res}->[$ir];
          if ( $ri->{des} eq "ALLAA" ) {
            foreach my $c2 ( @{ $self->activeChains() } ) {
              if ( ( ( !defined $cid ) || ( $cid eq $c2->{id} ) )
                && ( $#{ $c2->{res} } >= 0 ) )
              {
                for ( my $jr = 0 ; $jr <= $#{ $c2->{res} } ; $jr++ ) {
                  my $rj = $c2->{res}->[$jr];
                  if ( $rj->{des} eq "NATRO" ) {
                    if ( $self->minDistance( $ri, $rj, 1 ) <= $repackDist ) {
                      $rj->{des} = "NATAA";
                      push( @{$repackList}, $rj->{num} );
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return $repackList;
}

## method: setDesignInterface()
## Set residues near interface to be designed

sub setDesignInterface {
  my $self = shift;
  my $desDist = shift;

  $desDist = $self->{par}->{interfaceDist} if ( !defined $desDist );
  my $desHash = {};

  for ( my $cinx1 = 0 ; $cinx1 < $#{ $self->{chain} } ; $cinx1++ ) {
    my $c1 = $self->{chain}->[$cinx1];
    for ( my $cinx2 = $cinx1 + 1 ; $cinx2 <= $#{ $self->{chain} } ; $cinx2++ ) {
      my $c2 = $self->{chain}->[$cinx2];

      if ( ( $#{ $c1->{res} } >= 0 ) && ( $#{ $c2->{res} } >= 0 ) ) {
        my $cid = $c1->{id};
        for ( my $ir = 0 ; $ir <= $#{ $c1->{res} } ; $ir++ ) {
          my $ri = $c1->{res}->[$ir];
          for ( my $jr = 0 ; $jr <= $#{ $c2->{res} } ; $jr++ ) {
            my $rj = $c2->{res}->[$jr];
            if ( ( $ri->{des} ne "ALLAA" ) || ( $rj->{des} ne "ALLAA" ) ) {
              if ( $self->minDistance( $ri, $rj, 0, 1 ) <= $desDist ) {
                $ri->{des} = "ALLAA";
                $desHash->{ $c1->{id} . $ri->{num} } = 1;
              }
              if ( $self->minDistance( $rj, $ri, 0, 1 ) <= $desDist ) {
                $rj->{des} = "ALLAA";
                $desHash->{ $c2->{id} . $rj->{num} } = 1;
              }
            }
          }
        }
      }
    }
  }

  my $desList = ();
  foreach my $k ( keys %{$desHash} ) {
    push( @{$desList}, $k );
  }

  return $desList;
}

## method: limitDesBuried()
## Restrict current set of designed residues to those which are buried
## Note: B-factor field (aux2) is assumed to contain SASA!!

sub limitDesBuried {
  my $self = shift;
  my $scOnly = shift;

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      foreach my $r ( @{ $c->{res} } ) {
        if ( $r->{des} ne "NATRO" ) {
          if ( $self->residueSASA( $r, $scOnly ) > 0.5 ) {
            $r->{des} = "NATRO";
          }
        }
      }
    }
  }

  return;
}

## method: SASAatomicInterface(refMol)
## return an array of atoms which comprise the interface, defined
## as atoms which are more buried in the current structure than the reference
## and are within 5 A of an atom on a different chain
## Note: assumes the B-factor field contains the atomic SASA

sub SASAatomicInterface {
  my $self = shift;
  my $refMol = shift;
  my $SASAdiffCut = shift;
  my $scOnly = shift;

  $SASAdiffCut = 0.5 if ( !defined $SASAdiffCut );
  $scOnly = 0   if ( !defined $scOnly );

  my $alist = ();

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      my $refC = $refMol->getChain( $c->{id} );
      foreach my $r ( @{ $c->{res} } ) {
        my $start = $r->{start};
        my $end = $r->{end};
        my $refres = $refMol->getResidueInChain( $r->{num}, $c );
        my $refstart = $refres->{start};
        my $refend = $refres->{end};
        for ( my $ai = $start ; $ai <= $end ; $ai++ ) {
          my $atom = $c->{atom}->[$ai];
          if ( ( $atom->{aux1} >= 0.0 )
            && ( !$atom->{hyd} )
            && ( ( !$scOnly ) || ( !$atom->{bb} ) ) )
          {
            my $atomname = $atom->{atomname};
            for ( my $rai = $refstart ; $rai <= $refend ; $rai++ ) {
              if ( $refC->{atom}->[$rai]->{atomname} eq $atomname ) {
		my $sasaA=$atom->{aux2};
	        $sasaA=0.0 if ($sasaA < 0.0);
		my $sasaC=$refC->{atom}->[$rai]->{aux2};
	        $sasaC=0.0 if ($sasaC < 0.0);
                my $SASAdiff = $sasaC-$sasaA;
                if ( $SASAdiff >= $SASAdiffCut ) {
                  push( @{$alist}, $atom )
                    if ( $self->interfaceAtom( $atom, 5 ) );
                }
              }
            }
          }
        }
      }
    }
  }

  return $alist;
}

## method: limitDesBuriedVsRef()
## Restrict current set of designed residues to those which are more buried
## in the current structure relative to a reference structure (eg. unbound state)

sub limitDesBuriedVsRef {
  my $self = shift;
  my $refMol = shift;
  my $SASAdiffCut = shift;
  my $scOnly = shift;

  $SASAdiffCut = 20 if ( !defined $SASAdiffCut );

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      foreach my $r ( @{ $c->{res} } ) {
        if ( $r->{des} ne "NATRO" ) {
          my $refRes = $refMol->getResidueInChain( $r->{num}, $c );
          my $SASAdiff =
            $refMol->residueSASA( $refRes, $scOnly ) -
            $self->residueSASA( $r,        $scOnly );
          $r->{des} = "NATRO" if ( $SASAdiff < $SASAdiffCut );
        }
      }
    }
  }

  return;
}

## method: limitDesInterfaceContacts()
## Restrict current set of designed residues to those which have at least
## reqNum contacts to a different chain within distCut

sub limitDesInterfaceContacts {
  my $self = shift;
  my $distCut = shift;
  my $reqNum = shift;
  my $scOnly = shift;

  $distCut = 5 if ( !defined $distCut );
  my $SqDistCut = $distCut * $distCut;

  $reqNum = 1 if ( !defined $reqNum );
  $scOnly = 0 if ( !defined $scOnly );

  my @ac = @{ $self->activeChains() };

  for ( my $i = 0 ; $i <= $#ac ; $i++ ) {
    my $c1 = $ac[$i];
    for ( my $j = $i + 1 ; $j <= $#ac ; $j++ ) {
      my $c2 = $ac[$j];
      if ( ( $#{ $c1->{res} } >= 0 ) && ( $#{ $c2->{res} } >= 0 ) ) {

        my $ca1 = $c1->{atom};
        my $ca2 = $c2->{atom};
        foreach my $r1 ( @{ $c1->{res} } ) {
          if ( $r1->{des} ne "NATRO" ) {
            my $numC = 0;
            for ( my $ai1 = $r1->{start} ; $ai1 <= $r1->{end} ; $ai1++ ) {
              my $a1 = $ca1->[$ai1];
              if ( ( $a1->{aux1} >= 0.0 )
                && ( !$a1->{hyd} )
                && ( ( !$scOnly ) || ( !$a1->{bb} ) ) )
              {
                my $x1 = $a1->{xcoor};
                my $y1 = $a1->{ycoor};
                my $z1 = $a1->{zcoor};
                foreach my $r2 ( @{ $c2->{res} } ) {
                  for ( my $ai2 = $r2->{start} ; $ai2 <= $r2->{end} ; $ai2++ ) {
                    my $a2 = $ca2->[$ai2];
                    if ( ( $a2->{aux1} >= 0.0 ) && ( !$a2->{hyd} ) ) {
                      my $x2 = $a2->{xcoor};
                      my $y2 = $a2->{ycoor};
                      my $z2 = $a2->{zcoor};
                      my $SqDist = ( $x1 - $x2 ) * ( $x1 - $x2 );
                      $SqDist += ( $y1 - $y2 ) * ( $y1 - $y2 );
                      $SqDist += ( $z1 - $z2 ) * ( $z1 - $z2 );
                      $numC++ if ( $SqDist < $SqDistCut );
                    }
                  }
                }
              }
            }
            $r1->{des} = "NATRO" if ( $numC < $reqNum );
          }
        }

      }
    }
  }

  return;
}

## function: _randID()
## return a random protein ID (a number followed by three letters)

sub _randID {

  # Give the protein a random four-letter code, making
  # sure it doesn't already exist in the current dir
  # (destructor will remove it)
  my $randID;
  my $validID = 0;
  while ( !$validID ) {
#    $randID = int( rand(7) ) + 3;    # get an int between 3 and 9
    # note: reserve those starting with 1 or 2, these are likely real PDBs!
                                     # get three uppercase chars between A-Z
    $randID = 9;    # change of plans - start with 9 every time
    $randID .= chr( int( rand(26) ) + 65 );
    $randID .= chr( int( rand(26) ) + 65 );
    $randID .= chr( int( rand(26) ) + 65 );
    my @f = glob( $randID . "*" );
    my @g = glob( lc($randID) . "*" );
    $validID = 1 if ( ( $#f == -1 ) && ( $#g == -1 ) );
  }

  return $randID;
}

## method: _convertWater()
## Rename *.?W?? waters in Rosetta output to HETATM waters

sub _convertWater {
  my $self = shift;

  my $chainrec;
  my $crec = $self->{chainlookup}->{"+"};
  if ( defined $crec ) {
    $chainrec = $crec;
  }
  else {
    $chainrec = $self->_newChain("+");
  }

  my $nextres = 0;
  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      my $lastinx = $#{ $c->{res} };
      my $lastres = $c->{res}->[$lastinx]->{num};
      $nextres = $lastres if ( $lastres > $nextres );
    }
  }
  $nextres++;

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{atom} } >= 0 ) {
      my $lastatm = $#{ $c->{atom} };
      for ( my $i = 0 ; $i <= $lastatm ; $i++ ) {
        my $a = $c->{atom}->[$i];
        if ( ( $a->{atomname} =~ "W" )
          && ( defined $Sequence::_seqabbrev{ $a->{resname} } ) )
        {

          my $resnum = $a->{resnum};
          foreach my $r ( @{ $c->{res} } ) {
            if ( $r->{num} > $resnum ) {
              $r->{start}--;
              $r->{end}--;
            }
            elsif ( $r->{num} == $resnum ) {
              $r->{end}--;
            }
          }

          $a->{atomname} = "O";
          $a->{chain} = "+";
          $a->{resname} = "HOH";
          $a->{resnum} = $nextres;
          push( @{ $chainrec->{atom} }, $a );

          splice( @{ $c->{atom} }, $i, 1 );

          my $resrec = {};
          $resrec->{name} = "HOH";
          $resrec->{num} = $nextres;
          $resrec->{chain} = "+";
          $resrec->{start} = $#{ $chainrec->{atom} };
          $resrec->{end} = $resrec->{start};
          $resrec->{valid} = 1;
          push( @{ $chainrec->{res} }, $resrec );

          $nextres++;
          $i--;
          $lastatm--;

        }
      }
    }
  }

  undef $self->{par}->{currFile};

  return;
}

## method: unsSet = collectUnsSet()
## build an unsSet from the current molecule

sub collectUnsSet {
  my $self = shift;

  my $unsSet = {};

  if ( ( defined $self->{uns} ) && ( $#{ $self->{uns} } >= 0 ) ) {
    foreach my $uns ( @{ $self->{uns} } ) {
      $unsSet->{ $uns->{resnum} }->{ $uns->{restype} }->{ $uns->{atomname} } =
        1;
    }
  }

  return $unsSet;
}

## method: unsArr = listUns()
## return the UNS array

sub listUns {
  my $self = shift;

  my $uArr=();
  foreach my $uns ( @{$self->{uns}} ) {
      push(@{$uArr},$uns);
  }

  return $uArr;
}

## method: groupUnsArr = listGroupUns()
## return the group UNS array

sub listGroupUns {
  my $self = shift;

  my $gArr=();
  foreach my $gu ( @{$self->{gu}} ) {
      push(@{$gArr},$gu);
  }

  return $gArr;
}

## method: writeUns(fname)
## write UNS in the current molecule to a file

sub writeUns {
  my $self = shift;
  my $fname = shift;

  $fname = $self->{par}->{protein} . ".uns" if ( !defined $fname );
  my $outF = &GenUtil::getOutputFile($fname);

  if ( ( defined $self->{uns} ) && ( $#{ $self->{uns} } >= 0 ) ) {
    foreach my $uns ( @{ $self->{uns} } ) {
      printf $outF "%d %s %s\n", $uns->{resnum}, $uns->{restype},
        $uns->{atomname};
    }
  }

  return;
}

## function: unsSet = readUns(fname)
## read UNS from a file, return an unsSet

sub readUns {
  my $fname = shift;

  my $inF = &GenUtil::getInputFile($fname);

  my $unsSet = {};

  while (<$inF>) {
    chomp;
    my $inline = $_;
    my @p = split( ' ', $inline );
    if ( $#p != 2 ) {
      printf STDERR "Unexpected line in UNS file: %s\n", $inline;
      exit(1);
    }

    my $resnum = $p[0];
    my $restype = $p[1];
    my $atomname = $p[2];
    $unsSet->{$resnum}->{$restype}->{$atomname} = 1;
  }

  return $unsSet;
}

## method: countFunctionalContacts()
## Count the number of inter-chain contacts involving
## a subset of residues

sub countFunctionalContacts {
  my $self = shift;
  my $funcList = shift;

  my $desDist = $self->{par}->{interfaceDist};
  my $numFunc = 0;

  return 0 if ( !defined $self->activeChains() );

  foreach my $fnum ( @{$funcList} ) {
    foreach my $fc ( @{ $self->activeChains() } ) {
      my $fres = $self->getResidueInChain( $fnum, $fc );
      if ( defined $fres ) {
        my $found = 0;
        foreach my $c ( @{ $self->activeChains() } ) {
          if ( ( $#{ $c->{res} } >= 0 ) && ( $c->{id} ne $fc->{id} ) ) {
            for (
              my $jr = 0 ;
              ( ( $jr <= $#{ $c->{res} } ) && ( !$found ) ) ;
              $jr++
              )
            {
              $found = 1
                if ( $self->minDistance( $fres, $c->{res}->[$jr], 0, 1 ) <=
                $desDist );
            }
          }
        }
        $numFunc++ if ($found);
      }
    }
  }

  return $numFunc;
}

## method: countInterfaceDesignContacts([distCut][,chainA][,chainB][,CACBonly][,SConly])
## Count the number of atomic contacts across an interface
## Operates on all pairs of chains if no chain IDs are specified
## Operates on all partners for given chain if one chain ID is specified
## Operates only on one pair if two chain IDs are specified
## Returns an array of hashes with {AAcontacts}, {cidA}, {cidB}

sub countInterfaceDesignContacts {
  my $self = shift;
  my $distCut = shift;
  my $cidA = shift;
  my $cidB = shift;
  my $CACBonly = shift;
  my $SConly = shift;

  $distCut = 5 if ( !defined $distCut );
  my $SqDistCut = $distCut * $distCut;

  $CACBonly = 0 if ( !defined $CACBonly );
  $SConly = 0 if ( !defined $SConly );

  my @ac = @{ $self->activeChains() };
  my $AAlist = ();

  for ( my $i = 0 ; $i <= $#ac ; $i++ ) {
    my $c1 = $ac[$i];
    for ( my $j = $i + 1 ; $j <= $#ac ; $j++ ) {
      my $c2 = $ac[$j];
      if (
           ( $#{ $c1->{res} } >= 0 )
        && ( $#{ $c2->{res} } >= 0 )
        && ( ( !defined $cidA )
          || ( $cidA eq $c1->{id} )
          || ( $cidA eq $c2->{id} ) )
        && ( ( !defined $cidB )
          || ( $cidB eq $c1->{id} )
          || ( $cidB eq $c2->{id} ) )
        )
      {

        my $AAcontacts = 0;
        my $ca1 = $c1->{atom};
        my $ca2 = $c2->{atom};
        foreach my $r1 ( @{ $c1->{res} } ) {
          if ( $r1->{des} ne "NATRO" ) {
            for ( my $ai1 = $r1->{start} ; $ai1 <= $r1->{end} ; $ai1++ ) {
              my $a1 = $ca1->[$ai1];
              if (
                   ( $a1->{aux1} >= 0.0 )
                && ( !$a1->{hyd} )
                && ( ( !$SConly ) || ( !$a1->{bb} ) )
                && ( ( !$CACBonly )
                  || ( $a1->{atomname} eq "CA" )
                  || ( $a1->{atomname} eq "CB" ) )
                )
              {
                my $x1 = $a1->{xcoor};
                my $y1 = $a1->{ycoor};
                my $z1 = $a1->{zcoor};
                foreach my $r2 ( @{ $c2->{res} } ) {
                  if ( $r2->{des} ne "NATRO" ) {
                    for ( my $ai2 = $r2->{start} ; $ai2 <= $r2->{end} ; $ai2++ )
                    {
                      my $a2 = $ca2->[$ai2];
                      if (
                           ( $a2->{aux1} >= 0.0 )
                        && ( !$a2->{hyd} )
                        && ( ( !$SConly ) || ( !$a2->{bb} ) )
                        && ( ( !$CACBonly )
                          || ( $a2->{atomname} eq "CA" )
                          || ( $a2->{atomname} eq "CB" ) )
                        )
                      {
                        my $x2 = $a2->{xcoor};
                        my $y2 = $a2->{ycoor};
                        my $z2 = $a2->{zcoor};
                        my $SqDist = ( $x1 - $x2 ) * ( $x1 - $x2 );
                        $SqDist += ( $y1 - $y2 ) * ( $y1 - $y2 );
                        $SqDist += ( $z1 - $z2 ) * ( $z1 - $z2 );
                        $AAcontacts++ if ( $SqDist < $SqDistCut );
                      }
                    }
                  }
                }
              }
            }
          }
        }

        my $h = {};
        $h->{AAcontacts} = $AAcontacts;
        $h->{cidA} = $c1->{id};
        $h->{cidB} = $c2->{id};
        push( @{$AAlist}, $h );

      }
    }
  }

  return $AAlist;
}

## method: setDesignLoops()
## make residues defined as loops designable (ie. ALLAA)

sub setDesignLoops {
  my $self = shift;

  die "No chains defined for setDesignLoops"
    if ( !defined $self->activeChains() );

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
        my $ri = $c->{res}->[$ir];
        $ri->{des} = "ALLAA" if ( $ri->{loop} );
      }
    }
  }

  return;
}

## method: setAllDesType(cid,type)
## Set all residues to a particular type of design

sub setAllDesType {
  my $self = shift;
  my $cid = shift;
  my $type = shift;

  if ( defined $self->activeChains() ) {
    foreach my $c ( @{ $self->activeChains() } ) {
      if ( ( !defined $cid ) || ( $c->{id} eq $cid ) ) {
        if ( defined $c->{res} ) {
          for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
            $c->{res}->[$ir]->{des} = $type;
          }
        }
      }
    }
  }

  return;
}

## method: changeDesType(oldType,newType,[cid])
## Change all residues of a particular design type to a new type
## Operates on all chains if no chain id is given

sub changeDesType {
  my $self = shift;
  my $oldT = shift;
  my $newT = shift;
  my $cid = shift;

  die "No chains defined for changeDesType"
    if ( !defined $self->activeChains() );

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( ( ( !defined $cid ) || ( $cid eq $c->{id} ) )
      && ( $#{ $c->{res} } >= 0 ) )
    {
      for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
        $c->{res}->[$ir]->{des} = $newT
          if ( $c->{res}->[$ir]->{des} =~ $oldT );
      }
    }
  }

  return;
}

## method: setResDesType(cid,resnum,type)
## Define type of design to do with a particular residue

sub setResDesType {
  my $self = shift;
  my $cid = shift;
  my $ires = shift;
  my $type = shift;

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      if ( ( !defined $cid ) || ( $cid eq $c->{id} ) ) {
        my $r = $self->getResidueInChain( $ires, $c );
        if ( defined $r ) {
          $r->{des} = $type;
        }
      }
    }
  }

  return;
}

## method: type = reportResDesType(cid,resnum)
## Define type of design to do with a particular residue

sub reportResDesType {
  my $self = shift;
  my $cid = shift;
  my $ires = shift;

  my $type;

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      if ( ( !defined $cid ) || ( $cid eq $c->{id} ) ) {
        my $r = $self->getResidueInChain( $ires, $c );
        if ( defined $r ) {
          $type = $r->{des};
        }
      }
    }
  }

  return $type;
}

## method: setResnameDesType(cid,resname,type)
## Define type of design to do for all positions which are currently resname

sub setResnameDesType {
  my $self = shift;
  my $cid = shift;
  my $resname = shift;
  my $type = shift;

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      if ( ( !defined $cid ) || ( $cid eq $c->{id} ) ) {
        my $ra = $c->{res};
        for ( my $ir = 0 ; $ir <= $#{$ra} ; $ir++ ) {
          my $r = $c->{res}->[$ir];
          $r->{des} = $type
            if ( lc( $r->{name} ) eq lc($resname) );
        }
      }
    }
  }

  return;
}

## method: clearDesignRes()
## set all residues to NATRO

sub clearDesignRes {
  my $self = shift;

  die "No chains defined for clearDesignRes"
    if ( !defined $self->activeChains() );

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
        $c->{res}->[$ir]->{des} = "NATRO";
      }
    }
  }

  return;
}

## method: removeChain()
## polymorphism

sub removeChain {
  my $self = shift;
  my $cid = shift;

  $self->SUPER::removeChain($cid);
  $self->renumberAcrossChains() if ($self->{par}->{doRenumber});
  undef $self->{par}->{currFile};

  return;
}

## method: addInfoLine()
## polymorphism

sub addInfoLine {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};
  return $self->SUPER::addInfoLine(@args);
}

## method: enantiomer()
## polymorphism

sub enantiomer {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};
  return $self->SUPER::enantiomer(@args);
}

## method: fillCoorFromPDB()
## polymorphism

sub fillCoorFromPDB {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};
  return $self->SUPER::fillCoorFromPDB(@args);
}

## method: merge()
## polymorphism

sub merge {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};
  return $self->SUPER::merge(@args);
}

## method: move(dx,dy,dz,cid)
## polymorphism

sub move {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};
  return $self->SUPER::move(@args);
}

## method: matrixOperation(dx,dy,dz,cid)
## polymorphism

sub matrixOperation {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};
  return $self->SUPER::matrixOperation(@args);
}

## method: center()
## polymorphism

sub center {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};
  return $self->SUPER::center(@args);
}

## method: setChain()
## polymorphism

sub setChain {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};
  return $self->SUPER::setChain(@args);
}

## method: substituteAla()
## polymorphism

sub substituteAla {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};
  return $self->SUPER::substituteAla(@args);
}

## method: setSSbonds()
## polymorphism

sub setSSbonds {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};
  return $self->SUPER::setSSbonds(@args);
}

## method: findSSbonds()
## polymorphism

sub findSSbonds {
  my $self = shift;
  my @args = @_;

  undef $self->{par}->{currFile};
  return $self->SUPER::findSSbonds(@args);
}

## method: getLoopScore()
## return the total score of loop residues

sub getLoopScore {
  my $self = shift;

  return undef
    if ( !defined $self->activeChains() );

  my $E = 0.0;
  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
        $E += $c->{res}->[$ir]->{score} if ( $c->{res}->[$ir]->{loop} );
      }
    }
  }

  return $E;
}

## method: getDesScore()
## return a score hash for designable residues

sub getDesScore {
  my $self = shift;

  return undef
    if ( !defined $self->activeChains() );

  my $s = {};
  $s->{reslist} = ();
  my $numDesRes = 0;
  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      foreach my $r ( @{ $c->{res} } ) {
        if ( $r->{des} eq "ALLAA" ) {
          $numDesRes++;
          push( @{ $s->{reslist} }, ( $r->{num} . $c->{id} ) );
          if ( defined $r->{score} ) {
            foreach my $k ( keys %{ $r->{score} } ) {
              $s->{score}->{$k} = 0.0 if ( !defined $s->{score}->{$k} );
              $s->{score}->{$k} += $r->{score}->{$k};
              $s->{max}->{$k} = $r->{score}->{$k}
                if ( ( !defined $s->{max}->{$k} )
                || ( $s->{max}->{$k} < $r->{score}->{$k} ) );
              $s->{min}->{$k} = $r->{score}->{$k}
                if ( ( !defined $s->{min}->{$k} )
                || ( $s->{min}->{$k} > $r->{score}->{$k} ) );
            }
          }
          if ( defined $r->{VsPDB} ) {
            $s->{maxVsPDB} = {} if ( !defined $s->{maxVsPDB} );
            foreach my $k ( keys %{ $r->{VsPDB} } ) {
              $s->{VsPDB}->{$k} = 0.0 if ( !defined $s->{VsPDB}->{$k} );
              $s->{VsPDB}->{$k} += $r->{VsPDB}->{$k};
              $s->{maxVsPDB}->{$k} = $r->{VsPDB}->{$k}
                if ( ( !defined $s->{maxVsPDB}->{$k} )
                || ( $s->{maxVsPDB}->{$k} < $r->{VsPDB}->{$k} ) );
              $s->{minVsPDB}->{$k} = $r->{VsPDB}->{$k}
                if ( ( !defined $s->{minVsPDB}->{$k} )
                || ( $s->{minVsPDB}->{$k} > $r->{VsPDB}->{$k} ) );
            }
          }
        }
      }
    }
  }

  $s->{numDesRes} = $numDesRes;
  return $s;
}

## method: donateLoopResidues()
## Put the current structure's residue coordinates onto the
## coordinates of the RosettaMolecule passed in, for residues where
## self's "loop" is set. Note that the number of atoms
## in corresponding residues need not match

sub donateLoopResidues {
  my $self = shift;
  my $inMol = shift;

  die "No chains defined for donateLoopResidues"
    if ( !defined $self->activeChains() );

  foreach my $sc ( @{ $self->activeChains() } ) {
    my $dc = $inMol->getChain( $sc->{id} );
    die "Cannot lookup chain %s", $sc->{id}
      if ( !defined $dc->{res} );

    # Save the destination residues first, then
    # write over some of them with the source
    my %newcoors;
    foreach my $r ( @{ $dc->{res} } ) {
      my $trec = { res => $r, atom => $dc->{atom} };
      $newcoors{ $r->{num} } = $trec;
    }
    foreach my $r ( @{ $sc->{res} } ) {
      if ( $r->{loop} ) {
        my $trec = { res => $r, atom => $sc->{atom} };
        $newcoors{ $r->{num} } = $trec;
      }
    }

    my $newatom = ();
    my $newres = ();
    my $ainx = 1;
    foreach my $k ( sort { $a <=> $b } keys %newcoors ) {
      my $r = $newcoors{$k}->{res};
      my $a = $newcoors{$k}->{atom};
      my $rrec = {};
      %{$rrec} = %{$r};
      push( @{$newres}, $rrec );
      $rrec->{start} = $#{$newatom} + 1;
      for ( my $ia = $r->{start} ; $ia <= $r->{end} ; $ia++ ) {
        my $arec = {};
        %{$arec} = %{ $a->[$ia] };
        $arec->{atominx} = $ainx++;
        push( @{$newatom}, $arec );
      }
      $rrec->{end} = $#{$newatom};
    }
    $dc->{atom} = $newatom;
    $dc->{res} = $newres;
  }

  $inMol->{segmentlist} = undef;
  $inMol->_coorCache();
  undef $inMol->{par}->{currFile};

  return;
}

## method: countLoopScontacts([startRes] [,stopRes])
## return the number of sidechain-sidechain contacts
## between loops and a given region (set to all
## non-loop regions, if undefined). Distance cutoff
## is "interfaceDist".

sub countLoopScontacts {
  my $self = shift;
  my $startRes = shift;
  my $stopRes = shift;

  return undef
    if ( !defined $self->activeChains() );

  my $distCut = $self->{par}->{interfaceDist};
  my $distCutSq = $distCut * $distCut;

  my $nc = 0;
  foreach my $Lc ( @{ $self->activeChains() } ) {
    my $Lra = $Lc->{res};
    my $Laa = $Lc->{atom};
    if ( $#{$Lra} >= 0 ) {
      foreach my $Lr ( @{$Lra} ) {
        if ( ( defined $Lr->{loop} ) && ( $Lr->{loop} ) ) {
          for ( my $Lai = $Lr->{start} ; $Lai <= $Lr->{end} ; $Lai++ ) {
            my $La = $Laa->[$Lai];
            if ( ( !$La->{hyd} ) && ( !$La->{bb} ) ) {

              foreach my $Rc ( @{ $self->activeChains() } ) {
                my $Rra = $Rc->{res};
                my $Raa = $Rc->{atom};
                if ( $#{$Rra} >= 0 ) {
                  foreach my $Rr ( @{$Rra} ) {
                    if ( ( defined $Rr->{loop} ) && ( !$Rr->{loop} ) ) {
                      for (
                        my $Rai = $Rr->{start} ;
                        $Rai <= $Rr->{end} ;
                        $Rai++
                        )
                      {
                        my $Ra = $Raa->[$Rai];
                        if ( ( !$Ra->{hyd} ) && ( !$Ra->{bb} ) ) {

                          my $distSq =
                            ( $La->{xcoor} - $Ra->{xcoor} ) *
                            ( $La->{xcoor} - $Ra->{xcoor} );
                          $distSq +=
                            ( $La->{ycoor} - $Ra->{ycoor} ) *
                            ( $La->{ycoor} - $Ra->{ycoor} );
                          $distSq +=
                            ( $La->{zcoor} - $Ra->{zcoor} ) *
                            ( $La->{zcoor} - $Ra->{zcoor} );
                          $nc++ if ( $distSq < $distCutSq );

                        }
                      }
                    }
                  }
                }
              }
            }    # matches if not hyd or backbone....
          }
        }
      }
    }
  }

  return $nc;
}

## method: minLoopDistToRegion(startRes,stopRes)
## return the minimum distance from a certain region
## to the closest loop region

sub minLoopDistToRegion {
  my $self = shift;
  my $startRes = shift;
  my $stopRes = shift;

  return undef
    if ( !defined $self->activeChains() );

  $startRes = 1 if ( !defined $startRes );
  $stopRes = $self->numres() if ( !defined $stopRes );

  my $closest = 9999;
  my $loopList = ();
  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      foreach my $rL ( @{ $c->{res} } ) {
        push( @{$loopList}, $rL )
          if ( ( defined $rL->{loop} ) && ( $rL->{loop} ) );
      }
    }
  }

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      foreach my $rR ( @{ $c->{res} } ) {
        if ( ( $rR->{num} >= $startRes ) && ( $rR->{num} <= $stopRes ) ) {
          foreach my $rL ( @{$loopList} ) {
            my $dist = $self->minDistance( $rL, $rR );
            $closest = $dist if ( $dist < $closest );
          }
        }
      }
    }
  }

  return $closest;
}

## method: setLoop(cid,resnum,[type])
## Set (or clear) the loop field for a particular residue

sub setLoop {
  my $self = shift;
  my $cid = shift;
  my $ires = shift;
  my $type = shift;

  $type = 1 if ( !defined $type );

  my $c = $self->getChain($cid);
  if ( defined $c->{res} ) {
    my $r = $self->getResidueInChain( $ires, $c );
    if ( defined $r ) {
      $r->{loop} = $type;
    }
  }

  return;
}

## method: zapCoordinates([chain])
## sets all coordinates in the current structure
## to zero and occupancy to -1

sub zapCoordinates {
  my $self = shift;
  my $chainid = shift;
  my $start = shift;
  my $stop = shift;

  die "No chains defined for zapCoordinates"
    if ( !defined $self->activeChains($chainid) );

  undef $self->{par}->{currFile};

  $self->SUPER::zapCoordinates($chainid,$start,$stop);

  foreach my $c ( @{ $self->activeChains($chainid) } ) {
    foreach my $a ( @{ $c->{atom} } ) {
      if ( ( ( !defined $start ) || ( $a->{resnum} >= $start ) )
        && ( ( !defined $stop ) || ( $a->{resnum} <= $stop ) ) )
      {
        $a->{aux1} = -1.0;
      }
    }
    $self->_coorCache( $c->{id} );
  }

  return;
}

## method: _rosettaSetup([parameters])
## General Rosetta setup using the current parameters
## Note: should be called AFTER _setupOutfiles, since _setupOutfiles can change nstruct

sub _rosettaSetup {
  my $self = shift;
  my $tempInfiles = shift;

  my $setupCmd = "";

#  $setupCmd .= " -profile";

  $setupCmd .= " -nstruct " . $self->{par}->{nstruct};

  $setupCmd .= " -use_pdbseq";
  $setupCmd .= " -read_all_chains";
  $setupCmd .= " -try_both_his_tautomers";
  $setupCmd .= " -no_his_his_pairE";
  $setupCmd .= " -use_electrostatic_repulsion";
#  $setupCmd .= " -flip_symmetric_sidechains";
#  $setupCmd .= " -vary_bond_geometry_no_hydrogens";

  # jk doesn't work for idealize due to unusual energy landscape
  $setupCmd .= " -use_non_monotone_line_search" if ( $self->{par}->{monotone_line_search} );
#  $setupCmd .= " -use_inexact_line_search";

#  $setupCmd .= " -termini";

  #  $setupCmd.=" -find_disulf";
  #  $setupCmd.=" -norepack_disulf";
  #  $setupCmd.=" -no_filters";
  $setupCmd .= " -scorefile scorefile";
  $setupCmd .= " -ds_outpdb_only";
  $setupCmd .= " -preserve_header";
  $setupCmd .= " -no_option_display";

  # gzipped output appears to work, at least for dockpert on syd
  $setupCmd .= " -output_pdb_gz";

  # write a resfile if needed
  my $resfile = $self->{par}->{resfile};
  if ( ! &GenUtil::zexists($resfile) ) {
    $self->writeResfile($resfile);
    printf $tempInfiles "%s\n", $resfile;
  }
  $setupCmd .= " -resfile " . $resfile if ( &GenUtil::zexists($resfile) );

# note: output arrays will only be filled with validGzip files
#       if this command-line option is removed, remove that requirement as well!

  if ( $self->{par}->{multi_chain} ) {
      $setupCmd .= " -multi_chain";
  }

  if ( $self->{par}->{decoystats} ) {
      $setupCmd .= " -decoystats";
      if (defined $self->{par}->{gu_exempt}) {
	  $setupCmd .= " -gu_exempt ".$self->{par}->{gu_exempt};
      }
  }

  $setupCmd .= " -interface_ds"
    if ( $self->{par}->{interface_ds} );

  $setupCmd .= " -tight_core_analysis -tight_core_thres 5.0"
    if ( $self->{par}->{tight_core} );

  $setupCmd .= " -soft_rep" if ( $self->{par}->{soft_rep} );

  if ( $self->{par}->{exp_water} ) {
    $setupCmd .= " -explicit_h2o -solvate";
  } else {
    $setupCmd .= " -gen_born"
      if ( $self->{par}->{gb} && ! $self->{par}->{incorporate_hotspot} );
  }

  my $fastaFile = $self->{par}->{protein} . "A.fasta";
  $self->writeFasta($fastaFile);
  printf $tempInfiles "%s\n", $fastaFile;

  $setupCmd .= " -use_input_sc"    if ( $self->{par}->{use_input_sc} );
  $setupCmd .= " -read_hetero_h2o" if ( $self->{par}->{use_het_water} );
  $setupCmd .= " -jk_interface"    if ( $self->{par}->{jk_interface} );
  $setupCmd .= " -fake_native"    if ( $self->{par}->{fake_native} );

  if ( ! $self->{par}->{identifyInterfaceUNS} ) {

    $setupCmd .= " -rot_pert_input";
    if ( $self->{par}->{ex_rot} >= 1 ) {
      $setupCmd .= " -ex1";
      $setupCmd .= " -ex1aro";
      if ( $self->{par}->{ex_rot} < 2 ) {
	$setupCmd .= " -ex2aro_only";
      } else {
	$setupCmd .= " -ex2";
#	$setupCmd .= " -minimize_rot";
#	$setupCmd .= " -extrachi_cutoff -1";
	if ( $self->{par}->{ex_rot} > 2 ) {
	  $setupCmd .= " -ex3";
	  if ( $self->{par}->{ex_rot} == 4 ) {
	    $setupCmd .= " -ex4";
	  }
	}
      }
    } elsif ( $self->{par}->{ex_rot} > 0 ) {
      $setupCmd .= " -rot_pert -pert_size 5";
      $setupCmd .= " -pert_acc_prob " . $self->{par}->{ex_rot};
    }

  } else {
    $setupCmd .= " -rot_pert_input";
    $setupCmd .= " -ex1";
    $setupCmd .= " -ex1aro";
    $setupCmd .= " -ex2";
    $setupCmd .= " -ex3";
    $setupCmd .= " -rot_pert -pert_size 5";
    $setupCmd .= " -pert_acc_prob 0.9";
    $setupCmd .= " -extrachi_cutoff -1";
  }

  my $pathfile = $self->{par}->{pathfile};
  $self->writePathFile();
  printf $tempInfiles "%s\n", $pathfile;
  $setupCmd .= " -paths " . $pathfile;

  printf $tempInfiles "%s\n", &UWutil::reportJobFilename()
    if ( &UWutil::reportClusterCpus() > 0 );
  printf $tempInfiles "%s\n", $self->{par}->{stdout};
  printf $tempInfiles "%s\n", $self->{par}->{stderr};

  return $setupCmd;
}

## method: _setupOutfiles([$startArr],[noSeriesCode])

sub _setupOutfiles {
  my $self = shift;
  my $startArr = shift;
  my $noSeriesCode = shift;
  my $noSubdirs = shift;

  my $seriesCode=$self->{par}->{series};
  $seriesCode="" if ((defined $noSeriesCode) && ($noSeriesCode == 1));

  # Note: condor output file location will be confused if we don't have subdirs.
  # Always use subdirs unless parameter says otherwise (ie. a "curr" step, which doesn't use condor)
  $noSubdirs=0 if (! defined $noSubdirs);
#      my $nStart=$#{$startArr}+1;
#      $nStart=1 if ( $nStart < 1 );
#      if (($nStart * $self->{par}->{nstruct}) <= 200) {
#	  $noSubdirs=1;
#      } else {
#	  $noSubdirs=0;
#      }
#  }

  my $num_processes = &UWutil::reportClusterCpus();
  $num_processes = 100 if ( $num_processes == 1 );

  # If the number of input files is small and nstruct is large,
  # make multiple copies of the (few) input files so that we can
  # distribute input files across cpus (and directories)
  my @effective_startArr = @{$startArr};

#  if ( ! $noSubdirs ) {
#    while ( ($#effective_startArr + 1) < $num_processes ) {
#      @effective_startArr = (@effective_startArr, @effective_startArr);
#      $self->{par}->{nstruct} = int(($self->{par}->{nstruct} / 2.0) + 0.9);
#    }
#  }

  $num_processes = $#effective_startArr + 1 if ( $num_processes > $#effective_startArr );
  $num_processes = 1 if ( $noSubdirs );

  my $filelist = "infiles.list";
  my $flist = &GenUtil::getOutputFile($filelist);
  my $ename = $self->{par}->{outdir} . $self->{par}->{outputPDBs};
  &GenUtil::makeDir($self->{par}->{outdir}) if ( !-d $self->{par}->{outdir} );
  my $elist = &GenUtil::getOutputFile($ename);
  my $earr = ();
  my $fnum = 0;
  my $numOutfiles = 0;

  # If ename exists, fill earr by reading it, rather than from scratch
  if ( &GenUtil::zexists($ename ) ) {
    my $ein = &GenUtil::getInputFile($ename);
    while (<$ein>) {
      my $inline=$_;
      chomp($inline);
      my @p = split( ' ', $inline );
      push( @{$earr}, $p[1] );
    }
    undef $ein;
    return ($filelist,$earr,$num_processes);
  }

  if (! $noSubdirs) {
    for ( my $i = 0 ; $i < $num_processes ; $i++ ) {
      my $fulldir = $self->{par}->{outdir} . $i . "/";
      &GenUtil::makeDir($fulldir) if ( !-d $fulldir );
      &GenUtil::remove("infiles.list".$i) if ( &GenUtil::exists("infiles.list".$i) );
    }
  }

  foreach my $startfile ( @effective_startArr ) {

    my $process_num = $fnum % $num_processes;
    my $process_filelist = "infiles.list".$process_num;
    my $process_flist = &GenUtil::getAppendFile($process_filelist);

    my $outdir="";
    $outdir = $process_num . "/" if (! $noSubdirs);
    my $fulldir = $self->{par}->{outdir}.$outdir;
#    &GenUtil::makeDir($fulldir) if ( !-d $fulldir );
    my $fname = $outdir. "temp" . $fnum . ".pdb";

    my $startfile_name = $startfile->{filename};
    ( $startfile_name = $startfile_name ) =~ s/.gz$//gi;
    printf $flist "%s %s\n", $startfile_name, $fname;
    printf $process_flist "%s %s\n", $startfile_name, $fname;

    my $efile;

    for ( my $i = 1 ; $i <= $self->{par}->{nstruct} ; $i++ ) {

      $efile = $fulldir . $seriesCode . "temp" . $fnum . "_" . &GenUtil::zPad( $i, 4 ) . ".pdb";

      # Support for checkpointing via "in_progress" files
      &GenUtil::remove($efile.".in_progress") if ( &GenUtil::exists($efile.".in_progress") );

      # Support for checkpointing system (used in hotspot incorporation, maybe elsewhere later)
      my $checkpoint_file = substr( $efile, 0, -4) . ".checkpoint";
      if ( -e $checkpoint_file ) {
	&GenUtil::remove($efile) if ( &GenUtil::exists($efile) );
      }

      my $cfile = $fulldir . $seriesCode;
      $cfile .= "temp" . $fnum . "_" . &GenUtil::zPad( $i, 4 ) . ".clean.pdb";
      printf $elist "%s %s %s\n", $cfile, $efile, $startfile->{filename};
      $numOutfiles++;

      # print each name only once, regardless of nstruct
      last if ( $self->{par}->{incorporate_hotspot} );

    }
    push( @{$earr}, $efile );
    $fnum++;
  }
  undef $flist;
  undef $elist;

  # Note: handled above on the fly, instead of retroactively
  #  &clean_infile_list($filelist);

  return ($filelist,$earr,$num_processes);
}

## method: _setupFromStartArr([$startArr])
## Setup

sub _setupFromStartArr {
  my $self = shift;
  my $startArr = shift;

  $self->writePDB( $self->{par}->{protein} . ".pdb", "coor" );
  if ( ( !defined $startArr ) || ( $#{$startArr} == -1 ) ) {
    my $f = {};
    $f->{filename} = $self->{par}->{protein} . ".pdb";
    push( @{$startArr}, $f );
  }

  return $startArr;
}

## method: _makeFragLib([outStdErr])
## uses make_fragments.pl to build a fragment library

sub _makeFragLib {
  my $self = shift;
  my $outStdErr = shift;

  # Rename the PDB file so make_fragments won't see it
  my $currPDB = $self->{par}->{protein} . ".pdb";
  my $tempPDB = "TEMPFILE.pdb";
  &GenUtil::rename( $currPDB, $tempPDB );

  my $fastaFile = $self->{par}->{protein} . "A.fasta";
  $self->writeFasta($fastaFile);

  &GenUtil::log("calling makeFragments.pl");

  my $cmd = $nnmakeExec;
  $cmd .= " -nohoms " . $fastaFile;

  &UWutil::runCmd($cmd);

  &GenUtil::rename( $tempPDB, $currPDB );

  return;
}

## method: buildRandFragLib()
## build a random fragment library for loop modeling

sub buildRandFragLib {
  my $self = shift;

  &GenUtil::log("reading Vall");

  my $vall_in = &GenUtil::getInputFile($Vall_name);
  my @vall_lines = ();
  while (<$vall_in>) {
    chomp;
    push( @vall_lines, $_ );
  }
  undef $vall_in;

  &GenUtil::log("selecting random 3mer fragments for library");

  my $des = {};
  $des->{"H"} = 60;
  $des->{"E"} = 60;
  $des->{"L"} = 80;
  my $frags_needed = $des->{"H"} + $des->{"E"} + $des->{"L"};
  my $fname3 = "aa" . $self->{par}->{protein} . "A03_05.200_v1_3";
  my $outlib3 = &GenUtil::getOutputFile($fname3);
  my $c = $self->getChain();
  die "No chains defined for buildRandFragLib" if ( !defined $c );

  for ( my $ir = 2 ; $ir <= $#{ $c->{res} } ; $ir++ ) {

    printf $outlib3 " position:          %3d neighbors:          %3d\n\n",
      $ir - 1, $frags_needed;
    my $found = {};
    $found->{"H"} = 0;
    $found->{"E"} = 0;
    $found->{"L"} = 0;
    my $fragnum = 0;
    while ( $fragnum < $frags_needed ) {
      my $r = int( rand($#vall_lines) ) + 1;
      my $inline = $vall_lines[$r];
      my @p = split( ' ', $inline );
      if ( $p[4] >= 2 ) {
        my $SS = $p[2];
        if ( $found->{$SS} < $des->{$SS} ) {
          $found->{$SS}++;
          $fragnum++;
          for ( my $i = 0 ; $i < 3 ; $i++ ) {
            my @p = split( ' ', $vall_lines[$r] );
            printf $outlib3 " %4s %1s  %4d %s %s", substr( $p[0], 0, -1 ),
              substr( $p[0], -1, 1 ), $p[3], $p[1], $p[2];
            printf $outlib3 " %8.3f %8.3f %8.3f", $p[9], $p[10], $p[11];
            printf $outlib3 "    0.000    0.000     0.000 0     0.000";
            printf $outlib3 " P%3d F%3d\n", $ir - 1, $fragnum;
            $r++;
          }
          printf $outlib3 "\n";
        }
      }
    }
  }
  undef $outlib3;

  $des->{"H"} = 7;
  $des->{"E"} = 7;
  $des->{"L"} = 11;
  $frags_needed = $des->{"H"} + $des->{"E"} + $des->{"L"};
  my $fname9 = "aa" . $self->{par}->{protein} . "A09_05.200_v1_3";
  my $outlib9 = &GenUtil::getOutputFile($fname9);
  for ( my $ir = 8 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
    printf $outlib9 " position:          %3d neighbors:          200\n\n",
      $ir - 7, $frags_needed;
    my $found = {};
    $found->{"H"} = 0;
    $found->{"E"} = 0;
    $found->{"L"} = 0;
    my $fragnum = 0;
    while ( $fragnum < $frags_needed ) {
      my $r = int( rand($#vall_lines) ) + 1;
      my $inline = $vall_lines[$r];
      my @p = split( ' ', $inline );
      if ( $p[4] >= 8 ) {
        my $SS = $p[2];
        if ( $found->{$SS} < $des->{$SS} ) {
          $found->{$SS}++;
          $fragnum++;
          for ( my $i = 0 ; $i < 9 ; $i++ ) {
            my @p = split( ' ', $vall_lines[$r] );
            printf $outlib9 " %4s %1s  %4d %s %s", substr( $p[0], 0, -1 ),
              substr( $p[0], -1, 1 ), $p[3], $p[1], $p[2];
            printf $outlib9 " %8.3f %8.3f %8.3f", $p[9], $p[10], $p[11];
            printf $outlib9 "    0.000    0.000     0.000 0     0.000";
            printf $outlib9 " P%3d F%3d\n", $ir - 7, $fragnum;
            $r++;
          }
          printf $outlib9 "\n";
        }
      }
    }
    for ( my $i = $frags_needed + 1 ; $i <= 200 ; $i++ ) {

      # Put in 175 "filler" 9mer fragments
      for ( my $j = 1 ; $j <= 9 ; $j++ ) {
        printf $outlib9 " 1jnk _  %4d A L", $j;
        printf $outlib9 " %8.3f %8.3f %8.3f", 0.0, 0.0, 0.0;
        printf $outlib9 "    0.000    0.000     0.000 0     0.000";
        printf $outlib9 " P%3d F%3d\n", $ir - 7, $i;
      }
      printf $outlib9 "\n";
    }
  }
  undef $outlib9;

  return;
}

## method: buildDirectedRangFragLib()
## build a random fragment library with only one type of SSdes
## Note: SSdes should be "H", "E", or "L"

sub buildDirectedRandFragLib {
  my $self = shift;
  my $SSdes = shift;

  &GenUtil::log("reading Vall");

  my $vall_in = &GenUtil::getInputFile($Vall_name);
  my @vall_lines = ();
  while (<$vall_in>) {
    chomp;
    push( @vall_lines, $_ );
  }
  undef $vall_in;

  &GenUtil::log( "selecting random 3mer fragments of type " . $SSdes );

  my $frags_needed = 200;
  my $fname3 = "aa" . $self->{par}->{protein} . "A03_05.200_v1_3";
  my $outlib3 = &GenUtil::getOutputFile($fname3);
  my $c = $self->getChain();
  die "No chains defined for buildRandFragLib" if ( !defined $c );

  for ( my $ir = 2 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
    printf $outlib3 " position:          %3d neighbors:          %3d\n\n",
      $ir - 1, $frags_needed;
    my $found = 0;
    my $fragnum = 0;
    while ( $fragnum < $frags_needed ) {
      my $r = int( rand($#vall_lines) ) + 1;
      my $SSmatch = 1;
      my $inline = $vall_lines[$r];
      my @p = split( ' ', $inline );
      $SSmatch = 0 if ( $p[4] < 2 );
      for ( my $i = 0 ; ( ( $i < 3 ) && ($SSmatch) ) ; $i++ ) {
        my $inline = $vall_lines[ $r + $i ];
        my @p = split( ' ', $inline );
        $SSmatch = 0 if ( $p[2] ne $SSdes );
      }
      if ($SSmatch) {
        $found++;
        $fragnum++;
        for ( my $i = 0 ; $i < 3 ; $i++ ) {
          my @p = split( ' ', $vall_lines[$r] );
          printf $outlib3 " %4s %1s  %4d %s %s", substr( $p[0], 0, -1 ),
            substr( $p[0], -1, 1 ), $p[3], $p[1], $p[2];
          printf $outlib3 " %8.3f %8.3f %8.3f", $p[9], $p[10], $p[11];
          printf $outlib3 "    0.000    0.000     0.000 0     0.000";
          printf $outlib3 " P%3d F%3d\n", $ir - 1, $fragnum;
          $r++;
        }
        printf $outlib3 "\n";
      }
    }
  }
  undef $outlib3;

  $frags_needed = 25;
  my $fname9 = "aa" . $self->{par}->{protein} . "A09_05.200_v1_3";
  my $outlib9 = &GenUtil::getOutputFile($fname9);
  for ( my $ir = 8 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
    printf $outlib9 " position:          %3d neighbors:          200\n\n",
      $ir - 7, $frags_needed;
    my $found = 0;
    my $fragnum = 0;
    while ( $fragnum < $frags_needed ) {
      my $r = int( rand($#vall_lines) ) + 1;
      my $SSmatch = 1;
      my $inline = $vall_lines[$r];
      my @p = split( ' ', $inline );
      $SSmatch = 0 if ( $p[4] < 8 );
      for ( my $i = 0 ; ( ( $i < 9 ) && ($SSmatch) ) ; $i++ ) {
        my $inline = $vall_lines[ $r + $i ];
        my @p = split( ' ', $inline );
        $SSmatch = 0 if ( $p[2] ne $SSdes );
      }
      if ($SSmatch) {
        $found++;
        $fragnum++;
        for ( my $i = 0 ; $i < 9 ; $i++ ) {
          my @p = split( ' ', $vall_lines[$r] );
          printf $outlib9 " %4s %1s  %4d %s %s", substr( $p[0], 0, -1 ),
            substr( $p[0], -1, 1 ), $p[3], $p[1], $p[2];
          printf $outlib9 " %8.3f %8.3f %8.3f", $p[9], $p[10], $p[11];
          printf $outlib9 "    0.000    0.000     0.000 0     0.000";
          printf $outlib9 " P%3d F%3d\n", $ir - 7, $fragnum;
          $r++;
        }
        printf $outlib9 "\n";
      }
    }
    for ( my $i = $frags_needed + 1 ; $i <= 200 ; $i++ ) {

      # Put in 175 "filler" 9mer fragments
      for ( my $j = 1 ; $j <= 9 ; $j++ ) {
        printf $outlib9 " 1jnk _  %4d A L", $j;
        printf $outlib9 " %8.3f %8.3f %8.3f", 0.0, 0.0, 0.0;
        printf $outlib9 "    0.000    0.000     0.000 0     0.000";
        printf $outlib9 " P%3d F%3d\n", $ir - 7, $i;
      }
      printf $outlib9 "\n";
    }
  }
  undef $outlib9;

  return;
}

## method: writeFasta(outfile,[comment])
## write a Fasta file

sub writeFasta {
  my $self = shift;
  my $fastaFname = shift;
  my $comment = shift;

  my $fasta = &GenUtil::getOutputFile($fastaFname);
  if ( defined $comment ) {
    printf $fasta "> %s\n", $comment;
  }
  else {
    printf $fasta "> HEADER\n";
  }

  my $seq = Sequence::new($self);
  printf $fasta "%s\n", $seq->abbrevSeq();
  undef $fasta;

  return;
}

## method: writeLoopFile()
## write a loop file

sub writeLoopFile {
  my $self = shift;
  my $loopFname = shift;

  die "No chains defined for writeLoopFile"
    if ( !defined $self->activeChains() );

  $loopFname = $self->{par}->{protein} . ".loops" if ( !defined $loopFname );
  my $lfile = &GenUtil::getOutputFile($loopFname);

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
        if ( $c->{res}->[$ir]->{loop} == 1 ) {
          my $ri = $c->{res}->[$ir];
          my $Lstart = $ri->{num};
          my $Lend = $ri->{num};
          my $Llen = 1;
          $ir++;
          while ( ( defined $c->{res}->[$ir] )
            && ( $c->{res}->[$ir]->{loop} == 1 ) )
          {
            $Lend = $c->{res}->[$ir]->{num};
            $Llen++;
            $ir++;
          }
          printf $lfile "%3d  %3d  %3d\n", $Llen, $Lstart, $Lend;
        }
      }
    }
  }
  undef $lfile;

  return;
}

## method: readLoopFile()
## read a loop file

sub readLoopFile {
  my $self = shift;
  my $loopFname = shift;

  die "No chains defined for readLoopFile"
    if ( !defined $self->activeChains() );

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
        $c->{res}->[$ir]->{loop} = 0;
      }
    }
  }

  $loopFname = $self->{par}->{protein} . ".loops" if ( !defined $loopFname );
  if ( &GenUtil::exists($loopFname) ) {
    my $lfile = &GenUtil::getInputFile($loopFname);
    while (<$lfile>) {
      chomp;
      my ( $Lnumres, $LStart, $LStop ) = split( ' ', $_ );
      foreach my $c ( @{ $self->activeChains() } ) {
        if ( $#{ $c->{res} } >= 0 ) {
          for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
            my $ri = $c->{res}->[$ir];
            my $rnum = $ri->{num};
            if ( ( $rnum >= $LStart ) && ( $rnum <= $LStop ) ) {
              $c->{res}->[$ir]->{loop} = 1;
            }
          }
        }
      }
    }
    undef $lfile;
  }
  else {
    &GenUtil::log( "Warning: loop file " . $loopFname . " does not exist" );
    printf STDERR "Warning: loop file %s does not exist\n", $loopFname;
  }

  return;
}

## method: writeRegionFile()
## write a region file

sub writeRegionFile {
  my $self = shift;
  my $regionFname = shift;

  die "No chains defined for writeRegionFile"
    if ( !defined $self->activeChains() );

  $regionFname = $self->{par}->{protein} . ".regions"
    if ( !defined $regionFname );
  my $lfile = &GenUtil::getOutputFile($regionFname);

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
        if ( $c->{res}->[$ir]->{loop} == 1 ) {
          my $ri = $c->{res}->[$ir];
          my $Lstart = $ri->{num};
          my $Lend = $ri->{num};
          my $Llen = 1;
          $ir++;
          while ( ( defined $c->{res}->[$ir] )
            && ( $c->{res}->[$ir]->{loop} == 1 ) )
          {
            $Lend = $c->{res}->[$ir]->{num};
            $Llen++;
            $ir++;
          }
          printf $lfile "  %3d   %3d   %3d   %3d\n", $Lstart, $Lend, $Lstart,
            $Lend;
        }
      }
    }
  }
  undef $lfile;

  return;
}

## method: readRegionFile()
## read a region file

sub readRegionFile {
  my $self = shift;
  my $regionFname = shift;

  die "No chains defined for readRegionFile"
    if ( !defined $self->activeChains() );

  $regionFname = $self->{par}->{protein} . ".regions"
    if ( !defined $regionFname );
  my $lfile = &GenUtil::getInputFile($regionFname);
  while (<$lfile>) {
    chomp;
    my ( $LStart, $LStop, $junkA, $junkB ) = split( ' ', $_ );
    foreach my $c ( @{ $self->activeChains() } ) {
      if ( $#{ $c->{res} } >= 0 ) {
        for ( my $ir = 0 ; $ir <= $#{ $c->{res} } ; $ir++ ) {
          my $ri = $c->{res}->[$ir];
          my $rnum = $ri->{num};
          if ( ( $rnum >= $LStart ) && ( $rnum <= $LStop ) ) {
            $c->{res}->[$ir]->{loop} = 1;
          }
        }
      }
    }
  }
  undef $lfile;

  return;
}

## method: _readOutputDir(pathfile)
## Set the output directory by reading a pathfile

sub _readOutputDir {
  my $self = shift;
  my $pathfile = shift;

  my $inpath = &GenUtil::getInputFile($pathfile);

  # The 15th line contains the PDB output dir
  my $inline;
  for ( my $i = 1 ; $i <= 15 ; $i++ ) {
    $inline = <$inpath>;
  }
  chomp($inline);
  undef $inpath;
  my @p = split( ' ', $inline );

  return pop(@p);
}

## method: _processRosettaOutput()
## Move temporary files into the output dir,
## return an unsorted array of hashes with filenames

sub _processRosettaOutput {
  my $self = shift;

  &GenUtil::log("processing Rosetta output");

  # Move files listed in tempInfiles
  if ( &GenUtil::exists( $self->{par}->{tempInfiles} ) ) {
    my $tempInfiles = &GenUtil::getInputFile( $self->{par}->{tempInfiles} );
    while (<$tempInfiles>) {
      chomp;
      my $fname = $_;
      if ( &GenUtil::exists($fname) ) {
	# Move the expected file into the output directory
	&GenUtil::rename( $fname, $self->{par}->{outdir} );
      }

      # Also, try appending a zero. If this exists, keep incrementing the number
      # and moving over files until no additional files with these names exists.
      my $i=0;
      my $numbered_fname=$fname.$i;
      while ( &GenUtil::exists($numbered_fname) ) {
	&GenUtil::rename( $numbered_fname, $self->{par}->{outdir} );
	$i++;
	$numbered_fname=$fname.$i;
      }
    }

    undef $tempInfiles;
    &GenUtil::rename( $self->{par}->{tempInfiles}, $self->{par}->{outdir} );
  }
  &GenUtil::remove("checkpoint") if ( &GenUtil::exists("checkpoint") );

  # Also, move over any Rosetta MPI output files
  # (these are not listed in tempInfiles)
  my $mpinum=0;
  my $found=1;
  while ( $found ) {
    my $fname = "rosetta.mpi.out".$mpinum;
    if ( &GenUtil::exists( $fname ) ) {
      # Move the file into the output directory
      &GenUtil::rename( $fname, $self->{par}->{outdir} );
      $mpinum++;
    } else {
      $found=0;
    }
  }

  for ( my $i = 1 ; $i <= 6 ; $i++ ) {
    my $junkfile = "";
    $junkfile .= $self->{par}->{series} . $self->{par}->{protein};
    $junkfile .= "Aatom_set" . $i . ".pdb";
    &GenUtil::remove($junkfile);
  }
  my $junkfile = "";
  $junkfile .= $self->{par}->{series} . $self->{par}->{protein};
  $junkfile .= "A.chi_required.pdb";
  &GenUtil::remove($junkfile);

  $self->_collectHotspotOutnames() if ( $self->{par}->{incorporate_hotspot} );

  $self->_cleanOutput();

  my $ename = $self->{par}->{outdir} . $self->{par}->{outputPDBs};
  my $elist = &GenUtil::getInputFile($ename);
  my $output = ();
  while (<$elist>) {
    my $inline = $_;
    chomp($inline);
    my $f = {};
    ( $f->{filename}, $f->{Rosout}, $f->{origfile} ) =
      split( ' ', $inline );
    $f->{filename} = $f->{Rosout} if ( !$self->{par}->{doClean} );

    if ( $f->{filename} =~ /.gz$/ ) {
      push( @{$output}, $f )
	if ( &GenUtil::zexists( $f->{filename} ) &&
	     &GenUtil::validGzip( $f->{filename} ) );
    } else {
      push( @{$output}, $f )
	if ( &GenUtil::zexists( $f->{filename} ) );
    }
  }

  return $output;
}

## method: _collectHotspotOutnames()
## Since "incorporate_hotspot" generates an unpredictable number of output files
## (with an unpredictable number of names), update the "outputPDBs" file once
## Rosetta finishes.

sub _collectHotspotOutnames {
  my $self = shift;

  &GenUtil::log("updating hotspot output file names");

  my $ename = $self->{par}->{outdir} . $self->{par}->{outputPDBs};
  my $new_ename = $self->{par}->{outdir} . $self->{par}->{outputPDBs} . ".new";
  my $in_elist = &GenUtil::getInputFile($ename);
  my $out_elist = &GenUtil::getOutputFile($new_ename);
  while (<$in_elist>) {
    my $inline = $_;
    chomp($inline);
    my $f = {};
    ( $f->{filename}, $f->{Rosout}, $f->{origfile} ) =
      split( ' ', $inline );

    # Remove "_0001.pdb" from the filename, replace it with "_hs*_index*.des_*.pdb"
    # Note: leave a * at the end to match possbile .gz suffix
    my @new_fname_list = glob( substr( $f->{Rosout}, 0, -9) . "_hs*_index*.des_*.pdb*" );

    # List only the output files which conform to this new pattern in the output file list
    while ( $#new_fname_list >= 0 ) {
      my $outpdb = shift( @new_fname_list );
      if ( $outpdb !~ /premin/) {
	$outpdb = substr( $outpdb, 0, -3) if (&GenUtil::validGzip($outpdb));
	printf $out_elist "%s %s %s\n", $outpdb, $outpdb, $f->{origfile};
      }
    }
  }

  &GenUtil::rename( $new_ename, $ename );

  &GenUtil::log("done updating hotspot output file names");

  return;
}

## method: _cleanOutput()
## Create clean output files

sub _cleanOutput {
  my $self = shift;

  if ( $self->{par}->{doClean} ) {

    &GenUtil::log("cleaning Rosetta output files");

    my $ename = $self->{par}->{outdir} . $self->{par}->{outputPDBs};
    my $elist = &GenUtil::getInputFile($ename);
    while (<$elist>) {
      my $inline = $_;
      chomp($inline);
      my $f = {};
      ( $f->{filename}, $f->{Rosout}, $f->{origfile} ) =
        split( ' ', $inline );
      if ( &GenUtil::zexists( $f->{Rosout} )
        && &GenUtil::zexists( $f->{filename} )
        && &GenUtil::validGzip( $f->{filename} ) )
      {
        my $template = $self->{par}->{protein} . ".pdb";
        &_fitToTemplate( $template, $f );
      }
    }
  }

  &GenUtil::log("done cleaning Rosetta output files");

  return;
}

## function: _fitToTemplate()
## Match the input file to the template PDB based
## on residue numbering, NOT chain IDs. Write the output
## file using the coors from the input file and the chain
## IDs from the template. Residue names (and hence number
## of atoms) comes from the input coors, in case it's
## the output of design (ie. sequence mismatches).

sub _fitToTemplate {
  my $template = shift;
  my $f = shift;

  my $inPDB = $f->{Rosout};
  my $outPDB = $f->{filename};

  my $outMol = RosettaMolecule->new($template);

  my $infile = &GenUtil::getInputFile($inPDB);

  my $currChain = $outMol->{chain}->[0];
  my $currResInx = 0;
  my $currRes = $currChain->{res}->[0];

  # Clear all atom arrays and the atom indexes of residues
  foreach my $c ( @{ $outMol->activeChains() } ) {
    $c->{atom} = ();
    foreach my $r ( @{ $c->{res} } ) {
      $r->{start} = undef;
      $r->{end} = -1;
    }
  }

  while (<$infile>) {
    if (/^ATOM/) {
      my ( $atomname, $resname, $resnum, $iresnum, $alt, $chain, $seg );
      if ( ( $alt = substr( $_, 16, 1 ) ) =~ /[ A0-9]/ ) {
        ( $atomname = substr( $_, 12, 4 ) ) =~ s/ //g;
        $atomname .= $alt if ( $alt =~ /[0-9]/ );
        ( $resname = substr( $_, 17, 4 ) ) =~ s/ //g;
        $resnum = substr( $_, 22, 5 );
        ( $iresnum = $resnum ) =~ s/[A-Z]+//g;
        $iresnum += 0;
        $chain = substr( $_, 21, 1 );
        ( $seg = substr( $_, 72, 4 ) ) =~ s/[ \n]//g;

        # If the current atom doesn't belong here,
        # figure out where it belongs
        if ( $currRes->{num} != $iresnum ) {

          # First try the next residue, if there is one
          $currResInx++;
          if ( defined $currChain->{res}->[$currResInx] ) {
            $currRes = $currChain->{res}->[$currResInx];
          }

          # If this isn't the right place either,
          # start looking from scratch
          if ( $currRes->{num} != $iresnum ) {

            my $found = 0;
            foreach my $c ( @{ $outMol->activeChains() } ) {
              for (
                my $i = 0 ;
                ( ( $i <= $#{ $c->{res} } ) && ( !$found ) ) ;
                $i++
                )
              {
                if ( $c->{res}->[$i]->{num} == $iresnum ) {
                  $found = 1;
                  $currChain = $c;
                  $currResInx = $i;
                  $currRes = $currChain->{res}->[$i];
                }
              }
            }
            die "Cannot fit to template" if ( !$found );
          }
        }

        my $pdbrec = {};
        $pdbrec->{atominx} = substr( $_, 4, 7 ) + 0;
        $pdbrec->{atomname} = $atomname;
        $pdbrec->{resname} = $resname;
        $pdbrec->{resnum} = $iresnum;
        $pdbrec->{chain} = $currChain->{id};
        $pdbrec->{xcoor} = substr( $_, 30, 8 ) + 0.0;
        $pdbrec->{ycoor} = substr( $_, 38, 8 ) + 0.0;
        $pdbrec->{zcoor} = substr( $_, 46, 8 ) + 0.0;
        $pdbrec->{hyd} = ( $atomname =~ /^[0-9]*H.*/ ) ? 1 : 0;
        $pdbrec->{bb} = 1
          if ( ( $atomname eq "C" )
          || ( $atomname eq "O" )
          || ( $atomname eq "N" )
          || ( $atomname eq "H" )
          || ( $atomname eq "CA" ) );
        $pdbrec->{aux1} = substr( $_, 55, 6 ) + 0.0;
        $pdbrec->{aux2} = substr( $_, 61, 6 ) + 0.0;
        push( @{ $currChain->{atom} }, $pdbrec );

        my $ainx = $#{ $currChain->{atom} };
        $currRes->{start} = $ainx
          if ( !defined $currRes->{start} );
        $currRes->{end} = $ainx;
        $currRes->{name} = $resname;

      }
    }
  }

  undef $infile;
  $outMol->updateScore($inPDB);
  $outMol->writePDB( $outPDB, "full" );

  return;
}

# Return a random AA
sub genRandAA {

  my %AAlist = (
    1,  'GLY', 2,  'PRO', 3,  'ALA', 4,  'VAL', 5,  'LEU',
    6,  'ILE', 7,  'MET', 8,  'CYS', 9,  'PHE', 10, 'TYR',
    11, 'TRP', 12, 'HIS', 13, 'LYS', 14, 'ARG', 15, 'GLN',
    16, 'ASN', 17, 'GLU', 18, 'ASP', 19, 'SER', 20, 'THR'
  );

  return $AAlist{ int( rand(20) ) + 1 };
}

# Report the current value of a parameter
sub reportParameter {
  my $self = shift;
  my $k = shift;
  return $self->{par}->{$k};
}

sub setParameter {
  my $self = shift;
  $self->_getpar(@_);
  return;
}

sub _getpar {
  my $self = shift;
  my %arg = @_;
  foreach my $n ( keys %arg ) {
      $self->{par}->{$n} = $arg{$n};
      if ( ( $n eq "outdir" ) && ( chop( $arg{$n} ) ne "/" ) ) {
	  $self->{par}->{$n} .= "/";
	  &GenUtil::makeDir($self->{par}->{$n})
	      if ( !-d $self->{par}->{$n} );
      }
  }

  return;
}

1;

# INCOMPLETE OR DEPRECIATED (OUTDATED) SUBROUTINES

## method: robustDonateNatroResidues()
## Put the current structure's residue coordinates onto the
## coordinates of the RosettaMolecule passed in, for residues where
## self's des is set to NATRO. Note that missing residues in
## in either structure is not a problem

sub robustDonateNatroResidues {
  my $self = shift;
  my $inMol = shift;

  die "No chains defined for robustDonateNatroResidues"
    if ( !defined $self->activeChains() );

  foreach my $sc ( @{ $self->activeChains() } ) {
    my $dc = $inMol->getChain( $sc->{id} );
    die "Cannot lookup chain %s", $sc->{id}
      if ( !defined $dc->{res} );

    # Save the destination residues first, then
    # write over some of them with the source
    my %newcoors;
    foreach my $r ( @{ $dc->{res} } ) {
      my $trec = { res => $r, atom => $dc->{atom} };
      $newcoors{ $r->{num} } = $trec;
    }
    foreach my $r ( @{ $sc->{res} } ) {
      if ( $r->{des} eq "NATRO" ) {
        my $trec = { res => $r, atom => $sc->{atom} };
        $newcoors{ $r->{num} } = $trec;
      }
    }

    my $newatom = ();
    my $newres = ();
    my $ainx = 1;
    foreach my $k ( sort { $a <=> $b } keys %newcoors ) {
      my $r = $newcoors{$k}->{res};
      my $a = $newcoors{$k}->{atom};
      my $rrec = {};
      %{$rrec} = %{$r};
      push( @{$newres}, $rrec );
      $rrec->{start} = $#{$newatom} + 1;
      for ( my $ia = $r->{start} ; $ia <= $r->{end} ; $ia++ ) {
        my $arec = {};
        %{$arec} = %{ $a->[$ia] };
        $arec->{atominx} = $ainx++;
        push( @{$newatom}, $arec );
      }
      $rrec->{end} = $#{$newatom};
    }
    $dc->{atom} = $newatom;
    $dc->{res} = $newres;
  }

  $inMol->{segmentlist} = undef;
  $inMol->_coorCache();
  undef $inMol->{par}->{currFile};

  return;
}

## method: donateNatroResidues()
## Put the current structure's residue coordinates onto the
## coordinates of the RosettaMolecule passed in, for residues where
## self's des is set to NATRO. Note that the number of atoms
## in corresponding residues (and residue names) need not
## match, but ALL RESIDUES must exist in both cases

sub donateNatroResidues {
  my $self = shift;
  my $inMol = shift;

  die "No chains defined for donateNatroResidues"
    if ( !defined $self->activeChains() );

  foreach my $sc ( @{ $self->activeChains() } ) {
    my $dc = $inMol->getChain( $sc->{id} );
    die "Cannot lookup chain %s", $sc->{id}
      if ( !defined $dc->{res} );

    my $newatom = ();
    for ( my $nr = 0 ; $nr <= $#{ $sc->{res} } ; $nr++ ) {
      my $c = $dc;
      $c = $sc if ( $sc->{res}->[$nr]->{des} eq "NATRO" );
      my $aArr = $c->{atom};

      my $dres = $dc->{res}->[$nr];
      my $sres = $c->{res}->[$nr];
      my $start = $sres->{start};
      my $end = $sres->{end};
      $dres->{name} = $sres->{name};

      $dres->{start} = $#{$newatom} + 1;
      for ( my $ai = $start ; $ai <= $end ; $ai++ ) {
        push( @{$newatom}, $aArr->[$ai] );
      }
      $dres->{end} = $#{$newatom};
    }
    $dc->{atom} = $newatom;
  }
  $inMol->{segmentlist} = undef;
  $inMol->_coorCache();
  undef $inMol->{par}->{currFile};

  return;
}

## method: mergeChains()
## Incorporate residues from cB into cA. cA is kept in case of duplicates.

sub mergeChains {
  my $self = shift;
  my $cAid = shift;
  my $cBid = shift;

  my $cA = $self->getChain($cAid);
  die "Cannot find chain " . $cAid if ( !defined $cA->{res} );
  my $cB = $self->getChain($cBid);
  die "Cannot find chain " . $cBid if ( !defined $cB->{res} );

  # Save the destination residues first, then
  # write over some of them with the source
  my %newcoors;
  foreach my $r ( @{ $cB->{res} } ) {
    my $trec = { res => $r, atom => $cB->{atom} };
    $newcoors{ $r->{num} } = $trec;
  }
  foreach my $r ( @{ $cA->{res} } ) {
    my $trec = { res => $r, atom => $cA->{atom} };
    $newcoors{ $r->{num} } = $trec;
  }

  my $Acid = $cA->{id};
  my $newatom = ();
  my $newres = ();
  my $ainx = 1;
  foreach my $k ( sort { $a <=> $b } keys %newcoors ) {
    my $r = $newcoors{$k}->{res};
    my $a = $newcoors{$k}->{atom};
    my $rrec = {};
    %{$rrec} = %{$r};
    $rrec->{chain} = $Acid;
    push( @{$newres}, $rrec );
    $rrec->{start} = $#{$newatom} + 1;
    for ( my $ia = $r->{start} ; $ia <= $r->{end} ; $ia++ ) {
      my $arec = {};
      %{$arec} = %{ $a->[$ia] };
      $arec->{chain} = $Acid;
      $arec->{atominx} = $ainx++;
      push( @{$newatom}, $arec );
    }
    $rrec->{end} = $#{$newatom};
  }
  $cA->{atom} = $newatom;
  $cA->{res} = $newres;

  $self->removeChain($cBid);
  undef $self->{par}->{currFile};

  return;
}

## method: writeMutlistFile($WtMol,[parameters])
## Write a mutlist file

sub writeMutlistFile {
  my $self = shift;
  my $WtMol = shift;

  $self->_getpar(@_);

  my $mutLines = ();
  foreach my $dc ( @{ $self->activeChains() } ) {
    my $wc = $WtMol->getChain( $dc->{id} );
    if ( defined $wc ) {
      for ( my $ir = 0 ; $ir <= $#{ $dc->{res} } ; $ir++ ) {
        if ( $dc->{res}->[$ir]->{name} ne $wc->{res}->[$ir]->{name} ) {
          if ( !$dc->{res}->[$ir]->{loop} ) {
            my $desres = $dc->{res}->[$ir]->{name};
            my $WTres = $wc->{res}->[$ir]->{name};
            my $mutString = sprintf(
              "%13d  %s  %s  %s",
              $dc->{res}->[$ir]->{num},
              $dc->{id},
              $Sequence::_seqabbrev{$desres},
              $Sequence::_seqabbrev{$WTres}
            );
            push( @{$mutLines}, $mutString );
          }
        }
      }
    }
  }

  my $mutout = &GenUtil::getOutputFile( $self->{par}->{mutfile} );
  printf $mutout "START\n";
  printf $mutout "    1\n";
  printf $mutout "%3d\n", ( $#{$mutLines} ) + 1;
  foreach my $mut ( @{$mutLines} ) {
    printf $mutout "%s\n", $mut;
  }
  undef $mutout;

  return;
}

