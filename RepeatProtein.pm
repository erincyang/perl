# RepeatProtein package
# a derived class of RosettaMolecule,
# specific for systems in which one chain is a repeat protein
# (eg. ARs, TPRs, etc.)
#
# 2005, John Karanicolas, Baker group, UW
#
# derived from RosettaMolecule

package RepeatProtein;

require 5.004;

use strict;

use FileHandle;
use IPC::Open2;

use RosettaMolecule;
use Sequence;

use vars qw ( @ISA );

@ISA = ("RosettaMolecule");

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

## data: par -> { moduleContactRank, swapModule, etc.... }
## Extra Rosetta parameters specific to repeat proteins

## data: score -> { score, env, pair, vdw, hs, ss, sheet,
## data:            r-sigma, cb, rg, co, rama, hb_srbb,
## data:            hb_lrbb, etc.... }
## Rosetta score terms

## data: repeat -> { cid, positionOffset, module }
## specification of the repeat

## data: module[] -> { start, end, score, numC }
## specification of the repeat



## constructor: new([PDBfile],[parameters])
## creates a new RepeatProtein object and reads a PDB structures
## if a file name is given

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $farg=shift;
  my @pArr=@_;

  my $self=$class->SUPER::new($farg,@pArr);

  $self->{repeat}={};
  $self->{repeat}->{cid}=undef;
  $self->{repeat}->{module}=();
  $self->{repeat}->{positionOffset}=0;

  # Additional Rosetta parameter for repeat proteins
  $self->{par}->{moduleContactRank}=undef;
  $self->{par}->{swapModule}=undef;

  return $self;
}

# Destructor
DESTROY {
  my $self=shift;
  $self->SUPER::DESTROY();
  return;
}

## method: resetRepeat()
## reset all the information regarding repeats, EXCEPT the positionOffset

sub resetRepeat {
  my $self=shift;

  $self->{repeat}->{positionOffset}=0
      if (! defined $self->{repeat}->{positionOffset});
  $self->{repeat}->{cid}=undef;
  $self->{repeat}->{module}=();

  return;
}


## method: readPDB(file[,translate[,ignoreseg]])
## read an input PDB

sub readPDB {
  my $self=shift;
  my @args=@_;

  $self->SUPER::readPDB(@args);

  $self->resetRepeat();

  return;
}


## method: clone()
## clone the current molecule (ie. copy constructor)

sub clone {
  my $self = shift;
  my $valid = shift;

  my $n = $self->SUPER::clone($valid);

  # Clone the repeat hash
  $n->{repeat}={};
  $n->{repeat}->{cid}=$self->{repeat}->{cid};
  $n->{repeat}->{positionOffset}=$self->{repeat}->{positionOffset};
  $n->{repeat}->{module}=();
  if ($#{$self->{repeat}->{module}} >= 0) {
      foreach my $module ( @{$self->{repeat}->{module}} ) {
	  push(@{$n->{repeat}->{module}},$module);
      }
  }

  return $n;
}

## method: removeChain()
## polymorphism

sub removeChain {
  my $self = shift;
  my $cid = shift;

  $self->resetRepeat() if ($self->{repeat}->{cid} eq $cid);
  $self->SUPER::removeChain($cid);

  return;
}

## method: setChain([newcid])
## polymorphism

sub setChain {
  my $self = shift;
  my $newcid = shift;
  my @args = @_;

  my $c=$self->activeChains()->[0];
  $self->{repeat}->{cid}=$c->{id}
    if ($self->{repeat}->{cid} eq $c->{id});

  return $self->SUPER::setChain($newcid,@args);
}

## method: numModules()
## return the number of modules

sub numModules {
  my $self=shift;

  return 0 if (! defined $self->{repeat}->{module});
  return ($#{$self->{repeat}->{module}} + 1);
}

## method: defineRepeat(cid)
## define the bounds of each repeat

sub defineRepeat {
  my $self=shift;
  my $cid=shift;
  my $resBounds=shift;

  $self->resetRepeat();
  $self->defineRepeatChain($cid);
  $self->defineResBounds($resBounds);
  $self->setPositionOffset($self->{repeat}->{positionOffset});

  return;
}

## method: defineRepeatChain(cid)
## specify which chain is the repeat protein

sub defineRepeatChain {
  my $self=shift;
  my $cid=shift;

  $self->{repeat}->{cid}=$cid;
  return;
}

## method: getRepeatChain()
## return the chain containing the repeat protein

sub getRepeatChain {
  my $self=shift;

  return $self->getChain($self->{repeat}->{cid});
}

## method: getModule(moduleInx)
## return the module with the given moduleInx

sub getModule {
  my $self=shift;
  my $moduleInx=shift;

  return $self->{repeat}->{module}->[$moduleInx];
}

## method: defineResBounds(resBounds)
## define the residue indices which comprise each repeat

sub defineResBounds {
  my $self=shift;
  my $resBounds=shift;

  die "No resBounds defined"
      if ((! defined $resBounds) || ($#{$resBounds} < 0));

  my $cid = $self->{repeat}->{cid};
  die "Repeat chain must be defined before repeat bounds"
      if (! defined $cid);

  my $rc=$self->getRepeatChain();

  foreach my $module ( @{ $resBounds } ) {
      my $mstart = $module->{start};
      my $mend = $module->{end};
      push(@{$self->{repeat}->{module}},$module);
  }

  return;
}

## method: printRepeats( [outFname] )
## print the repeat definition (for verification)

sub printRepeats {
  my $self=shift;
  my $outFname=shift;

  return if ( ! defined $self->{repeat}->{cid});

  $outFname="-" if (! defined $outFname);
  my $outF=&GenUtil::getOutputFile($outFname);

  printf $outF "Repeat protein is chain %s\n",
    $self->{repeat}->{cid};

  my $c=$self->getRepeatChain();

  my $NcapEndInx;
  my $CcapStartInx;
  if ($#{$self->{repeat}->{module}} >= 0) {
      $NcapEndInx=$self->{repeat}->{module}->[0]->{start}-1;
      $CcapStartInx=$self->{repeat}->{module}->[$self->numModules()-1]->{end}+1;
  }

  if ((defined $NcapEndInx) && ($NcapEndInx >= 0)) {
      printf $outF "N-term cap: %d-%d\n",$c->{res}->[0]->{num},
        $c->{res}->[$NcapEndInx]->{num};
  }

  if ($self->numchains(1) > 1) {
      $self->computeRepeatContacts();
  }

  for (my $i=0; $i <= $#{$self->{repeat}->{module}}; $i++) {
      my $module=$self->getModule($i);
      my $startInx=$module->{start};
      my $endInx=$module->{end};
      printf $outF "Repeat %d: %d-%d",$i+1,
        $c->{res}->[$startInx]->{num}, $c->{res}->[$endInx]->{num};
      printf $outF " (module interface contacts: %d)",$module->{numC}
        if (defined $module->{numC});
      printf $outF " (module score: %f)",$module->{score}
        if (defined $module->{score});
      printf $outF "\n";
  }

  my $CcapEndInx=$self->lastResInx($c);
  if ((defined $CcapStartInx) && ($CcapStartInx <= $CcapEndInx)) {
      printf $outF "C-term cap: %d-%d\n",$c->{res}->[$CcapStartInx]->{num},
        $c->{res}->[$CcapEndInx]->{num};
  }

  undef $outF;

  return;
}

## method: detectRepeats()
## search the sequence of all chains for repeats

sub detectRepeats {
  my $self=shift;

  my $foundModule=0;

  my $cid;
  my $resBounds;

  # Detect across all chains, until we find a module
  foreach my $c ( @{$self->activeChains()} ) {
      if (! $foundModule) {
	  $cid=$c->{id};
	  if ($cid ne "+") {
	      ($foundModule,$resBounds)=$self->detectModuleInChain($cid);
	  }
      }
  }

  my $origPositionOffset=$self->{repeat}->{positionOffset};
  $self->{repeat}->{positionOffset}=0;

  $self->defineRepeat($cid,$resBounds) if ($foundModule);
  $self->setPositionOffset($origPositionOffset);

  return $foundModule;
}


## method: detectModuleInChain(chainId)
## search the sequence of the current chain for a repeat of the desired type
## (evaluation occurs in derived class, where >1 is a "hit")

sub detectModuleInChain {
  my $self=shift;
  my $cid=shift;

  my $foundModule=0;
  my $resBounds=();

  my $chain_score=0;

  my $c=$self->getChain($cid);
  my $lr=$self->lastResInx($c);

  my $scoreCut = 1.;

  for (my $inx=0; ($inx+$self->repeatLength()) < $lr; $inx++) {
      my $module_score=$self->evaluateModuleScore($c->{res},$inx);

      if ($module_score >= $scoreCut) {
	  my $module={};
	  $module->{start}=$inx;
	  $module->{end}=$inx+$self->repeatLength()-1;
	  $module->{score}=$module_score;
	  push(@{$resBounds},$module);
	  $chain_score+=$module_score;
      }

  }

  $foundModule=1 if ($chain_score > ($scoreCut * 2));

  # JK IF WE FOUND A MODULE,
  # GO BACK AND LOOK FOR ADJACENT MODULES USING A LOWER CUTOFF

  return ($foundModule,$resBounds);
}

## method: writePDBrepeat(outPDBname, exclude, pars)
## write a PDB file containing only the repeat protein
## (ie. throw out other chains)

sub writePDBrepeat {
  my $self=shift;
  my $outPDBname=shift;
  my $exclude=shift;
  my $mode=shift;
  my %wpar = @_;

  $exclude=0 if (! defined $exclude);

  # Note: we'll do a "clone" first, to avoid altering self
  my $mol=$self->clone(0);
  $mol->restrictToRepeat($exclude);
  $mol->writePDB($outPDBname,$mode,%wpar) if ($mol->numchains() > 0);
  undef $mol;

  return;
}

## method: writePDBmodule(outPDBname, moduleInx, pars)
## write a PDB file containing only the specified module

sub writePDBmodule {
  my $self=shift;
  my $outPDBname=shift;
  my $moduleInx=shift;
  my $mode=shift;
  my %wpar = @_;

  # Note: we'll do a "clone" first, to avoid altering self
  my $mol=$self->clone(0);
  $mol->restrictToModule($moduleInx);
  $mol->writePDB($outPDBname,$mode,%wpar) if ($mol->numchains() > 0);
  undef $mol;

  return;
}


## method: restrictToRepeat(exclude)
## write a PDB file containing only the repeat protein
## (ie. throw out other chains)

sub restrictToRepeat {
  my $self=shift;
  my $exclude=shift;

  $exclude=0 if (! defined $exclude);

  my $c=$self->getRepeatChain();
  if ($exclude) {
      $self->removeChain($c->{id});
  } else {
      $self->selectChain($c->{id});
      my $n=$self->clone(1);
      my $nc=$n->getChain();
      $c->{res}=$nc->{res};
      $c->{atom}=$nc->{atom};
  }

  return;
}


## method: restrictToModule(moduleInx)
## write a PDB file containing only the repeat protein
## (ie. throw out other chains)

sub restrictToModule {
  my $self=shift;
  my $moduleInx=shift;

  if (! defined $moduleInx) {
      return if ($self->numModules() > 1);
      $moduleInx=0;
  }

  $self->restrictToRepeat();
  my $c=$self->getRepeatChain();

  my $module=$self->getModule($moduleInx);
  my $startResInx=$module->{start};
  my $endResInx=$module->{end};

  # Set the module of interest to valid
  my $frag={};
  $frag->{chain}=$c->{id};
  $frag->{from}=$c->{res}->[$startResInx]->{num};
  $frag->{to}=$c->{res}->[$endResInx]->{num};
  my $farr=();
  push(@{$farr},$frag);
  $self->setValidResidues($farr,0);

  my $n=$self->clone(1);
  my $nc=$n->getChain();
  $c->{res}=$nc->{res};
  $c->{atom}=$nc->{atom};

  $self->resetRepeat();

  return;
}


## method: replaceModule(refMol, refDesc, moduleInx)
## replace a module in the current structure with the contents of refMol

sub replaceModule {
  my $self=shift;
  my $refMol=shift;
  my $moduleInx=shift;

  if (! defined $moduleInx) {
      return if ($self->numModules() > 1);
      $moduleInx=0;
  }

  die "refMol should only contain one chain when used for replacement"
      if ($refMol->numchains() > 1);

  my $c=$self->getRepeatChain();
  $refMol->setChain($c->{id});
  my $cRef=$refMol->getChain();

  # For now, require that the source module has the same number of residues as the destination module
  my $module=$self->getModule($moduleInx);
  my $moduleNumres=$module->{end}-$module->{start}+1;

  if ($refMol->numres() != $moduleNumres) {
      printf STDERR "refMol is of length %d\n",$refMol->numres();
      printf STDERR "Outgoing module is of length %d\n",$moduleNumres;
      die "Could not replace module due to length mismatch";
  }

  # Build just the outgoing module
  my $outmol=$self->clone(0);
  $outmol->restrictToModule($moduleInx);

  # Renumber and align refMol to match the outgoing module
  my $destStartnum=$c->{res}->[$module->{start}]->{num};
  my $destEndnum=$c->{res}->[$module->{end}]->{num};

  die "Outgoing module was not correctly built"
      if ($outmol->numres() != $moduleNumres);

  $refMol->renumber($destStartnum);
  $outmol->renumber($destStartnum);

  my $analyze=Analyze::new($outmol);
  $analyze->lsqfit($refMol,"ca",0,1);

  my $rmsd=$analyze->rmsd($refMol,0,undef,1);

  my $infoline="REPLACEMENT BACKBONE RMSD FOR RESIDUES ".$destStartnum."-".$destEndnum." WAS ".$rmsd->{"back"};
  $self->addInfoLine("MODULE",$infoline);

#  printf STDERR "CA rmsd for module swap is %f A\n",$rmsd->{"CA"};
#  printf STDERR "Backbone rmsd for module swap is %f A\n",$rmsd->{"back"};

  # Replace the atom and res arrays in self with those from refMol
  # First make a copy of the atom and res arrays in refMol
  # (otherwise we'll just transfer links, which will mess up refMol later)
  my $transferAtom=();
  foreach my $a ( @{$cRef->{atom}} ) {
      my $ta={};
      %{$ta}=%{$a};
      push(@{$transferAtom},$ta);
  }
  my $transferRes=();
  foreach my $r ( @{$cRef->{res}} ) {
      my $tr={};
      %{$tr}=%{$r};
      push(@{$transferRes},$tr);
  }

  # Now do the splice step
  my $startAtomInx=$c->{res}->[$module->{start}]->{start};
  my $endAtomInx=$c->{res}->[$module->{end}]->{end};
  my $origNumAtoms=$endAtomInx-$startAtomInx+1;
  splice(@{$c->{atom}},$startAtomInx,$origNumAtoms,@{$transferAtom});
  splice(@{$c->{res}},$module->{start},$moduleNumres,@{$transferRes});

  # At this point, residue numbers match but the atom indices to which they
  # refer will be wrong. We therefore need to traverse the res array of the
  # repeat chain and fix the "start" and "end" fields, based on what we find
  # in the atom array.
  my $resInx=0;
  $c->{res}->[$resInx]->{start}=0;
  my $lastResNum=$c->{atom}->[0]->{resnum};
  for (my $ainx=1; $ainx <= $#{$c->{atom}}; $ainx++) {
      $c->{atom}->[$ainx]->{atominx}=$ainx;
      if ($lastResNum != $c->{atom}->[$ainx]->{resnum}) {
	  $lastResNum = $c->{atom}->[$ainx]->{resnum};
	  $c->{res}->[$resInx]->{end}=$ainx-1;
	  $resInx++;
	  $c->{res}->[$resInx]->{start}=$ainx;
      }
  }
  $c->{res}->[$resInx]->{end}=$#{$c->{atom}};

  $self->_coorCache();
  undef $self->{par}->{currFile};

  return;
}


## method: transferModule(sourceMol, sourceModuleInx, destModuleInx)
## transfer a module from one repeat protein to
## the equivalent "frame" in the current repeat protein

sub transferModule {
  my $self=shift;
  my $sourceMol=shift;
  my $sourceModuleInx=shift;
  my $destModuleInx=shift;

  if (! defined $destModuleInx) {
      return if ($self->numModules() > 1);
      $destModuleInx=0;
  }
  if (! defined $sourceModuleInx) {
      return if ($sourceMol->numModules() > 1);
      $sourceModuleInx=0;
  }

  my $infoline="TRANSFERRING MODULE ".$sourceModuleInx+1;
  $self->addInfoLine("MODULE",$infoline);

  $sourceMol->restrictToModule($sourceModuleInx);
  $self->replaceModule($sourceMol,$destModuleInx);
  undef $sourceMol;

  return;
}

## method: designSwapModules([$startArr],[parameters])
## This polymorphic method is used instead of the one in RosettaMolecule
## if we're using "designSwapModules" method.
## This method is called by the "design" method

sub designSwapModules {
  my $self = shift;
  my $startArr = shift;
  $self->_getpar(@_);

  $self->{par}->{desSwapLoops}=1;
  my $desOut = $self->design($startArr);
  $self->{par}->{desSwapLoops}=0;

  return $desOut;
}

## method: _setupOutfiles($startArr)
## This polymorphic method is used instead of the one in RosettaMolecule
## if we're using "designSwapModules" method.
## This method is called by the "design" method

sub _setupOutfiles {
  my $self = shift;
  my $startArr = shift;
  my $noSeriesCode = shift;
  my $noSubdirs = shift;

  if ((! defined $self->{par}->{desSwapLoops}) ||
      ($self->{par}->{desSwapLoops} == 0)) {

      # Call regular subroutine in RosettaMolecule
      return $self->SUPER::_setupOutfiles($startArr,$noSeriesCode, $noSubdirs);

  } else {
      # Setup for designSwapLoops

      # Make sure required directory, files are present
      my $moduleDir = "alternate_modules/";
      die "Could not find directory".$moduleDir if (! -d $moduleDir);

      my @infoFiles = glob( $moduleDir."/*.info" );
      my $positionOffset=0;
      die "Not sure which info file to read" if ($#infoFiles > 0);
      if ($#infoFiles == 0) {
	  my $infoInfile=shift(@infoFiles);
	  my $inInfo = &GenUtil::getInputFile($infoInfile);
	  while (<$inInfo>) {
	      my $inline = $_;
	      chomp($inline);
	      if (uc($inline) =~ /^MODULE TYPE ([A-Za-z]+)/) {
		  if (uc($self->moduleName()) ne uc($1)) {
		      printf STDERR "Comparing module type %s",uc($self->moduleName());
		      printf STDERR " to module type %s\n",uc($1);
		      die "Wrong module type";
		  }
	      } elsif (uc($inline) =~ /^POSITION OFFSET ([0-9]+)/) {
		  $positionOffset=$1;
	      }
	  }
	  undef $inInfo;
	  &GenUtil::log("done reading module info file, using positionOffset=".$positionOffset);
      } else {
	  &GenUtil::log("warning: no module info file found");
	  &GenUtil::log("modules will be assumed to match");
	  &GenUtil::log("assuming positionOffset=0");
      }
      $self->setPositionOffset($positionOffset);

      # Setup an array of Molecules corresponding to the incoming modules
      my @moduleFiles = glob( $moduleDir."/*" );
      my $moduleArr=();
      foreach my $mfname ( @moduleFiles ) {
	  if (&PDButil::isaPDB($mfname)) {
	      my $newModule=Molecule->new($mfname);
	      # A slight departure from object-oriented coding below, but we'll allow it!
	      $newModule->{source}=$mfname;
	      push(@{$moduleArr},$newModule);
	  }
      }
      die "No alternate module PDBs found" if ($#{$moduleArr} < 0);

      # Write the output files, listing them for Rosetta as we go
      my $filelist = "infiles.list";
      my $flist = &GenUtil::getOutputFile($filelist);
      my $ename = $self->{par}->{outdir} . $self->{par}->{outputPDBs};
      my $elist = &GenUtil::getOutputFile($ename);
      my $earr = ();
      my $fnum = 1;
      my $numOutfiles = 0;

      # Log the method of swapping
      if (defined $self->{par}->{swapModule}) {
	  &GenUtil::log("swapping predefined module ".
			$self->{par}->{swapModule});
      } elsif  (defined $self->{par}->{moduleContactRank}) {
	  &GenUtil::log("swapping module based on ContactRank ".
			$self->{par}->{moduleContactRank});
      } else {
	  &GenUtil::log("warning: swapModule not found, moduleContactRank not found");
	  &GenUtil::log("replacing first module by default");
      }

      foreach my $startfile ( @{$startArr} ) {

	  # Build startMol via cloning to ensure the object type is correct
	  my $startMol=$self->clone(1);
	  $startMol->readPDB($startfile->{filename});

	  my $replaceModuleInx;
	  if (defined $startMol->{par}->{swapModule}) {
	      $replaceModuleInx=$startMol->{par}->{swapModule};
	  } elsif  (defined $startMol->{par}->{moduleContactRank}) {
	      # Compute number of interface contacts for each module
	      $startMol->computeRepeatContacts();
	      my $sMod=();
	      for (my $i=0; $i<$startMol->numModules(); $i++) {
		  my $f={};
		  my $module=$startMol->getModule($i);
		  $f->{mInx}=$i;
		  $f->{numC}=$module->{numC};
		  push(@{$sMod},$f);
	      }

	      # Sort in descending order of interface contacts
	      $sMod=[sort( {$b->{numC} <=> $a->{numC}} @{$sMod})];

	      # Pick the desired module
	      $replaceModuleInx=$sMod->[$startMol->{par}->{moduleContactRank}-1]->{mInx};

	  } else {
	      $replaceModuleInx=0;
	  }

	  foreach my $newModule ( @{$moduleArr} ) {

	      my $outdir = ($fnum % 100) . "/";
	      my $fulldir = $self->{par}->{outdir}.$outdir;
	      &GenUtil::makeDir($fulldir) if ( !-d $fulldir );

	      my $newMol=$startMol->clone(0);
	      $newMol->addInfoLine("MODULE SOURCE",$newModule->{source});
	      $newMol->replaceModule($newModule,$replaceModuleInx);
	      my $fname = "temp" . $fnum . ".pdb";
	      $newMol->writePDB($fname);

	      printf $flist "%s %s\n", $fname, $outdir.$fname;
	      my $efile;
	      for ( my $i = 1 ; $i <= $self->{par}->{nstruct} ; $i++ ) {
		  $efile = $fulldir . "temp" . $fnum . "_" . &GenUtil::zPad( $i, 4 ) . ".pdb";
		  my $cfile = $fulldir . "temp";
		  $cfile .= $fnum . "_" . &GenUtil::zPad( $i, 4 ) . ".clean.pdb";
		  printf $elist "%s %s %s %s\n", $cfile, $efile, $startfile->{filename};
		  $numOutfiles++;
	      }
	      push( @{$earr}, $efile );
	      $fnum++;
	  }
      }

      undef $flist;
      undef $elist;

      return ($filelist,$earr,$numOutfiles);

  }

  return;
}

## method: _processRosettaOutput()
## This polymorphic method is used instead of the one in RosettaMolecule
## if we're using "designSwapModules" method.
## This method is called by the "design" method

sub _processRosettaOutput {
  my $self = shift;

  if ((defined $self->{par}->{desSwapLoops}) &&
      ($self->{par}->{desSwapLoops} == 1)) {

      # We created "module-swapped" input files - remove these
      my $fnum = 1;
      my $fname = "temp" . $fnum . ".pdb";
      while (-e $fname) {
	  &GenUtil::remove($fname);
	  $fnum++;
	  $fname = "temp" . $fnum . ".pdb";
      }
  }

  # Finish by calling regular subroutine in RosettaMolecule
  return $self->SUPER::_processRosettaOutput();
}

## method: computeRepeatContacts()
## compute the number of C-beta - C-beta contacts between
## each module and other protein chains

sub computeRepeatContacts {
  my $self=shift;

  for (my $i=0; $i<$self->numModules(); $i++) {
      $self->computeModuleContacts($i);
  }

  return;
}

## method: computeModuleContacts(moduleInx)
## compute the number of C-beta - C-beta contacts between
## this module and other protein chains

sub computeModuleContacts {
  my $self=shift;
  my $moduleInx=shift;

  my $sqDistCut=49.;

  if (! defined $moduleInx) {
      if ($self->numModules() == 1) {
	  $moduleInx=0;
      } else {
	  return;
      }
  }

  my $rc=$self->getRepeatChain();
  my $module=$self->getModule($moduleInx);
  my $startResInx=$module->{start};
  my $endResInx=$module->{end};
  my $startAtomInx=$rc->{res}->[$startResInx]->{start};
  my $endAtomInx=$rc->{res}->[$endResInx]->{end};

  # Collect the Cbetas in the module first
  my $rCB=();
  for (my $ainx=$startAtomInx; $ainx<=$endAtomInx; $ainx++) {
      my $ra=$rc->{atom}->[$ainx];
      if (($ra->{atomname} eq "CB") ||
	  (($ra->{atomname} eq "CA") && ($ra->{resname} eq "GLY"))) {
	  push(@{$rCB},$ra);
      }
  }

  my $moduleContacts=0;
  foreach my $c ( @{ $self->activeChains() } ) {
      if (($c->{id} ne $rc->{id}) && ($#{$c->{atom}} >= 0)) {
	  foreach my $a ( @{$c->{atom}} ) {
	      if (($a->{atomname} eq "CB") ||
		  (($a->{atomname} eq "CA") && ($a->{resname} eq "GLY"))) {
		   my $x1=$a->{xcoor};
		   my $y1=$a->{ycoor};
		   my $z1=$a->{zcoor};

		   foreach my $ra ( @{$rCB} ) {
		       my $x2=$ra->{xcoor};
		       my $y2=$ra->{ycoor};
		       my $z2=$ra->{zcoor};
		       my $SqDist=($x1-$x2)*($x1-$x2);
		       $SqDist+=($y1-$y2)*($y1-$y2);
		       $SqDist+=($z1-$z2)*($z1-$z2);
		       $moduleContacts++ if ($SqDist < $sqDistCut);
		   }

	      }
	  }
      }
  }
  $module->{numC}=$moduleContacts;

  return;
}

## method: setResfileMSA()
## set the {des} fields given the MSA info (allowed AAs at each position)
## Note: we'll do this by disallowing all the other AAs at this position
## Note: we'll also set "NonBinding" residues to NATRO
##    (residues known not to be in the interface)

sub setResfileMSA {
  my $self=shift;

  return if ( $#{$self->{repeat}->{module}} < 0 );

  my $c=$self->getRepeatChain();
  for (my $i=0; $i <= $#{$self->{repeat}->{module}}; $i++) {
      my $module=$self->getModule($i);
      my $startInx=$module->{start};
      foreach my $pos ( @{$self->listWellConserved()} ) {
	  my $r=$c->{res}->[$startInx+$pos->{offset}];
	  $r->{des}="NOTAA  ".&Sequence::notRes($pos->{id});
      }
      foreach my $pos ( @{$self->listNonBinding()} ) {
	  my $r=$c->{res}->[$startInx+$pos];
	  $r->{des}="NATRO";
      }
  }

  return;
}

## method: printWellConservedPositions( [outFname] )
## print the residues which occupy the well-conserved positions
## (most popular type in parentheses)

sub printWellConservedPositions {
  my $self=shift;
  my $outFname=shift;

  return if ( $#{$self->{repeat}->{module}} < 0 );

  $outFname="-" if (! defined $outFname);
  my $outF=&GenUtil::getOutputFile($outFname);
  printf $outF "\nWell-conserved positions:\n";

  my $c=$self->getRepeatChain();
  for (my $i=0; $i <= $#{$self->{repeat}->{module}}; $i++) {
      my $module=$self->getModule($i);
      my $startInx=$module->{start};
      my @conslist=();
      foreach my $pos ( @{$self->listWellConserved()} ) {
	  my $r=$c->{res}->[$startInx+$pos->{offset}];
	  my $constr=$Sequence::_seqabbrev{$r->{name}}.$r->{num}." (".$pos->{id}.")";
	  push(@conslist,$constr);
      }
      printf $outF "Repeat %d: %s\n",$i+1,
        join(", ",@conslist);
  }

  printf $outF "\n";
  undef $outF;

  return;
}

## method: printSemiConservedPositions( [outFname] )
## print the residues which occupy the semi-conserved positions
## (most popular type in parentheses)

sub printSemiConservedPositions {
  my $self=shift;
  my $outFname=shift;

  return if ( $#{$self->{repeat}->{module}} < 0 );

  $outFname="-" if (! defined $outFname);
  my $outF=&GenUtil::getOutputFile($outFname);
  printf $outF "\nSemi-conserved positions:\n";

  my $c=$self->getRepeatChain();
  for (my $i=0; $i <= $#{$self->{repeat}->{module}}; $i++) {
      my $module=$self->getModule($i);
      my $startInx=$module->{start};
      my @conslist=();
      foreach my $pos ( @{$self->listSemiConserved()} ) {
	  my $r=$c->{res}->[$startInx+$pos->{offset}];
	  my $constr=$Sequence::_seqabbrev{$r->{name}}.$r->{num}." (".$pos->{id}.")";
	  push(@conslist,$constr);
      }
      printf $outF "Repeat %d: %s\n",$i+1,
        join(", ",@conslist);
  }

  printf $outF "\n";
  undef $outF;

  return;
}

## method: setPositionOffset()
## specify the position offset

sub setPositionOffset {
  my $self=shift;
  my $positionOffset=shift;

  $positionOffset=0 if (! defined $positionOffset);

  return if ($positionOffset == $self->{repeat}->{positionOffset});

  my $posShift=$positionOffset-$self->{repeat}->{positionOffset};

  my $c=$self->getRepeatChain();
  for (my $i=0; $i <= $#{$self->{repeat}->{module}}; $i++) {
      my $module=$self->getModule($i);
      $module->{start} += $posShift;
      $module->{end} += $posShift;
      if (($module->{start} < 0) || ($module->{end} > $self->lastResInx($c))) {
	  # Remove this module, as it exceeds defined residue bounds
	  splice(@{$self->{repeat}->{module}},$i,1);
	  $i--;
      }
  }

  $self->{repeat}->{positionOffset}=$positionOffset;

  return;
}

## method: reportPositionOffset()
## report the current position offset

sub reportPositionOffset {
  my $self=shift;

  return $self->{repeat}->{positionOffset};
}

## method: moduleName()
## return the name of this type of module

sub moduleName {
  my $self=shift;

  return "Non-repeating";
}

## method: repeating()
## if object is of this class, there may not be repeating units
## (derived classes have repeating units, this is just the base class)

sub repeating {
  my $self=shift;

  return 0;
}

## method: repeatLength()
## repeat length should be defined in a derived class (at least for now)

sub repeatLength {
  my $self=shift;

  die "repeat length should be defined in a derived class";
  return;
}


## method: evaluateModuleScore(resArr,startInx)
## rules for evaluating a module should be defined in a derived class

sub evaluateModuleScore {
  my $self=shift;

  die "evaluation of modules should be defined in a derived class";
  return;
}

## method: listWellConserved()
## list of well-conserved residues should be defined in a derived class

sub listWellConserved {
  my $self=shift;

  die "list of well-conserved residues should be defined in a derived class";
  return;
}

## method: listSemiConserved()
## list of semi-conserved residues should be defined in a derived class

sub listSemiConserved {
  my $self=shift;

  die "list of semi-conserved residues should be defined in a derived class";
  return;
}

