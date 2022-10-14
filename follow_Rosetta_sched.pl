#!/usr/bin/env perl

# Follow the schedule in the schedule file
# Note: chains A-M are designed and flexible, chains N-Z are fixed

sub usage {
  printf STDERR "usage:   follow_Rosetta_sched.pl [inPDB]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;
use GenUtil;
use KnownRepeatProteins;
use UWutil;
use PDButil;
use lib ".";
use autouse Custom => qw(customSetup);

# So I can run in the background
close(STDIN);

my $inPDB;
while ($#ARGV>=0) {
    if (lc($ARGV[0]) eq "-help" || lc($ARGV[0]) eq "-h") {
	&usage();
    } elsif ($ARGV[0] =~ /^-/) {
	printf STDERR "unknown option %s\n",shift @ARGV;
	&usage();
    } else {
	$inPDB=shift @ARGV;
    }
}

# The log file
&GenUtil::setLogFile("flex_des.log");
&GenUtil::log("starting flex_des");

# Hold if this is a restart and jobs are still running
&UWutil::waitForClusterCompletion();

# Load custom functions if available
my $customFunc={};
if (&GenUtil::zexists("./Custom.pm")) {
    &GenUtil::log("loading custom functions");
    &customSetup($customFunc);
} else {
    &GenUtil::log("no custom functions will be used");
}

# If there's one suitable PDB file, use that. Otherwise, if
# there's one suitable sched file, use that. Otherwise die.
if (! defined $inPDB) {
    my @pdbs=glob("????.pdb ????.pdb.gz");
    if ($#pdbs == 0) {
	$inPDB=$pdbs[0];
    } else {
	my @scheds=glob("????.sched");
        if ($#scheds == 0) {
	    $inPDB=substr($scheds[0],0,-6).".pdb";
	    die "Cannot find PDBfile for existing schedule file"
		if (! &GenUtil::zexists($inPDB));
	} else {
	    die "Need to specify inPDB";
	}
    }
} elsif (length($inPDB) != 8) {
    die "Please use a PDB file with a 4 character ID";
}

# Read in the structure
&GenUtil::log("setting up system from: ".$inPDB);
my $Rmol=KnownRepeatProteins::newProtein($inPDB);
my $prefix=&PDButil::prefix($inPDB);

# Setup the schedule
my $sched;
my $schedFile=$prefix.".sched";

if (&GenUtil::zexists($schedFile)) {
    $sched=&readSchedFile($schedFile);
} else {
    die "Cannot read ".$schedFile." - no schedule found";
}

# Set firstrun based on previous results
my $firstrun;
my $bestFile="best.01.pdb";
for ($firstrun=1; &GenUtil::zexists($bestFile);) {
    $firstrun++;
    $bestFile="best.".&GenUtil::zPad($firstrun,2).".pdb";
}

my $sArr;

# Update sArr from the latest "scores" file
if ($firstrun>1) {
    $sArr=&restart_from_scores("scores.".&GenUtil::zPad($firstrun-1,2).".out");
} else {
    $sArr=&restart_from_PDBs($inPDB);
}

for (my $runnum=$firstrun; ($runnum-1) <= $#{$sched}; $runnum++) {

     # Re-read the schedule file at every step, in case it has changed
     $sched=&readSchedFile($schedFile);
     my $step=$sched->[$runnum-1];

     &GenUtil::log("current cycle: ".$runnum);

     if ($runnum>1) {
       $Rmol->readPDB("best.".&GenUtil::zPad($runnum-1,2).".pdb");
     }

     # Reset parameters, preserving protein ID and currFile status
     my $protein=$Rmol->reportParameter("protein");
     my $currFile=$Rmol->reportParameter("currFile");
     $Rmol->resetPars(protein=>$protein,currFile=>$currFile);

     # Always use_input_sc and fake_native and jk_interface
     $Rmol->setParameter(use_input_sc=>1);
     $Rmol->setParameter(fake_native=>1);
     $Rmol->setParameter(jk_interface=>1);

     # If an appropriate resfile exists, point to this (always)
     my $resfile=$prefix.".resfile";
     if (&GenUtil::zexists($resfile)) {
       $Rmol->setParameter(resfile=>$resfile);
       $Rmol->readResfile();
     }

     # Turn on jk_interface flag
#     $Rmol->setParameter(jk_interface=>1);

     # The two-letter code
     my $tag=&GenUtil::zPad($runnum,2);
     $Rmol->setParameter(series=>$tag);

     # Process the instruction
     my $num_structures = 1;
     while ( ($num_structures > 0 ) && ( $#{$step}>=0) ) {
	 $sArr=&processInstruction($Rmol, shift(@{$step}), $sArr);
	 $num_structures = $#{$sArr}+1;
	 &GenUtil::log("current number of structures: ".$num_structures);
     }

     # Save the best one
     my $bestFile="best.".$tag.".pdb";
     my $newbest;
     if ( $num_structures > 0 ) {
       $newbest=$sArr->[0]->{filename};
     } else {
       $newbest="best.".&GenUtil::zPad($runnum-1,2).".pdb";
     }
     &GenUtil::copy($newbest,$bestFile);

     # Write the scores file
     my $scoreFile="scores.".$tag.".out";
     &write_scores($scoreFile,$sArr);

     # Write the origins file
     my $originFile="origins.".$tag.".out";
     &write_origins($originFile);

     # Tar results that are no longer needed, remove the directory
#     if ($runnum > 1) {
#       my $tartag=&GenUtil::zPad(($firstrun-1),2);
#       my @inlist=glob("*");
#       foreach my $flist ( @inlist ) {
#	 if ( (-d $flist) && ($flist =~ /_$tartag/) ) {
#	   if (! -e ($flist.".tar.gz")) {
#	     my $tarcmd="/bin/tar cfz old.tar.gz ".$flist;
#	     &GenUtil::log("tarring files from ".$flist);
#	     system($tarcmd);
#	     &GenUtil::rename("old.tar.gz",$flist.".tar.gz");
#	   }
#	   &GenUtil::remove($flist)
#	     if (&GenUtil::zexists($flist.".tar.gz"));
#	 }
#       }
#     }

     if ($#{$sArr} < 0) {
       &GenUtil::log("No structures remain.");
       die "No structures remain to process instruction";
     }

}

my $num_structures = $#{$sArr}+1;
if ( $num_structures > 0 ) {
  printf STDERR "Schedule complete, with %d output structures!\n", $num_structures;
} else {
  printf STDERR "Schedule complete, with no output structures!\n", $num_structures;
}

&GenUtil::log("EXIT - NORMAL TERMINATION");
exit(0);



#########################################################


sub readSchedFile {
    my $schedFile=shift;

    &GenUtil::log("reading schedule file");

    my $sched=();
    my $infile=&GenUtil::getInputFile($schedFile);
    while (<$infile>) {
	my $inline=$_;
	chomp($inline);
	# Remove comments
	($inline=$inline) =~ s/\#.*//g;
	my @p=split(' ',$inline);
	my $steps=();
        while ($#p>=0) {
	    my $inline=shift(@p);
	    # Remove leading and trailing whitespace
	    ($inline=$inline) =~ s/^ *//g;
	    ($inline=$inline) =~ s/ *$//g;
	    push(@{$steps},$inline)
		if ($inline ne "");
	}
        push(@{$sched},$steps) if ($#{$steps} >= 0);
    }

    undef $infile;
    return $sched;
}

sub restart_from_PDBs {
    my $inPDB=shift;

    # Don't start from any that have the current prefix in the filename
    # Only start from those that have "temp" somewhere in the filename
    my $sArr=();
    my $prefix=&PDButil::prefix($inPDB);
    my @pdbs=glob("*temp*");
    foreach my $pdbfile ( @pdbs ) {
      if (($pdbfile !~ $prefix) && &PDButil::isaPDB($pdbfile)) {
	my $r={};
	$r->{filename}=$pdbfile;
	push(@{$sArr},$r);
      }
    }

    # If there are none of these, start from those that have "start"
    # somewhere in the filename, unzipping them if they're zipped
    if ($#{$sArr} < 0) {
      my @pdbs=glob("*start*");
      foreach my $pdbfile ( @pdbs ) {
	if (&PDButil::isaPDB($pdbfile)) {
	  if (&GenUtil::validGzip($pdbfile)) {
	    system("gunzip ".$pdbfile);
	    $pdbfile=substr($pdbfile,0,-3);
	  }
	  my $r={};
	  $r->{filename}=$pdbfile;
	  push(@{$sArr},$r);
	}
      }
    }

    return $sArr;
}

sub restart_from_scores {
    my $fname=shift;

    my $sArr=();
    my $infile=&GenUtil::getInputFile($fname);

    # The first line is the header
    my $first_line=<$infile>;
    my @hinfo=split(' ',$first_line);

    # Read subsequent lines
    while(<$infile>) {
      my $inline=$_;
      chomp($inline);
      my @f=split(' ',$inline);
      my $r={};
      for (my $i=0; $i<=$#f; $i++) {
	$r->{$hinfo[$i]}=$f[$i];
      }
      push(@{$sArr},$r);

    }
    undef $infile;

    if ($#{$sArr} < 0) {
      &GenUtil::log("No structures remain.");
      printf STDERR "No structures remain to process instruction\n";
      exit(0);
    }

    return $sArr;
}

sub write_scores {
    my $fname=shift;
    my $sArr=shift;

    my $outfile=&GenUtil::getOutputFile($fname);

    # Write the header (based on fields present in the first sArr element)
    my $ele=$sArr->[0];
    my @flist=();
    printf $outfile "%s","filename";
    foreach my $k ( keys %{$ele} ) {
	if (($k ne "filename") && ($k ne "Rosout") &&
            ($k ne "Rosinp") && ($k ne "origfile")) {
	    push(@flist,$k);
	    printf $outfile "  %s",$k;
	}
    }
    printf $outfile "\n";

    # Write info for each sArr element
    for (my $i=0; $i<=$#{$sArr}; $i++) {
	printf $outfile "%s",$sArr->[$i]->{"filename"};
        for (my $j=0; $j<=$#flist; $j++) {
	     printf $outfile "  %f",$sArr->[$i]->{$flist[$j]};
	}
	printf $outfile "\n";
    }
    undef $outfile;

    return;
}

sub write_origins {
    my $fname=shift;

    my $outfile=&GenUtil::getOutputFile($fname);

    # Read the list of output PDBs for the most recent step, write a stripped-down
    # version of this
    my $infname=$Rmol->reportParameter("outdir");
    $infname.=$Rmol->reportParameter("outputPDBs");
    my $infile=&GenUtil::getInputFile($infname);
    while (<$infile>) {
	my $inline=$_;
	chomp($inline);
	my @p=split(' ',$inline);
	my $child=$p[1];
	my $parent=$p[2];
	printf $outfile "%s %s\n", $child, $parent;
    }
    undef $outfile;
    undef $infile;

    return;
}

sub setRosettaParams {
    my $Rmol=shift;
    my $args=shift;

    foreach my $k ( keys %{$args} ) {
	$Rmol->setParameter($k=>$args->{$k});
    }
#    $Rmol->setParameter(interface_ds=>1);

    my $gufile=$prefix.".gu";
    $Rmol->setParameter(gu_exempt=>$gufile)
	if (&GenUtil::exists($gufile));

    return;
}

sub processInstruction {
    my $Rmol=shift;
    my $instruction=shift;
    my $sArr=shift;

    if ($#{$sArr} < 0) {
      &GenUtil::log("No structures remain.");
      printf STDERR "No structures remain to process instruction\n";
      exit(0);
    }

    &GenUtil::log("Processing instruction ".$instruction);
    my @p=split(/:/,$instruction);
    my $cmd=shift(@p);
    my $args={};
    while ($#p>=0) {
	my ($key,$val)=split(/=/,shift(@p));
	$args->{$key}=$val;
    }
    &setRosettaParams($Rmol,$args);

    if ($cmd eq "curr") {

	my $outdir="curr_".$Rmol->reportParameter("series");
	$Rmol->setParameter(outdir=>$outdir);
	&GenUtil::log("starting curr step");

	&GenUtil::log("storing existing eval fields");
	my $prevEval={};
	foreach my $f ( @{$sArr}) {
	  my $fname=$f->{filename};
	  foreach my $k ( keys %{$f} ) {
	    $prevEval->{$fname}->{$k}=$f->{$k}
	      if ( ( $k ne "filename" ) && ( $k ne "Rosout" ) && ( $k ne "origfile" ) );
	  }
	}

	$sArr=$Rmol->null_call($sArr);

	&GenUtil::log("recovering existing eval fields");
	foreach my $f ( @{$sArr}) {
	  my $oldFname=$f->{origfile};
	  foreach my $k ( keys %{$prevEval->{$oldFname}} ) {
	    $f->{$k}=$prevEval->{$oldFname}->{$k};
	  }
	}

	&GenUtil::log("finished curr step");

    } elsif ($cmd eq "frag") {

	my $outdir="frag_".$Rmol->reportParameter("series");
	$Rmol->setParameter(outdir=>$outdir);
	&GenUtil::log("starting fragment insertions");
	$Rmol->readLoopFile();
	$Rmol->setParameter(nstruct=>$args->{tot}) if (defined $args->{tot});
	$sArr=$Rmol->abinitFrag($sArr);
	&GenUtil::log("finished fragment insertions");

    } elsif ($cmd eq "patchdock") {

	my $outdir="patchdock_".$Rmol->reportParameter("series");
	$Rmol->setParameter(dockConsFname=>($prefix.".cons"))
	  if (&GenUtil::exists($prefix.".cons"));
	$Rmol->setParameter(outdir=>$outdir);
	&GenUtil::log("starting patchdock");
	$sArr=$Rmol->patchdock();
	$Rmol->resetDockFlags();
	&GenUtil::log("finished patchdock");

    } elsif ($cmd eq "dockrand") {

	my $outdir="dockrand_".$Rmol->reportParameter("series");
	$Rmol->setParameter(outdir=>$outdir);
	&GenUtil::log("starting random docking");
	$Rmol->setParameter(decoystats=>0);
	$Rmol->setParameter(dockConsFname=>($prefix.".cons"))
	    if (&GenUtil::exists($prefix.".cons"));
	$sArr=$Rmol->dock($sArr,dockRand=>1);
	$Rmol->resetDockFlags();
	$Rmol->setParameter(decoystats=>1);
	&GenUtil::log("finished random docking");

    } elsif ($cmd eq "dockpert") {

	my $outdir="dockpert_".$Rmol->reportParameter("series");
	$Rmol->setParameter(outdir=>$outdir);

	$Rmol->setParameter(nstruct=>$args->{explode}) if (defined $args->{explode});
	if (defined $args->{tot}) {
	    if ($args->{tot} > ($#{$sArr}+1)) {
 	        $Rmol->setParameter(nstruct=>int($args->{tot}/($#{$sArr}+1)));
	    } else {
 	        $Rmol->setParameter(nstruct=>1);
	    }
	}

	&GenUtil::log("starting docking perturbations to generate diversity");
	$Rmol->setParameter(decoystats=>0);
#	$Rmol->setParameter(dockGenerateDiversity=>1);
	$Rmol->setParameter(dockConsFname=>($prefix.".cons"))
	    if (&GenUtil::exists($prefix.".cons"));
	$sArr=$Rmol->dock($sArr,dockPert=>"2 2 2");
	$Rmol->resetDockFlags();
	$Rmol->setParameter(decoystats=>1);
	&GenUtil::log("finished docking perturbations");

    } elsif ($cmd eq "dockmcm") {

	my $outdir="dockmcm_".$Rmol->reportParameter("series");
	$Rmol->setParameter(outdir=>$outdir);
	$Rmol->setParameter(interface_ds=>1);
	$Rmol->setParameter(nstruct=>$args->{explode}) if (defined $args->{explode});
	if (defined $args->{tot}) {
	    if ($args->{tot} > ($#{$sArr}+1)) {
 	        $Rmol->setParameter(nstruct=>int($args->{tot}/($#{$sArr}+1)));
	    } else {
 	        $Rmol->setParameter(nstruct=>1);
	    }
        }

	&GenUtil::log("starting docking mcm");
	$Rmol->setParameter(dockConsFname=>($prefix.".cons"))
	    if (&GenUtil::exists($prefix.".cons"));
	$sArr=$Rmol->dock($sArr,dockMcm=>1);
	$Rmol->resetDockFlags();
	&GenUtil::log("finished docking mcm");

    } elsif ($cmd eq "dockmin") {

	my $outdir="dockmin_".$Rmol->reportParameter("series");
	$Rmol->setParameter(outdir=>$outdir);
	$Rmol->setParameter(nstruct=>$args->{explode}) if (defined $args->{explode});
	&GenUtil::log("starting docking minimization");
	$Rmol->setParameter(dockConsFname=>($prefix.".cons"))
	    if (&GenUtil::exists($prefix.".cons"));
	$sArr=$Rmol->dock($sArr,dockMin=>1);
	$Rmol->resetDockFlags();
	&GenUtil::log("finished docking minimization");

    } elsif ($cmd eq "dockfunnel") {

	my $outdir="dockfunnel_".$Rmol->reportParameter("series");
	$Rmol->setParameter(outdir=>$outdir);
	if (defined $args->{explode}) {
	  $Rmol->setParameter(nstruct=>$args->{explode});
	} else {
	  $Rmol->setParameter(nstruct=>500);
	}
	&GenUtil::log("starting docking funnel step");
	$sArr=$Rmol->dockGenerateFunnel($sArr);
	$Rmol->resetDockFlags();
	&GenUtil::log("finished docking funnel step");

    } elsif ($cmd =~ "des") {

	# Support for commands "des", "des_int", "desock", "desrub", and "desSwapModules"
	my $outdir=$cmd."_".$Rmol->reportParameter("series");
	$Rmol->setParameter(outdir=>$outdir);
	$Rmol->setParameter(nstruct=>1);

#	$Rmol->setParameter(sidechain_entropy_penalty=>0.2);

	$Rmol->setParameter(mcmin=>0);
	$Rmol->setParameter(soft_rep=>0);
	if ((defined $args->{soft_rep}) && ($args->{soft_rep} != 0)) {
	    # If using soft_rep, use mcmin
	    $Rmol->setParameter(soft_rep=>1);
	    $Rmol->setParameter(mcmin=>1);
	}
	$Rmol->setParameter(mcmin=>1)
	  if ((defined $args->{mcmin}) && ($args->{mcmin} != 0));

#	my $nativePDB=$prefix.".native.pdb.gz";
#	$Rmol->setParameter(favored_PDB_seq=>$nativePDB);
#	$Rmol->setParameter(favor_native=>-1.0);

	foreach my $c ( @{$Rmol->activeChains()} ) {
	    if ($c->{id} =~ /[M-Z]/) {
		$Rmol->setParameter(fix_target_seq=>1);
	    }
	}

        # Use a predefined resfile if desired
	if ($cmd eq "desSwapModules") {
	    $Rmol->setParameter(desSwapLoops=>1);
	    $Rmol->setParameter(moduleContactRank=>$args->{moduleContactRank})
		if (defined $args->{moduleContactRank});
	    $Rmol->setParameter(des_int=>1);
	} elsif ($cmd eq "desock") {
	    $Rmol->setParameter(desock=>1);
	} elsif ($cmd eq "desrub") {
	    $Rmol->setParameter(desrub=>1);
	} else {
	    $Rmol->setParameter(des_int=>1);
	}

	# If an appropriate ARinfo file exists, point to this
	my @ARfilelist=glob("????.ARinfo");
        if ($#ARfilelist == 0) {
	    my $ARinfoFile=shift(@ARfilelist);
	    $Rmol->setParameter(ARinfo=>$ARinfoFile);
	}

	$Rmol->setParameter(no_new_CG=>1);
	$Rmol->setParameter(interface_ds=>1);
	$Rmol->setParameter(nstruct=>$args->{explode}) if (defined $args->{explode});

	if ((defined $args->{incorporate_hotspot}) && ($args->{incorporate_hotspot} != 0)) {
	  my $nativePDB=$prefix.".native.pdb";
	  $Rmol->setParameter(hotspot_restore_sidechains=>$nativePDB);
	  $Rmol->setParameter(dockConsFname=>($prefix.".cons"))
	    if (&GenUtil::exists($prefix.".cons"));
	  $Rmol->setParameter(dockPert=>"2 2 2");
	}

        # Run design mode
	&GenUtil::log("starting ".$cmd);
	$sArr=$Rmol->design($sArr);
	&GenUtil::log("finished ".$cmd);

	if ((defined $args->{incorporate_hotspot}) && ($args->{incorporate_hotspot} != 0)) {
	  $Rmol->resetDockFlags();
	}

    } elsif ($cmd eq "ref") {

	my $outdir="ref_".$Rmol->reportParameter("series");
	$Rmol->setParameter(outdir=>$outdir);
	$Rmol->setParameter(nstruct=>$args->{explode}) if (defined $args->{explode});

        # Run refinement
	&GenUtil::log("starting refinement");
	$sArr=$Rmol->loopFaRlx($sArr);
	&GenUtil::log("finished refinement");

    } elsif ($cmd eq "chamin") {

	my $outdir="chamin_".$Rmol->reportParameter("series");
	$Rmol->setParameter(outdir=>$outdir);

        # Call charmm
	&GenUtil::log("starting chamin");
	$Rmol->setParameter(chamin_int=>1) if ($Rmol->numchains()>1);
	$sArr=$Rmol->charmm_mini($sArr);
	&GenUtil::log("finished chamin");

    } elsif ($cmd eq "restoreNATRO") {

	&GenUtil::log("restoring native sidechains at NATRO positions");

	my $nativePDB=&UWutil::PDBfname($prefix.".native.pdb");
	my $Nmol=Molecule->new($nativePDB);
	my @chains=@{$Nmol->activeChains()};
	my $first_chain=shift(@chains);
	my $second_chain=shift(@chains);
	my $fixed_cid=$first_chain->{id};
	my $moving_cid=$second_chain->{id};

	foreach my $f ( @{$sArr}) {

	    my $infile = $f->{filename};
	    my $prefix=&PDButil::prefixAndPath($infile);
	    my $outfile = $prefix.".restored.pdb";
	    my $compress=&GenUtil::validGzip($infile);
	    $outfile .= ".gz" if ($compress);

	    if (! &GenUtil::exists($outfile)) {
	        my $Bmol=Molecule->new($infile);

		# copy NATRO (per Rmol) coors of the native fixed chain to Bmol
		my $Rc=$Rmol->getChain($fixed_cid);
		my $Nc=$Nmol->getChain($fixed_cid);
		my $Bc=$Bmol->getChain($fixed_cid);

		if ( ( $#{ $Rc->{res} } == 0 ) || ( $#{ $Nc->{res} } == 0 ) || ( $#{ $Bc->{res} } == 0 ) ) {
		  die "Error with chain IDs";
		}

		for ( my $ir = 0 ; $ir <= $#{ $Rc->{res} } ; $ir++ ) {
		  if ( $Rc->{res}->[$ir]->{des} eq "NATRO" ) {
		    my $Nr = $Nc->{res}->[$ir];
		    my $Br = $Bc->{res}->[$ir];
		    my $Ni=$Nr->{start};
		    for (my $Bi=$Br->{start}; $Bi<=$Br->{end}; $Bi++) {
		      if ( $Bc->{atom}->[$Bi]->{atomname} ne $Nc->{atom}->[$Ni]->{atomname} ) {
			my $resn = $Bc->{atom}->[$Bi]->{resname};
			if ( $Bc->{atom}->[$Bi]->{hyd} && $Nc->{atom}->[$Ni]->{hyd} &&
			     ( ( $resn eq "HIS" ) || ( $resn eq "HSD" ) || ( $resn eq "HSE" ) ) ) {
			  # Do nothing with mismatched His hydrogens - let optimize H handle these
			} else {
			  printf STDERR "Chain ID is: %s\n", $moving_cid;
			  printf STDERR "Residue is: %s%d\n", $Bc->{atom}->[$Bi]->{resname},
			    $Bc->{atom}->[$Bi]->{resnum};
			  printf STDERR "Atomnames are is: %d\n", $Bc->{atom}->[$Bi]->{atomname},
			    $Nc->{atom}->[$Ni]->{atomname};
			  die "Error: atom name mismatch";
			}
		      } else {
			$Bc->{atom}->[$Bi]->{xcoor}=$Nc->{atom}->[$Ni]->{xcoor};
			$Bc->{atom}->[$Bi]->{ycoor}=$Nc->{atom}->[$Ni]->{ycoor};
			$Bc->{atom}->[$Bi]->{zcoor}=$Nc->{atom}->[$Ni]->{zcoor};
		      }
		      $Ni++;
		    }
		  }
		}

		# orient the "moving" chain of the native onto the docked structure
		$Bmol->selectChain($moving_cid);
		my $oneB=$Bmol->clone(1);
		$Nmol->selectChain($moving_cid);
		my $oneN=$Nmol->clone(1);
		$Nmol->selectChain(undef);
		my $analyze=Analyze::new($oneB);
		$analyze->lsqfit($oneN,"ca",0,1);

		# copy NATRO (per Rmol) coors of the (aligned) native moving chain to Bmol
		$Rc=$Rmol->getChain($moving_cid);
		$Bc=$Bmol->getChain($moving_cid);
		$Nc=$oneN->getChain($moving_cid);

		for ( my $ir = 0 ; $ir <= $#{ $Rc->{res} } ; $ir++ ) {
		  if ( $Rc->{res}->[$ir]->{des} eq "NATRO" ) {
		    my $Nr = $Nc->{res}->[$ir];
		    my $Br = $Bc->{res}->[$ir];
		    my $Ni=$Nr->{start};
		    for (my $Bi=$Br->{start}; $Bi<=$Br->{end}; $Bi++) {
		      if ( $Bc->{atom}->[$Bi]->{atomname} ne $Nc->{atom}->[$Ni]->{atomname} ) {
			my $resn = $Bc->{atom}->[$Bi]->{resname};
			if ( $Bc->{atom}->[$Bi]->{hyd} && $Nc->{atom}->[$Ni]->{hyd} &&
			     ( ( $resn eq "HIS" ) || ( $resn eq "HSD" ) || ( $resn eq "HSE" ) ) ) {
			  # Do nothing with mismatched His hydrogens - let optimize H handle these
			} else {
			  printf STDERR "Chain ID is: %s\n", $moving_cid;
			  printf STDERR "Residue is: %s%d\n", $Bc->{atom}->[$Bi]->{resname},
			    $Bc->{atom}->[$Bi]->{resnum};
			  printf STDERR "Atomnames are is: %d\n", $Bc->{atom}->[$Bi]->{atomname},
			    $Nc->{atom}->[$Ni]->{atomname};
			  die "Error: atom name mismatch";
			}
		      } else {
			$Bc->{atom}->[$Bi]->{xcoor}=$Nc->{atom}->[$Ni]->{xcoor};
			$Bc->{atom}->[$Bi]->{ycoor}=$Nc->{atom}->[$Ni]->{ycoor};
			$Bc->{atom}->[$Bi]->{zcoor}=$Nc->{atom}->[$Ni]->{zcoor};
		      }
		      $Ni++;
		    }
		  }
		}

		# write the output PDB
		$Bmol->selectChain(undef);
		$Bmol->_coorCache();
		$Bmol->writePDB($outfile, undef, "compress"=>$compress);
	    }

	    $f->{filename}=$outfile;
	}

	&GenUtil::log("native NATRO sidechains restored");

    } elsif ($cmd eq "restore") {

	&GenUtil::log("restoring native sidechains");

	# Note: docking only moves the second chain, this is
	# the only one we'll replace with the native
	my $nativePDB=&UWutil::PDBfname($prefix.".native.pdb");
	my $Nmol=Molecule->new($nativePDB);
	my @chains=@{$Nmol->activeChains()};
	my $first_chain=shift(@chains);
	my $second_chain=shift(@chains);
	my $fixed_cid=$first_chain->{id};
	my $moving_cid=$second_chain->{id};

	foreach my $f ( @{$sArr}) {

	    my $infile = $f->{filename};
	    my $prefix=&PDButil::prefixAndPath($infile);
	    my $outfile = $prefix.".restored.pdb";
	    my $compress=&GenUtil::validGzip($infile);

	    if (! &GenUtil::exists($outfile)) {
	        my $Bmol=Molecule->new($infile);
		$Bmol->selectChain($moving_cid);
		my $oneB=$Bmol->clone(1);
		$Nmol->selectChain($moving_cid);
		my $oneN=$Nmol->clone(1);
		$Nmol->selectChain(undef);
		# Replace the existing PDB
		my $analyze=Analyze::new($oneB);
		$analyze->lsqfit($oneN,"ca",0,1);
		$Nmol->merge($oneN);
		$Nmol->writePDB($outfile, undef, "compress"=>$compress);
		$outfile .= ".gz" if ($compress);
		&GenUtil::remove($infile) if ( &GenUtil::exists($infile) );
	    }

	    $f->{filename}=$outfile;
	}

	&GenUtil::log("native sidechains restored");

    } elsif ($cmd eq "testhelices") {

	my $sel=uc($args->{sel});
	$sel="CA" if (! defined $sel);
	die "Unknown sel for helix filter"
	    if (($sel ne "CA") && ($sel ne "CB"));

	my $cutoff=$args->{cutoff};
	if (! defined $cutoff) {
	    if ($sel eq "CA") {
		$cutoff=4.7;
	    } elsif ($sel eq "CB") {
		$cutoff=3.8;
	    }
	}
	my $SQcut=$cutoff*$cutoff;

	my $helixFile=$prefix.".helices";

	if (&GenUtil::zexists($helixFile)) {

            # Read the .func file
	    &GenUtil::log("testing helices");
	    my @helices=();
	    my $inh=&GenUtil::getInputFile($helixFile);
	    while (<$inh>) {
		chomp;
		my @p=split(' ',$_);
		my $h={};
		$h->{cid}=$p[0];
		$h->{start}=$p[1]+0;
		$h->{stop}=$p[2]+0;
		push(@helices,$h);
	    }
	    undef $inh;

            # Keep only functional conformations
	    my $tmpArr=();
            while ($#{$sArr}>=0) {
	      my $f=shift(@{$sArr});

	      # Collect CA coors
	      my $ca_coors={};
	      my $inf=&GenUtil::getInputFile($f->{filename});
	      while (<$inf>) {
		if (/^ATOM/) {
		  if (/$sel/) {
		    chomp;
		    my $cid=substr($_,21,1);
		    (my $resnum=substr($_,22,5))=~s/ //g;
		    my $X=substr($_,30,8)+0.0;
		    my $Y=substr($_,38,8)+0.0;
		    my $Z=substr($_,46,8)+0.0;
		    $ca_coors->{$cid.$resnum."X"}=$X;
		    $ca_coors->{$cid.$resnum."Y"}=$Y;
		    $ca_coors->{$cid.$resnum."Z"}=$Z;
		  }
	        }
	      }

	      # Examine minimum helix-helix CA distances
	      my $minSQdist=999.9;
	      for (my $A=0; $A < $#helices; $A++) {
		my $hA=$helices[$A];
		for (my $B=$A+1; $B <= $#helices; $B++) {
		  my $hB=$helices[$B];
		  if ($hB->{cid} ne $hA->{cid}) {

		    for (my $iA=$hA->{start}; $iA <= $hA->{stop}; $iA++) {
		      if (defined $ca_coors->{$hA->{cid}.$iA."X"}) {
			for (my $iB=$hB->{start}; $iB <= $hB->{stop}; $iB++) {
			  if (defined $ca_coors->{$hB->{cid}.$iB."X"}) {
			    my $t=$ca_coors->{$hA->{cid}.$iA."X"}-
				$ca_coors->{$hB->{cid}.$iB."X"};
			    my $currdist=$t*$t;
			    $t=$ca_coors->{$hA->{cid}.$iA."Y"}-
				$ca_coors->{$hB->{cid}.$iB."Y"};
			    $currdist+=$t*$t;
			    $t=$ca_coors->{$hA->{cid}.$iA."Z"}-
				$ca_coors->{$hB->{cid}.$iB."Z"};
			    $currdist+=$t*$t;
			    if ($currdist < $minSQdist) {
				$minSQdist=$currdist;
			    }
			  }
		        }
		      }
		    }
		  }
		}
	      }

	      if ($minSQdist >= $SQcut) {
		  push(@{$tmpArr},$f);
	      } else {
#		    &GenUtil::remove($f->{filename});
#		    &GenUtil::remove($f->{Rosout});
	      }
	    }
	    $sArr=$tmpArr;
	}

    } elsif ($cmd eq "eval") {

	my $efield=$args->{key};
	# Test whether to evaluate the value of the sort field
	if (! defined $sArr->[0]->{$efield}) {

  	  # Get it from the latest scorefile if possible,
          # otherwise look to the PDB files (or generate them)
	  my $gotFromFasc=0;
	  my $outdir=$Rmol->reportParameter("outdir");
	  my $fascfile=$outdir."scorefile.fasc";
	  $fascfile=$outdir."scorefile.sc"
	      if (! &GenUtil::zexists($fascfile));

	  if (&GenUtil::zexists($fascfile)) {
	      my $infile=&GenUtil::getInputFile($fascfile);
	      my $inx=-1;
	      my $firstline=<$infile>;
	      my @p=split(' ',$firstline);
	      my $lastinx=$#p;
	      for (my $i=0; $i<$#p; $i++) {
		  if ($efield eq $p[$i]) {
		      $inx=$i;
		      $gotFromFasc=1;
		  }
	      }
	      if ($inx >= 0) {

                  # Create a hash relating filenames to elements of sArr
                  # Note: match to the Rosout field in sArr
		  my $fhash={};
		  foreach my $f (@{$sArr}) {
		      $fhash->{$f->{Rosout}}=$f;
		  }

		  my $foundMatch=0;
		  while (<$infile>) {
		      if (/output_decoy$/) {
			  my $inline=$_;
			  chomp($inline);
			  my @p=split(' ',$inline);
			  # Require lastinx matches, because some lines
			  # can accidentally contain data for two structures
			  if ($#p == $lastinx) {
			      my $fname=$outdir.&PDButil::prefix($p[0]).".pdb";
			      my $f=$fhash->{$fname};
			      $f->{$efield}=($p[$inx]+0.0) if (defined $f);
			      $foundMatch=1;
			  }
		      }
		  }
		  undef $infile;

		  # If we find no matches at all, the headers on the scorefile
		  # must not match the data.
		  # In this case, just reread it without this check.
		  if ( ! $foundMatch ) {
		    my $infile=&GenUtil::getInputFile($fascfile);
		    while (<$infile>) {
		      if (/output_decoy$/) {
			my $inline=$_;
			chomp($inline);
			my @p=split(' ',$inline);
			my $fname=$outdir.&PDButil::prefix($p[0]).".pdb";
			my $f=$fhash->{$fname};
			$f->{$efield}=($p[$inx]+0.0) if (defined $f);
			$foundMatch=1;
		      }
		    }
		    undef $infile;
		  }
	     }
	}

	if (! $gotFromFasc) {

	      my $outMol=RosettaMolecule->new($sArr->[0]->{Rosout});
	      $outMol->updateScore($sArr->[0]->{Rosout});
	      if (defined $outMol->{score}->{$efield}) {
		  foreach my $f (@{$sArr}) {
		      $outMol->updateScore($f->{Rosout});
		      my $shash=$outMol->getScoreHash();
		      $f->{$efield}=$shash->{$efield};
		  }

	      } elsif ($efield eq "sscont") {
		  foreach my $f (@{$sArr}) {
		      my $outMol=Molecule->new($f->{filename});
		      my $contacts=$outMol->countInterfaceContacts(undef,
						       undef,undef,0,1);
		      my $SScontacts=0;
  		      if ($#{$contacts} >= 0) {
			  foreach my $nc ( @{$contacts} ) {
			      $SScontacts+=$nc->{AAcontacts};
			  }
		      }
		      $f->{$efield}=$SScontacts;
		  }

	      } elsif ($efield eq "nonUBQcont") {
		  foreach my $f (@{$sArr}) {
		      my $outMol=Molecule->new($f->{filename});
		      $f->{$efield}=$outMol->count_nonUBQcont_Contacts();
		  }

	      } elsif ($efield eq "sc") {

		  # JK TEMP
#		  my $infile=&GenUtil::getInputFile("list_sc.txt");
		  foreach my $f (@{$sArr}) {
#		      my $inline;
#		      if ($inline = <$infile>) {
#			  chomp($inline);
#			  $f->{$efield}=$inline;
#			  &GenUtil::log("reusing Sc value ".$f->{$efield});
#		      } else {
			  my $outMol=Molecule->new($f->{filename});
			  $f->{$efield}=$outMol->computeSc();
			  &GenUtil::log("found Sc value ".$f->{$efield});
#		      }
		  }
#		  undef $infile;

	      } elsif ($efield eq "terminus_sqdist") {
		  foreach my $f (@{$sArr}) {
		      my $outMol=Molecule->new($f->{filename});
		      $f->{$efield}=$outMol->measure_UBQ_terminus_sqdist();
		  }

	      } elsif ($efield eq "UNS_AspAsn") {
		  foreach my $f (@{$sArr}) {
		    my $outMol=KnownRepeatProteins::newProtein($f->{filename});
		    my $unsval=0;
		    if ($outMol->moduleName() eq "AnkyrinRepeat") {
		      $unsval=$outMol->find_UNS_AspAsn();
		    }
		    $f->{$efield}=$unsval;
		  }

	      } else {
		  die "Unsupported eval field: ".$efield;
	      }

	  }
        }

    } elsif ($cmd eq "sort") {

	my $sfield=$args->{key};
	# Evaluate the value of the sort field if it's not present
	if (! defined $sArr->[0]->{$sfield}) {
	    $sArr=&processInstruction($Rmol, "eval:key=".$args->{key}, $sArr);
	}

	# Sort the following in descending order (ie. big values are better)
	if (($sfield eq "ncont") || ($sfield eq "sscont") ||
	    ($sfield eq "nonUBQcont") ||
	    ($sfield eq "sc") ||
	    ($sfield eq "interface_aromatic_sasa") ||
	    (lc($sfield) =~ "sasaprob") ||
	    (($sfield =~ /num_/) && ($sfield ne "num_interface_waters"))) {
	  $sArr=[sort( {$b->{$sfield} <=> $a->{$sfield}} @{$sArr})];
	} else {
	  # Sort anything else in ascending order
	  $sArr=[sort( {$a->{$sfield} <=> $b->{$sfield}} @{$sArr})];
	}

    } elsif ($cmd eq "cutval") {

        # Keep only elements in which "key" is better or equal to "val"
        my $field=$args->{key};
        my $val=$args->{val};

        # Sort first - the first element is better than the last,
        # ie. the sort determines whether a "good" value is big or small
	$sArr=&processInstruction($Rmol, "sort:key=".$args->{key}, $sArr);
	my $keepsmall=1;

	# Sort the following in descending order (ie. big values are better)
	my $sfield=$args->{key};
	if (($sfield eq "ncont") || ($sfield eq "sscont") ||
	    ($sfield eq "nonUBQcont") ||
	    ($sfield eq "sc") ||
	    ($sfield eq "interface_aromatic_sasa") ||
	    (lc($sfield) =~ "sasaprob") ||
	    (($sfield =~ /num_/) && ($sfield ne "num_interface_waters"))) {
	    $keepsmall=0;
	}

	my $didpop=1;

	# If all elements are identical, we've presumably filtered
	# on this binary property earlier - so just keep everything
#	if ($sArr->[0]->{$field} == $sArr->[$#{$sArr}]->{$field}) {
#	  $didpop=0;
#	}

	while ($didpop) {
	    my $last=pop(@{$sArr});
	    if ($keepsmall) {
		if ($last->{$field} <= $val) {
		    push(@{$sArr},$last);
		    $didpop=0;
		}
	    } else {
		if ($last->{$field} >= $val) {
		    push(@{$sArr},$last);
		    $didpop=0;
		}
	    }
	    $didpop=0 if ($#{$sArr}<0);
	}

    } elsif ($cmd eq "cutbest") {

        # Keep the best "keep" elements based on "key"
        my $field=$args->{key};
        my $keep=$args->{keep};
	$sArr=&processInstruction($Rmol, "sort:key=".$field, $sArr);
        $#{$sArr}=int($keep)-1 if ($#{$sArr} > int($keep)-1);

    } elsif ($cmd eq "cutbest_perstart") {

        # Keep the best "keep" elements based on "key"
        my $field=$args->{key};
        my $keep=$args->{keep};

	$keep=1 if (! defined $keep);

	my $startHash={};
	# Separate all the output structures by starting structure
	# Note: the format of each element is preserved
	while ($#{$sArr} >= 0) {
	    my $outdec=shift(@{$sArr});
	    my $startStruct=$outdec->{origfile};
	    $startHash->{$startStruct}=()
		if (! defined $startHash->{$startStruct});
	    push(@{$startHash->{$startStruct}},$outdec);
	}

	# For each set of starting structures, apply a "cutbest" step,
	# using the parameters passed in
	foreach my $startStruct ( keys %{$startHash} ) {
	    my $mArr=&processInstruction($Rmol,
	      "cutbest:key=".$field.":keep=".$keep,$startHash->{$startStruct});
	    while ($#{$mArr} >= 0) {
		push(@{$sArr},shift(@{$mArr}));
	    }
	}

	# Finish with a final "sort" of the complete set
	$sArr=&processInstruction($Rmol, "sort:key=".$field, $sArr);

    } elsif ($cmd eq "allscores") {
        # Write the scores file
	my $scoreFile="allscores.".$Rmol->reportParameter("series").".out";
	&write_scores($scoreFile,$sArr);

    } else {
	die "Unsupported instruction: ".$cmd;
    }

    return $sArr;
}


############################################################################



