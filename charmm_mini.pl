#!/usr/bin/env perl

# perform CHARMM minimization on a series of PDB files
# Note: PDB files may differ in sequence, etc.
#

sub usage {
  printf STDERR "usage:   charmm_mini.pl [-sidechain_only] [-interface] inPDB|PDBlist\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use Sys::Hostname;

use GenUtil;
use PDButil;
use RosettaMolecule;
use CHARMM;

my $inarg;
my $sc_only=0;
my $interface=0;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif (lc($ARGV[0]) eq "-sidechain_only") {
    shift @ARGV;
    $sc_only=1;
  } elsif (lc($ARGV[0]) eq "-interface") {
    shift @ARGV;
    $interface=1;
  } else {
    $inarg = shift @ARGV;
  }
}

&usage() if (! defined $inarg);

my %par;

&GenUtil::parsePar(\%par,"param=22x");
&GenUtil::parsePar(\%par,"cutnb=20.0");
&GenUtil::parsePar(\%par,"cutoff=18.0");
&GenUtil::parsePar(\%par,"cuton=16.0");
&GenUtil::parsePar(\%par,"shake=1");
&GenUtil::parsePar(\%par,"gb");
&GenUtil::parsePar(\%par,"gbmvsa=0.015");
&GenUtil::parsePar(\%par,"gbmvas=0.9");
&GenUtil::parsePar(\%par,"gbmvap7=0.0001");
&GenUtil::parsePar(\%par,"gbmvad=-0.1");
&GenUtil::parsePar(\%par,"gbmvaig=32");
&GenUtil::parsePar(\%par,"scalerad=mnina");
&GenUtil::parsePar(\%par,"cmap");
&GenUtil::parsePar(\%par,"xpar=par_cmap27.inp");
&GenUtil::parsePar(\%par,"xtop=top_all22_prot_cmap.inp");
&GenUtil::parsePar(\%par,"gbmvafrq=1");
&GenUtil::parsePar(\%par,"gbmvaemp=999999.0");

&GenUtil::parsePar(\%par,"sdsteps=250");
&GenUtil::parsePar(\%par,"minsteps=500");

my $charmm=&CHARMM::new();
$charmm->loadParameters(%par);

my $inlist=();

if (&PDButil::isaPDB($inarg)) {
    my $f={};
    $f->{in}=$inarg;
    my $outfile=lc(&PDButil::prefix($inarg).".mini.pdb");
    $f->{out}=$outfile;
    push(@{$inlist},$f);

} else {
    my $inp=GenUtil::getInputFile($inarg);
    while (<$inp>) {
	my $inline=$_;
	chomp($inline);
	my @p=split(' ',$inline);
	my $infile=$p[0];
	my $outfile;
	if ($#p == 1) {
	    $outfile=$p[1];
	} else {
	    $outfile=lc(&PDButil::prefix($infile).".mini.pdb");
	}
	my $f={};
	$f->{in}=$infile;
	$f->{out}=$outfile;
	push(@{$inlist},$f);
    }
}

foreach my $f ( @{$inlist} ) {

    my $fname=$f->{in};
    my $outfile=$f->{out};

    if ((&GenUtil::zexists($fname)) && (! &GenUtil::exists($outfile))) {

	system("touch ".$outfile);

	$charmm->setupFromPDB($fname,undef,undef,1);

	$charmm->clearEnergy();
	$charmm->setupEnergy();
	$charmm->verbose("cons fix sele type n .or. type ca .or. type c .or. type o end purg")
	    if ($sc_only);

# Example of "fix all but residues 71 and 78":
#	$charmm->verbose("cons fix sele (.not. resi 71) .and. (.not. resi 78) end purg");


	if ($interface) {
	    $charmm->verbose("define intA sele ((iseg 1) .and. ((.not. iseg 1) .around. 6.0) ) end");
	    $charmm->verbose("define intB sele ((iseg 2) .and. ((.not. iseg 2) .around. 6.0) ) end");
	    $charmm->verbose("define intC sele ((iseg 3) .and. ((.not. iseg 3) .around. 6.0) ) end");
	    $charmm->verbose("define intD sele ((iseg 4) .and. ((.not. iseg 4) .around. 6.0) ) end");
	    $charmm->verbose("cons fix sele .not. (intA .or. intB .or. intC .or. intD) show end purg");
	}

	$charmm->minimizeSD();
	$charmm->minimize();

	my $tmpfile=lc(hostname."-pdb$$");
	$charmm->writePDB($tmpfile);
	&GenUtil::rename($tmpfile,$outfile);
	$charmm->verbose("delete atom sele all end");

	my $f={};
	$f->{Rosout}=$outfile;
	$f->{filename}=$outfile;
	&RosettaMolecule::_fitToTemplate($fname,$f);

    }

}

$charmm->finish();
exit 0;

