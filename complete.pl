#!/usr/bin/env perl
#
# adds missing atoms to a PDB structure 
#
# http://mmtsb.scripps.edu/doc/complete.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   complete.pl [options] [PDBfile]\n";
  printf STDERR "options: [-hsd list] [-hse list] [-param 19|22] [-blocked]\n";
  printf STDERR "         [-log file] [-cmd file]\n";
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
use Molecule;
use CHARMM;

my %par = (
  hsd => "",
  hse => "",
  param => 22,
  blocked => 0
);

my $fname;

my $logfile;
my $cmdfile;

my $fraglist;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-hsd") {
    shift @ARGV;
    $par{hsd}=shift @ARGV;
  } elsif ($ARGV[0] eq "-hse") {
    shift @ARGV;
    $par{hse}=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-param") {
    shift @ARGV;
    $par{param}=shift @ARGV;
  } elsif ($ARGV[0] eq "-blocked") {
    shift @ARGV;
    $par{blocked}=1;
  } else {
    $fname=shift @ARGV;
  }
}

my $mol=Molecule->new();
$mol->readPDB($fname);

my $ss=$mol->getSSBonds();

$mol->completeResidue();
$mol->fixCOO() if (!$par{blocked});

$mol->translate("CHARMM22");
$mol->fixHistidine($par{hsd},$par{hse});
$mol->generateSegNames();

my $charmm=&CHARMM::new($logfile,$cmdfile);
$charmm->loadParameters(%par);
$charmm->setupFromMolecule($mol);

my $chmoutpdb="t$$.out.pdb";
$charmm->writePDB($chmoutpdb);
my $outmol=Molecule->new();
$outmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($par{param}),chainfromseg=>1);

$outmol->setSSBonds($ss);
$outmol->writePDB("-",translate=>"CHARMM22");
&GenUtil::remove($chmoutpdb);
  
$charmm->finish();

