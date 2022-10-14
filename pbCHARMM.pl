#!/usr/bin/env perl

# get Poisson-Boltzmann energy from PDB file 
#
# http://mmtsb.scripps.edu/doc/pbCHARMM.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   pbCHARMM.pl [options] [PDBfile]\n";
  printf STDERR "options: [-par param=19|22,hsd=list,hse=list,scalerad,\n";
  printf STDERR "               smooth,dcel=value,epsp=value,epsw=value,epsr=value]\n";
  printf STDERR "         [-psf PSFfile CRDfile]\n";
  printf STDERR "         [-partial file] [-threads n]\n";
  printf STDERR "         [-atomic] [-pairs] [-keepcharge]\n";
  printf STDERR "         [-log logFile] [-cmd logFile]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use IO::File;
use IO::Handle;

use GenUtil;
use Molecule;
use CHARMM;

my %par;

my $logFile;
my $cmdlog;

my $inpfile="-";
my $base="";

my $single;
my $double;
my $keepcharge;
my $threads=1;

my $partial;

my $psffile;
my $crdfile;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logFile=(shift @ARGV);
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdlog=(shift @ARGV);
  } elsif ($ARGV[0] eq "-atomic") {
    shift @ARGV;
    $single=1;
  } elsif ($ARGV[0] eq "-pairs") {
    shift @ARGV;
    $double=1;
  } elsif ($ARGV[0] eq "-threads") {
    shift @ARGV;
    $threads=shift @ARGV;
  } elsif ($ARGV[0] eq "-partial") {
    shift @ARGV;
    $partial=shift @ARGV;
  } elsif ($ARGV[0] eq "-keepcharge") {
    shift @ARGV;
    $keepcharge=1;
  } elsif ($ARGV[0] eq "-psf") {
    shift @ARGV;
    $psffile=shift @ARGV;
    $crdfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } else {
    die "Unknown option $ARGV[0]" if ($ARGV[0]=~/^-/);
    $inpfile=(shift @ARGV);
    $done=1;
  }
}

if (defined $single || defined $double) {
  my $tag=$$;
  my $mol=Molecule->new();

  if (defined $psffile) {
    $mol->readCRD($crdfile);
  } else {
    $mol->readPDB($inpfile);
    $mol->fixHistidine($par{hsd},$par{hse});
    $mol->translate(&CHARMM::getConvType($par{param}));
    $mol->generateSegNames();
  }

  my %pb;

  if (defined $partial && -r $partial) {
    open INP,"$partial";
    while (<INP>) {
      chomp;
      my ($inx,$pbe)=split(/ +/);
      $pb{$inx}=$pbe;
    }
    close INP;
  }

  my $natom=$#{$mol->{atom}}+1;

  my @todolist=();
  my @allist=();
  my $inx=0;

  if (defined $double) {
    my @atomlist=();
    foreach my $c ( @{$mol->{chain}} ) {
      foreach my $a ( @{$c->{atom}} ) { 
       my $jnx=$c->{id}.":".$a->{resnum}.":".$a->{atomname}.":".$inx;
       $inx++;
       push(@atomlist,$jnx);
      }
    }

    for (my $i=0; $i<$#atomlist; $i++) {
     for (my $j=$i; $j<=$#atomlist; $j++) {
      my $jnx=$atomlist[$i]."=".$atomlist[$j];  
      if (!exists $pb{$jnx}) {
	push(@todolist,$jnx);
      }
      push(@allist,$jnx);
     }
    }
  } elsif (defined $single) { 
   foreach my $c ( @{$mol->{chain}} ) {
    foreach my $a ( @{$c->{atom}} ) { 
      my $j=$c->{id}.":".$a->{resnum}.":".$a->{atomname}.":".$inx;
      $inx++;
      if (!exists $pb{$j}) {
	push(@todolist,$j);
      }
      push(@allist,$j);
    }
   }
  }
  
  my @pidlist;
  for (my $i=1; $i<=$threads; $i++) {
    my $pid=fork();
    if (!$pid) {
      $logFile="$i-$logFile" if (defined $logFile);
      my $outfile="$tag-$i-pbout";
      my $handle=&GenUtil::getOutputFile($outfile);
      $handle->autoflush(1);
      my $charmm=&CHARMM::new($logFile,$cmdlog);

      $charmm->loadParameters(%par);

      if (defined $psffile) {
	$charmm->setupFromPSF($psffile,$crdfile);
      } else {
	$charmm->setupFromMolecule($mol);
      }

      $charmm->orient();
      for (my $ia=$i-1; $ia<=$#todolist; $ia+=$threads) {
	my $inx=$todolist[$ia];
	my $pbener=$charmm->atomPoissonBoltzmann($inx,$keepcharge);
	printf $handle "%s %f\n",$inx,$pbener;
	printf STDERR "%s %f\n",$inx,$pbener;
      }
      $charmm->finish();
      close $handle;
      exit 0;
    } else {
      push (@pidlist,$pid);
    }
  } 
  foreach my $p (@pidlist) {
    waitpid($p,0);
  }

  for (my $i=1; $i<=$threads; $i++) {
    open INP,"$tag-$i-pbout";
    while (<INP>) {
      chomp;
      my ($inx,$pbe)=split(/ +/);
      $pb{$inx}=$pbe;
    }
    close INP;
    unlink "$tag-$i-pbout";
  }

  foreach my $pkey (@allist) { 
    printf "%s %f %f\n",$pkey,$pb{$pkey},-163.9603525/$pb{$pkey}
      if ($pb{$pkey} != 0.0);
  }
} else {
  my $charmm=&CHARMM::new($logFile,$cmdlog);
  $charmm->loadParameters(%par);

  if (defined $psffile) {
    $charmm->setupFromPSF($psffile,$crdfile);
  } else {
    $charmm->setupFromPDB($inpfile);
  }

  $charmm->orient();
  my $pbener=$charmm->poissonBoltzmann();
  $charmm->finish();

  printf "%f\n",$pbener;
}

exit 0;

