#!/usr/bin/env perl

# compute Sc (shape complementarity) for a PDB file
#
# http://mmtsb.scripps.edu/doc/computeSc.pl.html
# 2008, John Karanicolas, Baker lab, UW

sub usage {
  printf STDERR "usage:   computeSc.pl [PDBfile]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;
use Molecule;

my $fname="-";

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } else {
    $fname = shift @ARGV;
  }
}

my $mol=Molecule->new();
$mol->readPDB($fname);
printf STDOUT "Sc value is: %f\n", $mol->computeSc();

exit(0);

