#!/usr/bin/env perl

# count the number of heavy-atom contacts at the interface between chains
#
# http://mmtsb.scripps.edu/doc/countIntCont.html
# 2003, John Karanicolas, Baker group, UW

sub usage {
  printf STDERR "usage:   countIntCont.pl [options] [PDBfile]\n";
  printf STDERR "options: [-SConly] [-distcut num] [-chainA id] [-chainB id]\n";
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
my $distcut=5;
my $chainA;
my $chainB;
my $SConly=0;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif (lc($ARGV[0]) eq "-distcut") {
    shift @ARGV;
    $distcut=shift @ARGV;
  } elsif (lc($ARGV[0]) eq "-chaina") {
    shift @ARGV;
    $chainA=shift @ARGV;
  } elsif (lc($ARGV[0]) eq "-chainb") {
    shift @ARGV;
    $chainB=shift @ARGV;
  } elsif (lc($ARGV[0]) eq "-sconly") {
    shift @ARGV;
    $SConly=1;
  } else {
    $fname = shift @ARGV;
  }
}

my $mol=Molecule->new($fname);
my $contacts=$mol->countInterfaceContacts($distcut,$chainA,$chainB,undef,$SConly);

if ($#{$contacts} >= 0) {
    foreach my $nc ( @{$contacts} ) {
	printf STDOUT "%d atomic contacts between chains %s and %s\n",
   	  $nc->{AAcontacts}, $nc->{cidA}, $nc->{cidB};
    }
}

