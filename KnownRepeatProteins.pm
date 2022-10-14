# A wrapper for the specific repeat proteins available
# (ie. derived classes of RepeatProtein.pm)
# 
# http://mmtsb.scripps.edu/doc/KnownRepeatProteins.pm.html
# 2005, John Karanicolas, Baker group, UW

package KnownRepeatProteins;

require 5.004;

use strict;
use UWutil;

use AnkyrinRepeat;
use TPR;


## function: $outMol = newProtein(inPDB)
## Acts as a general wrapper for the repeat protein derived classes.
## Tries to setup input PDB using each derived class of RepeatProtein.
## If one looks good (ie. modules detected), an object of this type
## is returned. Otherwise, an object of generic type RepeatProtein is
## returned. This is essentially a wrapper for a constructor of unknown
## type.

sub newProtein {
  my $inPDB=shift;

  $inPDB=&UWutil::PDBfname($inPDB);

  die "inPDB is not defined" if (! defined $inPDB);

  my $outMol=AnkyrinRepeat->new($inPDB);

  if ($outMol->numModules <= 0) {
      $outMol=TPR->new($inPDB);
  }

  if ($outMol->numModules <= 0) {
      $outMol=RepeatProtein->new($inPDB);
  }

#  &GenUtil::log("setting up class: ".$outMol->moduleName());
  return $outMol;

}

1;
