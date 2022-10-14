# TPR package
# a derived class of RepeatProtein, specific for TPRs
#
# 2005, John Karanicolas, Baker group, UW
#
# derived from RepeatProtein

package TPR;

require 5.004;

use strict;

use FileHandle;
use IPC::Open2;

use RepeatProtein;
use GenUtil;
use Sequence;

use vars qw ( @ISA );

@ISA = ("RepeatProtein");

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

## data: repeat -> { chain, positionOffset, module }
## specification of the repeat

## data: module[] -> { start, end }
## specification of the repeat

use vars qw ( $moduleName $repeatLength );

BEGIN {
  $moduleName="TPR";
  $repeatLength=34;
}


## constructor: new([PDBfile],[parameters])
## creates a new TPR object and reads a PDB structures
## if a file name is given

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $farg=shift;
  my @pArr=@_;

  my $self = $class->SUPER::new($farg,@pArr);

  $self->detectRepeats() if (defined $farg);

  return $self;
}

# Destructor
DESTROY {
  my $self=shift;
  $self->SUPER::DESTROY();
  return;
}


## method: readPDB(file[,translate[,ignoreseg]])
## read an input PDB

sub readPDB {
  my $self=shift;
  my @args=@_;

  $self->SUPER::readPDB(@args);

  $self->detectRepeats();

  return;
}

## method: removeChain()

sub removeChain {
  my $self = shift;
  my $cid = shift;

  $self->resetRepeat() if ($self->{repeat}->{cid} eq $cid);
  $self->SUPER::removeChain($cid);
  $self->detectRepeats();

  return;
}

## method: evaluateModuleScore(resArr,startInx)
## search the sequence of the current chain for an TPR repeat

sub evaluateModuleScore {
  my $self=shift;
  my $resArr=shift;
  my $startInx=shift;

  # We'll define TPR modules starting at.... JK FILL THIS IN
#  printf STDERR "JK Warning: evaluateModuleScore for TPR repeats not yet written\n";
  return 0;

  my $module_score=0.;

  # We'll define TPR modules starting at the conserved Gly between
  # the helices (position 13 of the Mosavi et al. 2002 PNAS paper),
  # since this is where we wish to "swap out" modules

  # 3 points for each "in-frame" well-conserved residue
  foreach my $pos ( @{$self->listWellConserved()} ) {
      my $currRes=$resArr->[$startInx+$pos->{offset}]->{name};
      $module_score+=3
	  if ($pos->{id} =~ $Sequence::_seqabbrev{$currRes});
  }

  # 1 point for each "in-frame" semi-conserved residue
  foreach my $pos ( @{$self->listSemiConserved()} ) {
      my $currRes=$resArr->[$startInx+$pos->{offset}]->{name};
      $module_score+=1
	  if ($pos->{id} =~ $Sequence::_seqabbrev{$currRes});
  }

  $module_score /= 15.;

  return $module_score;
}


## method: listNonBinding()
## return a list residues which are not allowed to participate in the interface

sub listNonBinding {
  my $self=shift;

  my $nonBinding=();

  # We'll define TPR modules starting at.... JK FILL THIS IN
  die "listNonBinding for TPR repeats not yet written";

  return $nonBinding;
}


## method: listWellConserved()
## return a list well-conserved residues

sub listWellConserved {
  my $self=shift;

  my $wellConserved=();

  # We'll define TPR modules starting at.... JK FILL THIS IN
  die "listWellConserved for TPR repeats not yet written";

  return $wellConserved;
}


## method: listSemiConserved()
## return a list semi-conserved residues

sub listSemiConserved {
  my $self=shift;

  my $semiConserved=();

  # We'll define TPR modules starting at.... JK FILL THIS IN
  die "listSemiConserved for TPR repeats not yet written";

  return $semiConserved;
}


## method: moduleName()
## return the name of this type of module (ie. TPR)

sub moduleName {
  my $self=shift;

  return $moduleName;
}

## method: repeatLength()
## return list of residue numbers of highly conserved positions

sub repeatLength {
  my $self=shift;

  return $repeatLength;
}

## method: repeating()
## returns 1 to indicate that this object is of a class with repeats

sub repeating {
  my $self=shift;

  return 1;
}


1;

