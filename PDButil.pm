# Utilities for PDB manipulation
# 
# 2005, John Karanicolas, Baker group, UW

package PDButil;

require 5.004;

use strict;

use FileHandle;
use IPC::Open2;
use Sys::Hostname;

use GenUtil;
use Molecule;


## function: string = prefix(filename)
## chop off either .pdb or .pdb.gz, remove dir name, return the remainder

sub prefix {
  my $f=shift;

  # Remove suffix
  $f=substr($f,0,-2) if ($f =~ /gz$/);
  $f=substr($f,0,-1) if ($f =~ /\.$/);
  $f=substr($f,0,-3) if ($f =~ /pdb$/);
  $f=substr($f,0,-1) if ($f =~ /\.$/);

  # Remove path
  $f=substr($f,1) if ($f =~ /^\//);
  $f =~ s/.*\///;

  return $f;
}

## function: string = prefixAndPath(filename)
## chop off either .pdb or .pdb.gz, return the remainder (preserves dir name)

sub prefixAndPath {
  my $f=shift;

  # Remove suffix
  $f=substr($f,0,-2) if ($f =~ /gz$/);
  $f=substr($f,0,-1) if ($f =~ /\.$/);
  $f=substr($f,0,-3) if ($f =~ /pdb$/);
  $f=substr($f,0,-1) if ($f =~ /\.$/);

  return $f;
}

## function: bool = isaPDB(filename)
## return whether or not this filename is a PDB file

sub isaPDB {
  my $f=shift;

  return 1 if ($f =~ /\.pdb$/);
  return 1 if ($f =~ /\.pdb\.gz$/);
  return 0;
}

## function: print_header(PDBfname)
## print header lines of the specified PDB file to STDOUT

sub print_header {
  my $PDBfname=shift;
  my $outFname=shift;

  my $infile=&GenUtil::getInputFile($PDBfname);
  my $outfile=&GenUtil::getOutputFile($outFname);

  while (<$infile>) {
      if ((! /^(ATOM|HETATM|MODEL|CONECT|SCALE|ORIG|CRYST|MASTER|TER|END)/)
	  && (! /REMARK 290/)) {
	  my $inline=$_;
	  chomp($inline);
	  printf $outfile "%s\n",$inline;
      }
  }

  undef $infile;
  undef $outfile;
  return;
}


## function: shortenRasmolList(string)
## Given list of comma separated residues, use dashes where possible

sub shortenRasmolList {
  my $str=shift;

  my @rlist=split(/,/,$str);
  return if ($#rlist < 0);

  my @slist=sort { $a <=> $b } @rlist;

  my $prevMin;
  my $prevMax;
  my $currRes;
  my @rstr=();
  while ($#slist>=0) {
      my $currRes=shift(@slist);
      if (! defined $prevMin) {
	  $prevMin=$currRes;
	  $prevMax=$currRes;
      } elsif ($currRes == ($prevMax+1)) {
	  $prevMax=$currRes;
      } else {
	  my $currStr="";
	  if ($prevMin == $prevMax) {
	      $currStr=$prevMin;
	  } else {
	      $currStr=$prevMin."-".$prevMax;
	  }
	  push(@rstr,$currStr);
	  $prevMin=$currRes;
	  $prevMax=$currRes;
      }
  }

  my $currStr="";
  if ($prevMin == $prevMax) {
      $currStr=$prevMin;
  } else {
      $currStr=$prevMin."-".$prevMax;
  }
  push(@rstr,$currStr);

  return join(",",@rstr);
}


## function: cleanCat(outPDBname,currPDBfname,refPDBfname,[newCid])
## concatenates both PDBs into a single file,
## separating them with "MODEL" fields (for rasmol viewing)
## Note: chain ID of reference will be changed

sub cleanCat {
  my $outPDBfname=shift;
  my $currPDBfname=shift;
  my $refPDBfname=shift;
  my $newCid=shift;

  my $Cmol=Molecule->new($currPDBfname);
  my $Rmol=Molecule->new($refPDBfname);

  my $Chains_used={};

  if (defined $newCid) {
      my @Rchains=@{$Rmol->activeChains()};
      foreach my $Rc ( @Rchains ) {
	  $Rmol->selectChain($Rc->{id});
	  $Rmol->setChain($newCid);
      }

  } else {

      # Choose new chain IDs randomly; do this by picking
      # letters until we find one that's not already taken
      foreach my $Cc ( @{$Cmol->activeChains()} ) {
	  $Chains_used->{$Cc->{id}}=1;
      }
      my @Rchains=@{$Rmol->activeChains()};
      foreach my $Rc ( @Rchains ) {
	  if (defined $Chains_used->{$Rc->{id}}) {
	      $Rmol->selectChain($Rc->{id});
	      my $newCid=chr(int(rand(26)+65));
	      while (defined $Chains_used->{$newCid}) {
		  $newCid=chr(int(rand(26)+65));
	      }
	      $Chains_used->{$newCid}=1;
	      $Rmol->setChain($newCid);
	  }
      }
  }

  my $Cname=lc(hostname."-Cpdb$$");
  my $Rname=lc(hostname."-Rpdb$$");
  $Cmol->writePDB($Cname,"generic");
  $Rmol->writePDB($Rname,"generic");

  my $outPDB=&GenUtil::getOutputFile($outPDBfname);

  printf $outPDB "MODEL 1\n";
  my $inC=&GenUtil::getInputFile($Cname);
  while (<$inC>) {
      if (! /^END/) {
	  chomp;
	  printf $outPDB "%s\n",$_;
      }
  }
  undef $inC;
  printf $outPDB "ENDMDL 1\n";

  printf $outPDB "MODEL 2\n";
  my $inR=&GenUtil::getInputFile($Rname);
  while (<$inR>) {
      if (! /^END/) {
	  chomp;
	  printf $outPDB "%s\n",$_;
      }
  }
  undef $inR;
  printf $outPDB "ENDMDL 2\n";
  printf $outPDB "END\n";

  undef $outPDB;
  &GenUtil::remove($Cname);
  &GenUtil::remove($Rname);
  return;
}


1;
