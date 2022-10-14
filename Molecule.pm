# Molecule package
# read/write/convert structure info
#
# http://mmtsb.scripps.edu/doc/Molecule.pm.html
# 2000, Michael Feig, Brooks group, TSRI
# 2002-2005, John Karanicolas, Baker lab, UW
#

package Molecule;

require 5.004;

use strict;
no warnings "redefine";

use FileHandle;
use IPC::Open2;

use GenUtil;
use Sequence;
use SICHO;
use Analyze;

## data: chain[] -> { id atom[] res[] resinx[]
## data:              xcoor[] ycoor[] zcoor[] }
## data:   atom[] -> { atominx atomname resname resnum chain
## data:               xcoor ycoor zcoor hyd bb seg aux1 aux2
## data:               aux3 }
## data:   res[]  -> { name num start end seg valid chain }
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

## data: info -> { title[] header[] }
## info from PDB header

## constructor: new([PDBfile])
## creates a new Molecule object and reads a PDB structures
## if a file name is given

## constructor: new([PDBfile])
## creates a new Molecule object and reads a PDB structures
## if a file name is given
## (this form can act as a base class for derived classes)

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $farg = shift;

  my $self = {};
  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{segmentlist}=undef;
  $self->{ssbond}=();
  $self->{info}={};

  bless($self,$class);

  $self->readPDB($farg) 
    if (defined $farg);

  return $self;
}

## method: readPDB(file[,translate[,ignoreseg]])
## reads a protein structures from a PDB file.
## <mark>translate</mark> may be set to <mark>CHARMM19</mark>
## for proper recognition of histidine residues.
## If <mark>ignoreseg</mark> is set segment IDs from the
## PDB are not read.

sub readPDB {
  my $self=shift;
  my $fname=&GenUtil::getInputFile(shift);
  my %par=@_;

  my $translate=(defined $par{translate})?$par{translate}:"GENERIC";
  my $ignoreseg=(defined $par{ignoreseg} && $par{ignoreseg})?1:0;
  my $chainfromseg=(defined $par{chainfromseg} && $par{chainfromseg})?1:0;
  my $ignorehet=(defined $par{ignorehet} && $par{ignorehet})?1:0;
  my $ignoressbond=(defined $par{ignoressbond} && $par{ignoressbond})?1:0;
  my $activemodel;

  $translate=uc($translate);

  my $lastchain=".";
  my $chainrec;

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{ssbond}=();

  my $lastnum=-999;
  my $lastinum=-999;
  my $ignore;

  my $newchain=0;
  my $got_first_atom=0;

  # Clear the header info
  $self->{info}={};

 READPDB:
  while(<$fname>) {

    # don't include lines that start with these as part of the header
    if ((! /^(SSBOND|ATOM|HETATM|TER|END|COMPND|MODEL|KEYWDS|EXPDTA|AUTHOR|REVDAT|JRNL|REMARK)/) &&
	(! /^(DBREF|SEQADV|SEQRES|HET|FORMUL|HELIX|SHEET|TURN|LINK|SITE|CRYST|ORIG|SCALE|CONECT|MASTER|ANISOU)/)) {
	if (! $got_first_atom) {
	    chomp;
	    my @p=split(' ',$_);
	    my $k=shift(@p);
	    if (! defined $self->{info}->{$k}) {
		$self->{info}->{$k}=();
		push(@{$self->{info}->{$k}},join(" ",@p));
	    } else {
		shift(@p);
		push(@{$self->{info}->{$k}},join(" ",@p));
	    }
	}

    } elsif ((/^SSBOND/) && (! $ignoressbond)) {
      my $trec={};
      ($trec->{chain1}=substr($_,15,1))=~s/ //g;
      ($trec->{resnum1}=substr($_,16,5))=~s/ //g;
      ($trec->{chain2}=substr($_,29,1))=~s/ //g;
      ($trec->{resnum2}=substr($_,30,5))=~s/ //g;
      push(@{$self->{ssbond}},$trec);
    } elsif (/^MODEL +([0-9]+)/) {
      if (!defined $par{model} || $1 == $par{model}) {
	$activemodel=1;
      } else {
	$activemodel=0;
      }
    }

    # jk allow reading past "END" but not "ENDMDL"
    my $lc_line = lc($_);
    last READPDB if ($lc_line =~ /^(endmdl|\#endmdl)/ && (!defined $activemodel || $activemodel));

    if (/^ATOM/ && (!defined $activemodel || $activemodel) && (! /HOH/)) {
      my ($atomname, $resname, $resnum, $iresnum, $alt, $chain, $seg);

      if (($alt=substr($_,16,1))=~/[ A0-9]/) {
	($atomname=substr($_,12,4))=~s/ //g;
	$atomname.=$alt if ($alt=~/[0-9]/);
	($resname=substr($_,17,4))=~s/ //g;
	($resnum=substr($_,22,5))=~s/ //g;

	if ($resnum =~ /[A-Z]+/) {
	    ($iresnum=$resnum)=~s/([A-Z])+//g;
            # Change the sign of this residue number, apply an offset
	    $iresnum+=ord($1);
	    $iresnum*=-1;
	} else {
	    $iresnum=$resnum;
	}

	$got_first_atom=1;
	$iresnum+=0;
	$chain=substr($_,21,1);
	($seg=substr($_,72,4))=~s/[ \n]//g;

	if ($chain eq " " && $chainfromseg && !$ignoreseg && $seg=~/...([A-Z\+])/) {
	  $chain=$1;
	}

	if ($resnum ne $lastnum) {
	  $ignore=($iresnum == $lastinum);
	}

	if (!$ignore) {
	  $atomname="CD1" 
	    if ($resname eq "ILE" && $atomname eq "CD");
	  $atomname="O"
	    if (($atomname eq "OT1" || $atomname eq "O1" || $atomname eq "OCT1")
		&& $resname =~/ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HSP|HIS|HID|HIE|HIP|HSP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX/);
	  $atomname="OXT"
	    if (($atomname eq "OT2" || $atomname eq "O2" || $atomname eq "OCT2")
		&& $resname =~/ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HSP|HIS|HSP|HID|HIE|HIP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX/);
	  
	  if ($translate eq "CHARMM19") {
	    $resname=~s/HSD/HSE/;
	    $resname=~s/HIS/HSD/;
          } elsif ($translate eq "AMBER") {
	    $resname=~s/CYX/CYS/;
	    $resname=~s/HID/HSD/;
	    $resname=~s/HIE/HSE/;
	    $resname=~s/HIP/HSP/;
            $resname=~s/ASH/ASP/;
            $resname=~s/GLH/GLU/;
            $resname=~s/CTG/GLY/;

	    $atomname="HT1" if ($atomname eq "H1");
	    $atomname="HT2" if ($atomname eq "H2");
	    $atomname="HT3" if ($atomname eq "H3");
	    $atomname="HN" if ($atomname eq "H");
	    $atomname="HB1" if ($atomname eq "HB2" && $resname=~/MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HSD|HSE|HSP/);
	    $atomname="HB2" if ($atomname eq "HB3" && $resname=~/MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HSD|HSE|HSP/);
	    $atomname="HG1" if ($atomname eq "HG2" && $resname=~/MET|GLU|GLN|LYS|ARG|PRO/);
	    $atomname="HG2" if ($atomname eq "HG3" && $resname=~/MET|GLU|GLN|LYS|ARG|PRO/);
	    $atomname="HG1" if ($atomname eq "HSG" && $resname=~/CYS/);
	    $atomname="HG1" if ($atomname eq "HG" && $resname=~/CYS/);

	    $atomname="HD1" if ($atomname eq "HD2" && $resname=~/LYS|ARG|PRO/);
	    $atomname="HD2" if ($atomname eq "HD3" && $resname=~/LYS|ARG|PRO/);
	    $atomname="HE1" if ($atomname eq "HE2" && $resname=~/LYS/);
	    $atomname="HE2" if ($atomname eq "HE3" && $resname=~/LYS/);
	    $atomname="HA1" if ($atomname eq "HA2" && $resname=~/GLY/);
	    $atomname="HA2" if ($atomname eq "HA3" && $resname=~/GLY/);
	    $atomname="HG1"  if ($atomname eq "HG" && $resname=~/SER/);
	    $atomname="HD1" if ($atomname eq "HD11" && $resname=~/ILE/);
	    $atomname="HD2" if ($atomname eq "HD12" && $resname=~/ILE/);
	    $atomname="HD3" if ($atomname eq "HD13" && $resname=~/ILE/);
	    $atomname="HG11" if ($atomname eq "HG12" && $resname=~/ILE/);
	    $atomname="HG12" if ($atomname eq "HG13" && $resname=~/ILE/);
	    $atomname="HN1" if ($atomname eq "H2" && $resname=~/PRO/);
	    $atomname="HN2" if ($atomname eq "H3" && $resname=~/PRO/);
          } elsif ($translate eq "IMPACT") {
	    $resname="HSD" if ($resname eq "HID");
	    $resname="HSE" if ($resname eq "HIE");
	    $atomname="HT1" if ($atomname eq "H1");
	    $atomname="HT2" if ($atomname eq "H2");
	    $atomname="HT3" if ($atomname eq "H3");
	    $atomname="HN" if ($atomname eq "H");
            $atomname="HN1" if ($atomname eq "2H" && $resname eq "PRO");
            $atomname="HN2" if ($atomname eq "3H" && $resname eq "PRO");
	    $atomname="HA1" if ($atomname eq "1HA");
	    $atomname="HA2" if ($atomname eq "2HA");
	    $atomname="HB1" if ($atomname eq "1HB");
	    $atomname="HB2" if ($atomname eq "2HB");
	    $atomname="HB3" if ($atomname eq "3HB");
	    $atomname="HG1" if ($atomname eq "1HG");
	    $atomname="HG2" if ($atomname eq "2HG");
	    $atomname="HE1" if ($atomname eq "1HE");
	    $atomname="HE2" if ($atomname eq "2HE");
	    $atomname="HE3" if ($atomname eq "2HE");

	    $atomname="HD1" if ($atomname eq "1HD1" && $resname eq "ILE");
	    $atomname="HD2" if ($atomname eq "2HD1" && $resname eq "ILE");
	    $atomname="HD3" if ($atomname eq "3HD1" && $resname eq "ILE");

	    $atomname="HD11" if ($atomname eq "1HD1");
	    $atomname="HD12" if ($atomname eq "2HD1");
	    $atomname="HD13" if ($atomname eq "3HD1");
	    $atomname="HD21" if ($atomname eq "1HD2");
	    $atomname="HD22" if ($atomname eq "2HD2");
	    $atomname="HD23" if ($atomname eq "3HD2");
	    $atomname="HD1" if ($atomname eq "1HD");
	    $atomname="HD2" if ($atomname eq "2HD");
	    $atomname="HD3" if ($atomname eq "2HD");
	    $atomname="HZ1" if ($atomname eq "1HZ");
	    $atomname="HZ2" if ($atomname eq "2HZ");
	    $atomname="HZ3" if ($atomname eq "2HZ");
	    $atomname="HG11" if ($atomname eq "1HG1");
	    $atomname="HG12" if ($atomname eq "2HG1");
	    $atomname="HG13" if ($atomname eq "3HG1");
	    $atomname="HG21" if ($atomname eq "1HG2");
	    $atomname="HG22" if ($atomname eq "2HG2");
	    $atomname="HG23" if ($atomname eq "3HG2");
	    $atomname="HH11" if ($atomname eq "1HH1");
	    $atomname="HH12" if ($atomname eq "2HH1");
	    $atomname="HH21" if ($atomname eq "1HH2");
	    $atomname="HH22" if ($atomname eq "2HH2");
	    $atomname="HG1"  if ($atomname eq "HG" && $resname=~/SER|CYS/);
	  } else { 
#	    $resname=~s/CYX/CYS/;
#	    $resname=~s/HID/HSD/;
#	    $resname=~s/HIE/HSE/;
#	    $resname=~s/HIP/HSP/;
	  }

	  if ($chain ne $lastchain) {
            my $crec=$self->{chainlookup}->{$chain};
	    $chainrec=(defined $crec)?$crec:$self->_newChain($chain);
	    $newchain=1;
	  } else {
	    $newchain=0;
	  }

	  my $pdbrec={};

	  $pdbrec->{atominx}=substr($_,4,7)+0;
	  $pdbrec->{atomname}=$atomname;
	  $pdbrec->{resname}=$resname;
	  $pdbrec->{resnum}=$iresnum;

	  $pdbrec->{chain}=$chainrec->{id};
	  $pdbrec->{xcoor}=substr($_,30,8)+0.0;
	  $pdbrec->{ycoor}=substr($_,38,8)+0.0;
	  $pdbrec->{zcoor}=substr($_,46,8)+0.0;
	  $pdbrec->{hyd}=($atomname=~/^[0-9]*H.*/)?1:0;
	  $pdbrec->{bb}=1 if (($atomname eq "C") || ($atomname eq "O") ||
			      ($atomname eq "N") || ($atomname eq "H") ||
			      ($atomname eq "HN") || ($atomname eq "CA"));
	  $pdbrec->{aux1}=substr($_,55,6)+0.0;
	  $pdbrec->{aux2}=substr($_,61,6)+0.0;

	  $pdbrec->{seg}=$seg
	    unless (defined $ignoreseg && $ignoreseg);
	  
	  push (@{$chainrec->{atom}}, $pdbrec);

#	  printf "%s %s %d %s %s %s %s\n",$chainrec->{id},$resnum,$iresnum,$lastnum,$resname,$seg,$atomname;

	  if ($resnum ne $lastnum || $newchain) {
	    my $resrec={};
	    $resrec->{name}=$resname;
	    $resrec->{num}=$iresnum;
	    $resrec->{chain}=$chainrec->{id};
	    $resrec->{start}=$#{$chainrec->{atom}};
	    $resrec->{end}=$resrec->{start};
	    $resrec->{valid}=1;
	    $resrec->{seg}=$seg
	      unless (defined $ignoreseg && $ignoreseg);
	    push(@{$chainrec->{res}},$resrec);
	  } else {
	    $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
	  }

	}
	$lastnum=$resnum;
	$lastinum=$iresnum;
	$lastchain=$chain;
      }
    } elsif (/^HETATM/ && (! $ignorehet) &&
              (!defined $activemodel || $activemodel) && (! /WSS/)) {
      my ($atomname, $resname, $resnum, $iresnum, $alt, $chain, $seg);

      ($atomname=substr($_,12,4))=~s/ //g;
      $atomname.=$alt if ($alt=~/[0-9]/);
      ($resname=substr($_,17,4))=~s/ //g;
      ($resnum=substr($_,22,5))=~s/ //g;

      if ($resnum =~ /[A-Z]+/) {
	  ($iresnum=$resnum)=~s/([A-Z])+//g;
          # Change the sign of this residue number, apply an offset
	  $iresnum+=ord($1);
	  $iresnum*=-1;
      } else {
	  $iresnum=$resnum;
      }
      $iresnum+=0;
      ($seg=substr($_,72,4))=~s/[ \n]//g;

      $chain="+";

      my $crec=$self->{chainlookup}->{$chain};
      if (defined $crec) {
        $chainrec=$crec;
        $newchain=0;
      } else {
        $chainrec=$self->_newChain($chain);
	$newchain=1;
      }

      my $pdbrec={};
	
      $pdbrec->{atominx}=substr($_,6,5)+0;
      $pdbrec->{atomname}=$atomname;
      $pdbrec->{resname}=$resname;
      $pdbrec->{resnum}=$iresnum;

      $pdbrec->{chain}=$chainrec->{id};
      $pdbrec->{xcoor}=substr($_,30,8)+0.0;
      $pdbrec->{ycoor}=substr($_,38,8)+0.0;
      $pdbrec->{zcoor}=substr($_,46,8)+0.0;
      $pdbrec->{hyd}=($atomname=~/^[0-9]*H.*/)?1:0;
      $pdbrec->{bb}=1 if (($atomname eq "C") || ($atomname eq "O") ||
			  ($atomname eq "N") || ($atomname eq "H") ||
			  ($atomname eq "HN") || ($atomname eq "CA"));
      $pdbrec->{aux1}=substr($_,55,6)+0.0;
      $pdbrec->{aux2}=substr($_,61,6)+0.0;
      $pdbrec->{seg}=$seg
	unless (defined $ignoreseg && $ignoreseg);
      
      push (@{$chainrec->{atom}}, $pdbrec);

#      printf "HA %s %d %d %s %s %s\n",$chainrec->{id},$resnum,$lastnum,$resname,$seg,$atomname;
      
      if ($resnum ne $lastnum || $newchain) {
	my $resrec={};
	$resrec->{name}=$resname;
	$resrec->{num}=$iresnum;
	$resrec->{chain}=$chainrec->{id};
	$resrec->{start}=$#{$chainrec->{atom}};
	$resrec->{end}=$resrec->{start};
	$resrec->{valid}=1;
	$resrec->{seg}=$seg
	  unless (defined $ignoreseg && $ignoreseg);
	push(@{$chainrec->{res}},$resrec);
      } else {
	$chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
      }
      $lastnum=$resnum;
      $lastinum=$iresnum;
    }
  }

  $self->_correctHisType();

  $self->{segmentlist}=undef;

  $self->_coorCache();

  undef $fname;
}

## method: readUnformattedPDB(file[,translate[,ignoreseg]])
## reads a protein structures from a PDB file,
## renumbering across chain breaks and merging across chain IDs and NMR models,
## reading across ENDs and TERs

sub readUnformattedPDB {
  my $self=shift;
  my $fname=&GenUtil::getInputFile(shift);

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{ssbond}=();

  my $last_inresnum=-999;
  my $outresnum=0;
  my $chainrec;

  my $got_first_atom=0;

  # Clear the header info
  $self->{info}={};

  while(<$fname>) {

    # don't include lines that start with these as part of the header
    if ((! /^(SSBOND|ATOM|HETATM|TER|END|COMPND|MODEL|KEYWDS|EXPDTA|AUTHOR|REVDAT|JRNL|REMARK)/) &&
	(! /^(DBREF|SEQADV|SEQRES|HET|FORMUL|HELIX|SHEET|TURN|LINK|SITE|CRYST|ORIG|SCALE|CONECT|MASTER|ANISOU)/)) {
      if (! $got_first_atom) {
	chomp;
	my @p=split(' ',$_);
	my $k=shift(@p);
	if (! defined $self->{info}->{$k}) {
	  $self->{info}->{$k}=();
	  push(@{$self->{info}->{$k}},join(" ",@p));
	} else {
	  shift(@p);
	  push(@{$self->{info}->{$k}},join(" ",@p));
	}
      }

    } else {
      my $inline = $_;

      # Replace seleno-Met with Met and methyl-Lys with Lys
      if ( $inline =~ / MSE / ) {
	$inline=~s/MSE/MET/;
	$inline=~s/HETATM/ATOM  /;
      } elsif ( $inline =~ / MLY / ) {
	$inline=~s/MLY/LYS/;
	$inline=~s/HETATM/ATOM  /;
      }

      if (/^ATOM/ && (! /HOH/)) {

	my ($atomname, $resname, $resnum, $iresnum, $chain, $seg);
	($atomname=substr($_,12,4))=~s/ //g;
	($resname=substr($_,17,4))=~s/ //g;
	($resnum=substr($_,22,5))=~s/ //g;

	if ($resnum =~ /[A-Z]+/) {
	  ($iresnum=$resnum)=~s/([A-Z])+//g;
	  # Change the sign of this residue number, apply an offset
	  $iresnum+=ord($1);
	  $iresnum*=-1;
	} else {
	  $iresnum=$resnum;
	}
	$iresnum+=0;

	my $newres=0;
	if ( $last_inresnum != $iresnum ) {
	  $last_inresnum = $iresnum;
	  $outresnum++;
	  $newres=1;
	}
	$iresnum = $outresnum;
	$chain="A";

	if ( ! $got_first_atom ) {
	  $got_first_atom=1;
	  my $crec=$self->{chainlookup}->{$chain};
	  $chainrec=(defined $crec)?$crec:$self->_newChain($chain);
	}

	my $pdbrec={};

	$pdbrec->{atominx}=substr($inline,4,7)+0;
	$pdbrec->{atomname}=$atomname;
	$pdbrec->{resname}=$resname;
	$pdbrec->{resnum}=$iresnum;

	$pdbrec->{chain}=$chainrec->{id};
	$pdbrec->{xcoor}=substr($_,30,8)+0.0;
	$pdbrec->{ycoor}=substr($_,38,8)+0.0;
	$pdbrec->{zcoor}=substr($_,46,8)+0.0;
	$pdbrec->{hyd}=($atomname=~/^[0-9]*H.*/)?1:0;
	$pdbrec->{bb}=1 if (($atomname eq "C") || ($atomname eq "O") ||
			    ($atomname eq "N") || ($atomname eq "H") ||
			    ($atomname eq "HN") || ($atomname eq "CA"));
	$pdbrec->{aux1}=substr($_,55,6)+0.0;
	$pdbrec->{aux2}=substr($_,61,6)+0.0;

	push (@{$chainrec->{atom}}, $pdbrec);

#	  printf "%s %s %d %s %s %s %s\n",$chainrec->{id},$resnum,$iresnum,$lastnum,$resname,$seg,$atomname;

	if ( $newres ) {
	  my $resrec={};
	  $resrec->{name}=$resname;
	  $resrec->{num}=$iresnum;
	  $resrec->{chain}=$chainrec->{id};
	  $resrec->{start}=$#{$chainrec->{atom}};
	  $resrec->{end}=$resrec->{start};
	  $resrec->{valid}=1;
	  push(@{$chainrec->{res}},$resrec);
	} else {
	  $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
	}

      }
    }
  }

  $self->_correctHisType();

  $self->{segmentlist}=undef;

  $self->_coorCache();

  undef $fname;
}

## method: readCRD(file)
## reads a protein structures from a CHARMM CRD file.

sub readCRD {
  my $self=shift;
  my $fname=&GenUtil::getInputFile(shift);

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;

  my $lastchain=".";
  my $chainrec;

  my $lastnum=-999;
  my $ignore;

  my $newchain=0;
  
  my $first=1;

  while(<$fname>) {
    if (!/^\*/) {
      if ($first) {
	$first=0;
      } else {
	my ($atomname, $resname, $resnum, $iresnum, $chain,$seg);
      
	($atomname=substr($_,16,4))=~s/ //g;
	($resname=substr($_,11,4))=~s/ //g;
	$resnum=substr($_,5,5);
	$iresnum=$resnum+0;
	($seg=substr($_,51,4))=~s/[ \n]//g;

	$chain=($resname eq "TIP3" || $resname eq "HOH")?"+":" ";
	  
	$atomname="CD1" 
	  if ($resname eq "ILE" && $atomname eq "CD");
	$atomname="O"
	  if (($atomname eq "OT1" || $atomname eq "O1" || $atomname eq "OCT1")
	      && $resname =~/ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HSP|HIS|HSP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX/);
	$atomname="OXT"
	  if (($atomname eq "OT2" || $atomname eq "O2" || $atomname eq "OCT2")
	      && $resname =~/ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HSP|HIS|HSP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX/);
	
	$resname=~s/HSD/HSE/;
	$resname=~s/HIS/HSD/;
	
	if ($chain ne $lastchain) {
	  my $crec=$self->{chainlookup}->{$chain};
	  $chainrec=(defined $crec)?$crec:$self->_newChain($chain);
	  $newchain=1;
	} else {
	  $newchain=0;
	}
	
	my $pdbrec={};
	
	$pdbrec->{atominx}=substr($_,0,5)+0;
	$pdbrec->{atomname}=$atomname;
	$pdbrec->{resname}=$resname;
	$pdbrec->{resnum}=$resnum+0;
	
	$pdbrec->{chain}=$chainrec->{id};
	$pdbrec->{xcoor}=substr($_,20,10)+0.0;
	$pdbrec->{ycoor}=substr($_,30,10)+0.0;
	$pdbrec->{zcoor}=substr($_,40,10)+0.0;
	$pdbrec->{hyd}=($atomname=~/^[0-9]*H.*/)?1:0;
	$pdbrec->{bb}=1 if (($atomname eq "C") || ($atomname eq "O") ||
			    ($atomname eq "N") || ($atomname eq "H") ||
			    ($atomname eq "HN") || ($atomname eq "CA"));
	$pdbrec->{seg}=$seg;
	  
	push (@{$chainrec->{atom}}, $pdbrec);

#	printf "%s %s %d %s %s %s %s\n",$chainrec->{id},$resnum,$iresnum,$lastnum,$resname,$seg,$atomname;
	
	if ($iresnum != $lastnum || $newchain) {
	  my $resrec={};
	  $resrec->{name}=$resname;
	  $resrec->{num}=$iresnum;
	  $resrec->{chain}=$chainrec->{id};
	  $resrec->{start}=$#{$chainrec->{atom}};
	  $resrec->{end}=$resrec->{start};
	  $resrec->{valid}=1;
	  $resrec->{seg}=$seg;
	  push(@{$chainrec->{res}},$resrec);
	} else {
	  $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
	}
	$lastnum=$iresnum;
	$lastchain=$chain;
      }
    }
  }

  $self->{segmentlist}=undef;

  $self->_coorCache();

  undef $fname;
}

## method: writePDB(file[,translate=>format, longaux2=>1, ssbond=>0])
## writes out the current structure in PDB format
## the output format may be specified through <mark>translate</mark>.
## Possible formats are <mark>CHARMM19</mark>, <mark>CHARMM22</mark>,
## <mark>AMBER</mark>, <mark>GENERIC</mark>, and <mark>GENERIC_NOH</mark>.
## A segment ID as required by CHARMM may be given as the third argument

sub writePDB {
  my $self=shift;
  my $PDBname=shift;
  my %wpar=@_;

  my $fname=&GenUtil::getOutputFile($PDBname);

  my $translate=$wpar{translate};
  my $longaux2=$wpar{longaux2};
  my $ssbond=$wpar{ssbond};
  my $reqHis=$wpar{reqHis};
  my $compress=$wpar{compress};
  my $allow_rosetta_water=$wpar{allow_rosetta_water};

  $translate="GENERIC" if (!defined $translate);
  $reqHis="" if (!defined $reqHis);
  $compress=1 if (!defined $compress);
  $allow_rosetta_water=1 if (!defined $allow_rosetta_water);

  die "empty molecule, nothing to write"
    if (!defined $self->{chain} || $#{$self->{chain}}<0);

  # These lines will go at the top, followed by "MODULE", 
  # followed by all other info lines (order unimportant)
  my @startFields=();
  push(@startFields,"HEADER");
  push(@startFields,"TITLE");
  push(@startFields,"MODULE");

  foreach my $k ( @startFields ) {
      if ( ( defined $self->{info}->{$k} ) && ($#{$self->{info}->{$k}} >= 0)) {
	  my $a=$self->{info}->{$k};
	  my $lastLineInx=$#{$a};
	  # The first line has extra spaces and no line number
	  printf $fname "%-9s %s\n",$k,$a->[0];
	  for (my $i=1; $i<=$lastLineInx; $i++) {
	      printf $fname "%-6s %3d %s\n",$k,$i+1,$a->[$i];
	  }
      }
  }

  foreach my $k ( keys %{$self->{info}} ) {
      my $already_printed=0;
      foreach my $nk ( @startFields ) {
	  $already_printed=1 if ($k eq $nk);
      }
      if ( ! $already_printed ) {
	  my $a=$self->{info}->{$k};
	  my $lastLineInx=$#{$a};
	  # The first line has extra spaces and no line number
	  printf $fname "%-9s %s\n",$k,$a->[0];
	  for (my $i=1; $i<=$lastLineInx; $i++) {
	      printf $fname "%-6s %3d %s\n",$k,$i+1,$a->[$i];
	  }
      }
  }

  if (!defined $ssbond || $ssbond) {
    my %havess;
    my $sinx=0;
    foreach my $s ( @{$self->{ssbond}}) {
      my $c1=$self->getChain($s->{chain1});
      my $c2=$self->getChain($s->{chain2});
      if (defined $c1 && defined $c2) {
	my $r1=$self->getResidueInChain($s->{resnum1},$c1);
	my $r2=$self->getResidueInChain($s->{resnum2},$c2);
	if (defined $r1 && defined $r2 && 
	    $r1->{name} eq "CYS" && $r2->{name} eq "CYS") {
	  printf $fname "SSBOND%4d CYS %1s%5d    CYS %1s%5d\n",
	    ++$sinx,$s->{chain1},$s->{resnum1},$s->{chain2},$s->{resnum2}
	      unless ($havess{"$s->{chain1}:$s->{resnum1}:$s->{chain2}:$s->{resnum2}"});

	  $havess{"$s->{chain1}:$s->{resnum1}:$s->{chain2}:$s->{resnum2}"}=1;
	  $havess{"$s->{chain2}:$s->{resnum2}:$s->{chain1}:$s->{resnum1}"}=1;
	}
      }
    }
  }

  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{atom}}>=0) {
      my $cterm=0;
      my $lastres=$c->{res}->[$#{$c->{res}}];
      for (my $i=$lastres->{start}; $i<=$lastres->{end}; $i++) {
	if ($c->{atom}->[$i]->{atomname} eq "OXT") {
	  $cterm=1;
	}
      }

      my $prevres;
      my $HisType="";
      foreach my $a ( @{$c->{atom}} ) {
	my $ta;
	%{$ta}=%{$a};

	if ($translate=~/CHA/) {
	  $ta->{atomname}="CD" 
	    if ($ta->{resname} eq "ILE" && $ta->{atomname} eq "CD1");
	  $ta->{atomname}="OT1"
	    if ($ta->{atomname} eq "O" && $ta->{resnum} eq $lastres->{num} && $cterm);
	  $ta->{atomname}="OT2"
	    if ($ta->{atomname} eq "OXT" && $ta->{resnum} eq $lastres->{num});
	  $ta->{resname}=~s/HOH/TIP3/;
	} 

	if (($ta->{resname}=~/HSD|HSE|HIS|HSP/) && ($translate =~ /CHARMM/)) {
	    my $hisName=$ta->{resname};
	    if ($translate =~ /CHARMM19/) {
		if ($hisName eq "HSD") {
		    $hisName="HIS";
		} elsif ($hisName eq "HSE") {
		    $hisName="HSD";
		} elsif ($hisName eq "HSP") {
		    $hisName="HSC";
		}
	    } elsif (($hisName eq "HIS") && ($translate =~ /CHARMM22/)) {
		$hisName="HSD";
	    }
	    $ta->{resname}=$hisName;

	} elsif ($translate =~ /AMBER/) {
	  $ta->{resname}=~s/HSD/HID/;
	  $ta->{resname}=~s/HSE/HIE/;
	  $ta->{resname}=~s/HSP/HIP/;
#	  $ta->{atomname}="CD1" 
#	    if ($ta->{resname} eq "ILE" && $ta->{atomname} eq "CD");

	  if ($translate !~ /CHAMBER/) {
	    $ta->{atomname}="H1" if ($ta->{atomname} eq "HT1");
	    $ta->{atomname}="H2" if ($ta->{atomname} eq "HT2");
	    $ta->{atomname}="H3" if ($ta->{atomname} eq "HT3");
	    $ta->{atomname}="H" if ($ta->{atomname} eq "HN");
	    $ta->{atomname}="HB3" if ($ta->{atomname} eq "HB2" && $ta->{resname}=~/MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HID|HIE/);
	    $ta->{atomname}="HB2" if ($ta->{atomname} eq "HB1" && $ta->{resname}=~/MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HID|HIE/);
	    $ta->{atomname}="HG3" if ($ta->{atomname} eq "HG2" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	    $ta->{atomname}="HG2" if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	    $ta->{atomname}=(($ta->{resnum} eq $c->{atom}->[0]->{resnum}) || ($ta->{resnum} eq $lastres->{num} && $cterm))
	      ?"HSG":"HG" if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/CYS/);
	    $ta->{atomname}="HD3" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/LYS|ARG|PRO/);
	    $ta->{atomname}="HD2" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/LYS|ARG|PRO/);
	    $ta->{atomname}="HE3" if ($ta->{atomname} eq "HE2" && $ta->{resname}=~/LYS/);
	    $ta->{atomname}="HE2" if ($ta->{atomname} eq "HE1" && $ta->{resname}=~/LYS/);
	    $ta->{atomname}="HA3" if ($ta->{atomname} eq "HA2" && $ta->{resname}=~/GLY/);
	    $ta->{atomname}="HA2" if ($ta->{atomname} eq "HA1" && $ta->{resname}=~/GLY/);
	    $ta->{atomname}="HG"  if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/SER/);
	    $ta->{atomname}="HD11" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/ILE/);
	    $ta->{atomname}="HD12" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/ILE/);
	    $ta->{atomname}="HD13" if ($ta->{atomname} eq "HD3" && $ta->{resname}=~/ILE/);
	    $ta->{atomname}="HG13" if ($ta->{atomname} eq "HG12" && $ta->{resname}=~/ILE/);
	    $ta->{atomname}="HG12" if ($ta->{atomname} eq "HG11" && $ta->{resname}=~/ILE/);
	    $ta->{atomname}="H2" if ($ta->{atomname} eq "HN1" && $ta->{resname}=~/PRO/);
	    $ta->{atomname}="H3" if ($ta->{atomname} eq "HN2" && $ta->{resname}=~/PRO/);
	  } 
	} elsif ($translate =~ /IMPACT/) {
	  $ta->{resname}="HID" if ($ta->{resname} =~/HSD/);
	  $ta->{resname}="HIE" if ($ta->{resname} =~/HSE/);
	  $ta->{atomname}="1H" if ($ta->{atomname} eq "HT1");
	  $ta->{atomname}="2H" if ($ta->{atomname} eq "HT2");
	  $ta->{atomname}="3H" if ($ta->{atomname} eq "HT3");
	  $ta->{atomname}="H" if ($ta->{atomname} eq "HN");
	  $ta->{atomname}="2H" if ($ta->{atomname} eq "HN1");
	  $ta->{atomname}="3H" if ($ta->{atomname} eq "HN2");
	  $ta->{atomname}="3HB" if ($ta->{atomname} eq "HB3" && $ta->{resname}=~/ALA/);
	  $ta->{atomname}="2HB" if ($ta->{atomname} eq "HB2" && $ta->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $ta->{atomname}="1HB" if ($ta->{atomname} eq "HB1" && $ta->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $ta->{atomname}="2HG" if ($ta->{atomname} eq "HG2" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $ta->{atomname}="1HG" if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $ta->{atomname}="2HA" if ($ta->{atomname} eq "HA2" && $ta->{resname}=~/GLY/);
	  $ta->{atomname}="1HA" if ($ta->{atomname} eq "HA1" && $ta->{resname}=~/GLY/);
	  $ta->{atomname}="3HE" if ($ta->{atomname} eq "HE3" && $ta->{resname}=~/MET/);
	  $ta->{atomname}="2HE" if ($ta->{atomname} eq "HE2" && $ta->{resname}=~/MET|LYS/);
	  $ta->{atomname}="1HE" if ($ta->{atomname} eq "HE1" && $ta->{resname}=~/MET|LYS/);
	  $ta->{atomname}="2HD" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/ARG|PRO/);
	  $ta->{atomname}="1HD" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/ARG|PRO/);
	  $ta->{atomname}="1HH1" if ($ta->{atomname} eq "HH11" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="2HH1" if ($ta->{atomname} eq "HH12" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="1HH2" if ($ta->{atomname} eq "HH21" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="2HH2" if ($ta->{atomname} eq "HH22" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="1HE2" if ($ta->{atomname} eq "HE21" && $ta->{resname}=~/GLN/);
	  $ta->{atomname}="2HE2" if ($ta->{atomname} eq "HE22" && $ta->{resname}=~/GLN/);
	  $ta->{atomname}="1HD1" if ($ta->{atomname} eq "HD11" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="2HD1" if ($ta->{atomname} eq "HD12" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="3HD1" if ($ta->{atomname} eq "HD13" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="1HD2" if ($ta->{atomname} eq "HD21" && $ta->{resname}=~/LEU|ASN/);
	  $ta->{atomname}="2HD2" if ($ta->{atomname} eq "HD22" && $ta->{resname}=~/LEU|ASN/);
	  $ta->{atomname}="3HD2" if ($ta->{atomname} eq "HD23" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="1HD" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/LYS|ARG|PRO/);
	  $ta->{atomname}="2HD" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/LYS|ARG|PRO/);
	  $ta->{atomname}="3HZ" if ($ta->{atomname} eq "HZ3" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="2HZ" if ($ta->{atomname} eq "HZ2" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="1HZ" if ($ta->{atomname} eq "HZ1" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="1HG1" if ($ta->{atomname} eq "HG11" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="2HG1" if ($ta->{atomname} eq "HG12" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="3HG1" if ($ta->{atomname} eq "HG13" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="1HG2" if ($ta->{atomname} eq "HG21" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="2HG2" if ($ta->{atomname} eq "HG22" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="3HG2" if ($ta->{atomname} eq "HG23" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="HG"  if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/SER|CYS/);
	  $ta->{atomname}="1HD1" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HD1" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="3HD1" if ($ta->{atomname} eq "HD3" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="1HG1" if ($ta->{atomname} eq "HG11" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HG1" if ($ta->{atomname} eq "HG12" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="1HG2" if ($ta->{atomname} eq "HG21" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HG2" if ($ta->{atomname} eq "HG22" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="3HG2" if ($ta->{atomname} eq "HG23" && $ta->{resname}=~/ILE/);
	} elsif ($translate =~ /ICM/) {
	  $ta->{resname}="HIS" if ($ta->{resname} =~/HSD/);
	  $ta->{resname}="HIS" if ($ta->{resname} =~/HID/);
	  $ta->{resname}="HIE" if ($ta->{resname} =~/HSE/);
	  $ta->{atomname}="1H" if ($ta->{atomname} eq "HT1");
	  $ta->{atomname}="2H" if ($ta->{atomname} eq "HT2");
	  $ta->{atomname}="3H" if ($ta->{atomname} eq "HT3");
	  $ta->{atomname}="HN" if ($ta->{atomname} eq "HN");
	  $ta->{atomname}="2H" if ($ta->{atomname} eq "HN1");
	  $ta->{atomname}="3H" if ($ta->{atomname} eq "HN2");
	  $ta->{atomname}="3HB" if ($ta->{atomname} eq "HB3" && $ta->{resname}=~/ALA/);
	  $ta->{atomname}="2HB" if ($ta->{atomname} eq "HB2" && $ta->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $ta->{atomname}="1HB" if ($ta->{atomname} eq "HB1" && $ta->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $ta->{atomname}="2HG" if ($ta->{atomname} eq "HG2" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $ta->{atomname}="1HG" if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $ta->{atomname}="2HA" if ($ta->{atomname} eq "HA2" && $ta->{resname}=~/GLY/);
	  $ta->{atomname}="1HA" if ($ta->{atomname} eq "HA1" && $ta->{resname}=~/GLY/);
	  $ta->{atomname}="3HE" if ($ta->{atomname} eq "HE3" && $ta->{resname}=~/MET/);
	  $ta->{atomname}="2HE" if ($ta->{atomname} eq "HE2" && $ta->{resname}=~/MET|LYS/);
	  $ta->{atomname}="1HE" if ($ta->{atomname} eq "HE1" && $ta->{resname}=~/MET|LYS/);
	  $ta->{atomname}="2HD" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/ARG|PRO/);
	  $ta->{atomname}="1HD" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/ARG|PRO/);
	  $ta->{atomname}="1HH1" if ($ta->{atomname} eq "HH11" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="2HH1" if ($ta->{atomname} eq "HH12" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="1HH2" if ($ta->{atomname} eq "HH21" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="2HH2" if ($ta->{atomname} eq "HH22" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="1HE2" if ($ta->{atomname} eq "HE21" && $ta->{resname}=~/GLN/);
	  $ta->{atomname}="2HE2" if ($ta->{atomname} eq "HE22" && $ta->{resname}=~/GLN/);
	  $ta->{atomname}="1HD1" if ($ta->{atomname} eq "HD11" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="2HD1" if ($ta->{atomname} eq "HD12" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="3HD1" if ($ta->{atomname} eq "HD13" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="1HD2" if ($ta->{atomname} eq "HD21" && $ta->{resname}=~/LEU|ASN/);
	  $ta->{atomname}="2HD2" if ($ta->{atomname} eq "HD22" && $ta->{resname}=~/LEU|ASN/);
	  $ta->{atomname}="3HD2" if ($ta->{atomname} eq "HD23" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="1HD" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/LYS|ARG|PRO/);
	  $ta->{atomname}="2HD" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/LYS|ARG|PRO/);
	  $ta->{atomname}="3HZ" if ($ta->{atomname} eq "HZ3" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="2HZ" if ($ta->{atomname} eq "HZ2" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="1HZ" if ($ta->{atomname} eq "HZ1" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="1HG1" if ($ta->{atomname} eq "HG11" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="2HG1" if ($ta->{atomname} eq "HG12" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="3HG1" if ($ta->{atomname} eq "HG13" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="1HG2" if ($ta->{atomname} eq "HG21" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="2HG2" if ($ta->{atomname} eq "HG22" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="3HG2" if ($ta->{atomname} eq "HG23" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="HG"  if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/SER|CYS/);
	  $ta->{atomname}="1HD1" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HD1" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="3HD1" if ($ta->{atomname} eq "HD3" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="1HG1" if ($ta->{atomname} eq "HG11" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HG1" if ($ta->{atomname} eq "HG12" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="1HG2" if ($ta->{atomname} eq "HG21" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HG2" if ($ta->{atomname} eq "HG22" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="3HG2" if ($ta->{atomname} eq "HG23" && $ta->{resname}=~/ILE/);
	} elsif (uc($translate) =~ /GENERIC/) {
	  $ta->{resname}=~s/HSD/HIS/;
	  $ta->{resname}=~s/HSE/HIS/;
	  $ta->{resname}=~s/HSP/HIS/;	  
	  $ta->{atomname}="CD1" 
	    if ($ta->{resname} eq "ILE" && $ta->{atomname} eq "CD");
	}

	if ( ( $translate !~ /NOH/ || $ta->{atomname}!~/^[0-9]*H.*/) &&
	     ( $allow_rosetta_water || $ta->{atomname}!~/[0-9]*W.*/ ) ) {
	  printf $fname "%s\n",&_pdbLine($ta,$translate=~/CHA/,$longaux2);
	}
	$prevres=$ta->{resnum};
      }
      printf $fname "TER\n";
    }
  }
  printf $fname "END\n";

  undef $fname;

  &GenUtil::compress($PDBname) if ($compress);

}

## method: readAmber(partopfile,coorfile)
## reads protein structure from Amber topology
## and coordinate files

sub readAmberPre6 {
  my $self=shift;
  my $partop=&GenUtil::getInputFile(shift);
  my $coorfile=&GenUtil::getInputFile(shift);

  my $mode="";

  my $natom;
  my $ntype;
  my $nres;

  my @aname;
  my @lbres;
  my @ipres;

  while (<$partop>) {
    chomp;
    if (/\%FLAG (.+)/) {
      ($mode=$1)=~s/ +//g;
    } elsif ($mode eq "POINTERS") {
      my @num=&GenUtil::readFormat($partop,8,30);
      $natom=$num[0];
      $ntype=$num[1];
      $nres=$num[11];
      $mode="";
    } elsif ($mode eq "ATOM_NAME") {
      @aname=&GenUtil::readFormat($partop,4,$natom);
      $mode="";
    } elsif ($mode eq "RESIDUE_LABEL") {
      @lbres=&GenUtil::readFormat($partop,4,$nres);
      $mode="";
    } elsif ($mode eq "RESIDUE_POINTER") {
      @ipres=&GenUtil::readFormat($partop,8,$nres);
      $mode="";
    }
  }
  undef $partop;

  <$coorfile>;
  my $coornatom=<$coorfile>;
  chomp $coornatom;
  
  die "number of atoms does not match: expected: $natom, found: $coornatom"
    if ($coornatom != $natom);

  my @x=&GenUtil::readFormat($coorfile,12,$coornatom*3);
  undef $coorfile;

  $ipres[$nres]=$natom+1;

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  
  my $chainrec=$self->_newChain("");

  my $ainx=1;
  for (my $ir=0; $ir<$nres; $ir++) {
    my $resrec={};
    ($resrec->{name}=$lbres[$ir])=~s/ +//g;
    $resrec->{num}=$ir+1;
    $resrec->{chain}="";
    $resrec->{start}=$ipres[$ir]-1;    
    $resrec->{end}=$ipres[$ir+1]-2;
    $resrec->{valid}=1;

    push (@{$chainrec->{res}},$resrec);
    
    for (my $ia=$ipres[$ir]; $ia<$ipres[$ir+1]; $ia++) {
      my $atomrec={};
      $atomrec->{atominx}=$ainx++;
      ($atomrec->{atomname}=$aname[$ia-1])=~s/ +//g;
      $atomrec->{resname}=$resrec->{name};
      $atomrec->{resnum}=$ir+1;
      $atomrec->{chain}="";
      $atomrec->{xcoor}=$x[($ia-1)*3];
      $atomrec->{ycoor}=$x[($ia-1)*3+1];
      $atomrec->{zcoor}=$x[($ia-1)*3+2];
      push (@{$chainrec->{atom}},$atomrec);
    }
  }

  $self->_coorCache();
}

## method: readAmber(partopfile,coorfile)
## reads protein structure from Amber topology
## and coordinate files
## Amber6

sub readAmber6 {
  my $self=shift;
  my $partop=&GenUtil::getInputFile(shift);
  my $coorfile=&GenUtil::getInputFile(shift);

  my $mode="";

  my $natom;
  my $ntype;
  my $nres;

  <$partop>;
  my @num=&GenUtil::readFormat($partop,6,30);
  $natom=$num[0];
  $ntype=$num[1];
  $nres=$num[11];
  my @aname=&GenUtil::readFormat($partop,4,$natom);
  my @junk=&GenUtil::readFormat($partop,16,$natom);
  @junk=&GenUtil::readFormat($partop,16,$natom);
  @junk=&GenUtil::readFormat($partop,6,$natom);
  @junk=&GenUtil::readFormat($partop,6,$natom);
  @junk=&GenUtil::readFormat($partop,6,$ntype*$ntype);
  my @lbres=&GenUtil::readFormat($partop,4,$nres);
  my @ipres=&GenUtil::readFormat($partop,6,$nres);
  undef $partop;

  <$coorfile>;
  my $coornatom=<$coorfile>;
  chomp $coornatom;
  
  die "number of atoms does not match: expected: $natom, found: $coornatom"
    if ($coornatom != $natom);

  my @x=&GenUtil::readFormat($coorfile,12,$coornatom*3);
  undef $coorfile;

  $ipres[$nres]=$natom+1;

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  
  my $chainrec=$self->_newChain("");

  my $ainx=1;
  for (my $ir=0; $ir<$nres; $ir++) {
    my $resrec={};
    ($resrec->{name}=$lbres[$ir])=~s/ +//g;
    $resrec->{num}=$ir+1;
    $resrec->{chain}="";
    $resrec->{start}=$ipres[$ir]-1;    
    $resrec->{end}=$ipres[$ir+1]-2;
    $resrec->{valid}=1;

    push (@{$chainrec->{res}},$resrec);
    
    for (my $ia=$ipres[$ir]; $ia<$ipres[$ir+1]; $ia++) {
      my $atomrec={};
      $atomrec->{atominx}=$ainx++;
      ($atomrec->{atomname}=$aname[$ia-1])=~s/ +//g;
      $atomrec->{resname}=$resrec->{name};
      $atomrec->{resnum}=$ir+1;
      $atomrec->{chain}="";
      $atomrec->{xcoor}=$x[($ia-1)*3];
      $atomrec->{ycoor}=$x[($ia-1)*3+1];
      $atomrec->{zcoor}=$x[($ia-1)*3+2];
      push (@{$chainrec->{atom}},$atomrec);
    }
  }

  $self->_coorCache();
}

## method: readAmber(partopfile,coorfile)
## reads protein structure from Amber topology
## and coordinate files
## Amber7

sub readAmber {
  my $self=shift;
  my $partop=&GenUtil::getInputFile(shift);
  my $coorfile=&GenUtil::getInputFile(shift);

  my $mode="";
  my $items=0;
  my $len=0;

  my @num=();
  my @aname=();
  my @lbres=();
  my @ipres=();

  my $format=0;

  while (<$partop>) {
    chomp;

    if (/\%FLAG ([A-Za-z0-9_]+)/) {
      $mode=$1;
    } elsif ( /\%FORMAT\(([0-9]+)[a-zA-Z]([0-9]+).*\)/) { 
      $items=$1;
      $len=$2;
     } elsif ($mode eq "POINTERS") {
      for (my $i=0; $i+$len<=length($_); $i+=$len) {
	push(@num,substr($_,$i,$len)+0);
      }      
    } elsif ($mode eq "ATOM_NAME") {
      for (my $i=0; $i+$len<=length($_); $i+=$len) {
	push(@aname,substr($_,$i,$len));
      }      
    } elsif ($mode eq "RESIDUE_LABEL") {
      for (my $i=0; $i+$len<=length($_); $i+=$len) {
	push(@lbres,substr($_,$i,$len));
      }      
    } elsif ($mode eq "RESIDUE_POINTER") {
      for (my $i=0; $i+$len<=length($_); $i+=$len) {
	push(@ipres,substr($_,$i,$len)+0);
      }      
    }
  }

  my $natom=$num[0];
  my $ntype=$num[1];
  my $nres=$num[11];

  undef $partop;

  <$coorfile>;
  my $coornatom=<$coorfile>;
  chomp $coornatom;
  
  die "number of atoms does not match: expected: $natom, found: $coornatom"
    if ($coornatom != $natom);

  my @x=&GenUtil::readFormat($coorfile,12,$coornatom*3);
  undef $coorfile;

  $ipres[$nres]=$natom+1;

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  
  my $chainrec=$self->_newChain("");

  my $ainx=1;
  for (my $ir=0; $ir<$nres; $ir++) {
    my $resrec={};
    ($resrec->{name}=$lbres[$ir])=~s/ +//g;
    $resrec->{num}=$ir+1;
    $resrec->{chain}="";
    $resrec->{start}=$ipres[$ir]-1;    
    $resrec->{end}=$ipres[$ir+1]-2;
    $resrec->{valid}=1;

    push (@{$chainrec->{res}},$resrec);
    
    for (my $ia=$ipres[$ir]; $ia<$ipres[$ir+1]; $ia++) {
      my $atomrec={};
      $atomrec->{atominx}=$ainx++;
      ($atomrec->{atomname}=$aname[$ia-1])=~s/ +//g;
      $atomrec->{resname}=$resrec->{name};
      $atomrec->{resnum}=$ir+1;
      $atomrec->{chain}="";
      $atomrec->{xcoor}=$x[($ia-1)*3];
      $atomrec->{ycoor}=$x[($ia-1)*3+1];
      $atomrec->{zcoor}=$x[($ia-1)*3+2];
      push (@{$chainrec->{atom}},$atomrec);
    }
  }

  $self->_coorCache();
}

## method: writeAmber(file)
## writes Amber coordinate file

sub writeAmber {
  my $self=shift;
  my $file=&GenUtil::getOutputFile(shift);
  
  my $c=$self->activeChains()->[0];

  print $file "\n";
  printf $file "%5d\n",$#{$c->{atom}}+1;
  my $i;
  my $a=$c->{atom};
  for ($i=0; $i<=$#{$a}; $i++) {
    printf $file "%12.7f%12.7f%12.7f",
    $a->[$i]->{xcoor}, $a->[$i]->{ycoor},$a->[$i]->{zcoor};
    printf $file "\n" if (($i%2)==1); 
  }
  printf $file "\n" if (($i%2)==1);
}

## method: slist = generateSegNames(segname)
## generates segment names for the current structure from
## the template given as the argument.

sub generateSegNames {
  my $self=shift;
  my $segname=shift;
  my $lastseg;

  $segname="PRO0" if (!defined $segname);

  die "empty molecule"
    if ($#{$self->activeChains()}<0);

  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{res}}>=0) {
      my @part;
      my $ipart=0;
      my $lastnum=-99999;
      for (my $ir=0; $ir<=$#{$c->{res}}; $ir++) {
	my $r=$c->{res}->[$ir];
	my $num=$r->{num};
        $ipart++ if ($num!=$lastnum+1); 
# && ($num!=$c->{res}->[0]->{num} || $lastnum!=-99999));
	$part[$ir]=$ipart;
	$lastnum=$num;
      }

      for (my $ir=0; $ir<=$#{$c->{res}}; $ir++) {
	my $r=$c->{res}->[$ir];
	if ($r->{name} eq "TIP3" || $r->{name} eq "HOH") {
	  $r->{seg}="TIP3";
	} elsif ($r->{name} =~ /GUA|ADE|URA|THY|CYT/) {
	  $r->{seg}="NA00";
	} else {
	  if ($#{$self->{chain}}>0) {
	    if ($ipart>1) {
	      $r->{seg}=sprintf("%1s%02d%1s",
		 substr($segname,0,1),$part[$ir],$c->{id});
	    } else {
	      $r->{seg}=sprintf("%3s%1s",
		 substr($segname,0,3),($c->{id} eq "" || $c->{id} eq " ")?"0":$c->{id});
	    }
	  } else {
	    if ($ipart>1) {
	      $r->{seg}=sprintf("%2s%02d",
                 substr($segname,0,2),$part[$ir]);
	    } else {
	      $r->{seg}=$segname;
	    }
	  }
	}

	for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	  $c->{atom}->[$i]->{seg}=$r->{seg};
	}
      }
    }
  }

  $self->{segmentlist}=undef;
}

## method: slist = generateSplitSegNames()
## generates segment names, splitting where chain
## breaks are found. Names go AA00-ZZ00. Returns a list
## of generated segment names.

sub generateSplitSegNames {
  my $self=shift;

  die "empty molecule" if ($#{$self->activeChains()}<0);

  my $segnum1=65;
  my $segnum2=65;
  my $segname=chr($segnum1).chr($segnum2)."00";
  my $Sqdist_cut=2.0*2.0;

  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{res}}>=0) {
	my $r=$c->{res}->[0];
	my $prevC;
	if ($r->{name} eq "TIP3" || $r->{name} eq "HOH") {
	  $r->{seg}="TIP3";
	} elsif ($r->{name} =~ /GUA|ADE|URA|THY|CYT/) {
	  $r->{seg}="NA00";
        } else {
	  $r->{seg}=$segname;
	  for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	      my $a=$c->{atom}->[$i];
	      $a->{seg}=$r->{seg};
	      $prevC=$a if ($a->{atomname} eq "C");
	  }
        }

        for (my $ir=1; $ir<=$#{$c->{res}}; $ir++) {
	   my $r=$c->{res}->[$ir];
	   if ($r->{name} eq "TIP3" || $r->{name} eq "HOH") {
	       $r->{seg}="TIP3";
	   } elsif ($r->{name} =~ /GUA|ADE|URA|THY|CYT/) {
	       $r->{seg}="NA00";
	   } else {
  	       # If the peptide bond length to the previous residue is
               # defined and is reasonable, keep the segid, otherwise
               # increment
	       my $currN;
	       my $currC;
	       for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
		   my $a=$c->{atom}->[$i];
		   $currN=$a if ($a->{atomname} eq "N");
		   $currC=$a if ($a->{atomname} eq "C");
	       }

	       if ((defined $prevC) && (defined $currN)) {
		   my $x1=$prevC->{xcoor};
		   my $x2=$currN->{xcoor};
		   my $y1=$prevC->{ycoor};
		   my $y2=$currN->{ycoor};
		   my $z1=$prevC->{zcoor};
		   my $z2=$currN->{zcoor};
 		   my $SqDist=($x1-$x2)*($x1-$x2);
		   $SqDist+=($y1-$y2)*($y1-$y2);
		   $SqDist+=($z1-$z2)*($z1-$z2);
		   if ($SqDist > $Sqdist_cut) {
		       $segnum2++;
		       if ($segnum2 == 91) {
			   $segnum1++;
			   $segnum2=65;
		       }
		       $segname=chr($segnum1).chr($segnum2)."00";
		   }
	       }
	       $prevC=$currC;
	       $r->{seg}=$segname;
	       for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
		   my $a=$c->{atom}->[$i];
		   $a->{seg}=$r->{seg};
		   $prevC=$a if ($a->{atomname} eq "C");
	       }

	   }
	}
    }

    # Increment the segid for the new chain
    $segnum2++;
    if ($segnum2 == 91) {
	$segnum1++;
	$segnum2=65;
    }
    $segname=chr($segnum1).chr($segnum2)."00";

  }

  $self->{segmentlist}=undef;
}

## method: slist = getSegNames()
## generates and returns segment name list for the current structure 

sub getSegNames {
  my $self=shift;

  my $lastseg;

  if (!defined $self->{segmentlist} || $#{$self->{segmentlist}}<0) {
    $self->{segmentlist}=();

    my $rec;

    foreach my $c ( @{$self->activeChains()} ) {
      if ($#{$c->{res}}>=0) {
	$lastseg=undef;
	foreach my $r ( @{$c->{res}} ) { 
	  if (!defined $lastseg || $r->{seg} ne $lastseg) {
	    $rec={};
	    $rec->{name}=$r->{seg};
	    $rec->{first}=(!defined $lastseg);
	    $rec->{last}=0;
	    $rec->{chain}=$c->{id};
	    $rec->{from}=$r->{num};
	    push(@{$self->{segmentlist}},$rec);
	    $lastseg=$r->{seg};
	  } 
	  $rec->{to}=$r->{num};
	}
	$rec->{last}=1;
      }
    }
  } 

  return $self->{segmentlist};
}


## method: copySegNames(mol)
## sets the segment names from another Molecule object

sub copySegNames {
  my $self=shift;
  my $mol=shift;

  foreach my $c ( @{$self->activeChains()} ) {
    my $mc=$mol->getChain($c->{id});
    $mc=$mol->getChain() if (!defined $mc && $#{$mol->{chain}}==0);
    if (defined $mc) {
      foreach my $r ( @{$c->{res}} ) {
	my $mr=$mol->getResidueInChain($r->{num},$mc);
	if (defined $mr) {
	  $r->{seg}=$mr->{seg};
	  for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	    $c->{atom}->[$i]->{seg}=$r->{seg};
	  }
	}
      }
    }
  }

  $self->{segmentlist}=undef;
}

## method: getChain([chainid])
## returns chain for given ID or first chain if
## no argument is given

sub getChain {
  my $self=shift;
  my $chainid=shift;

  return $self->{chain}->[0] 
    if (!defined $chainid || $chainid eq "");

  return $self->{chainlookup}->{uc $chainid};
}

## method: activeChains([chainid])
## returns a list of active chains

sub activeChains {
  my $self=shift;
  my $chainid=shift;

  my $list=();

  if (defined $chainid) {
    push(@{$list},$self->getChain($chainid));
    return $list;
  } elsif (defined $self->{defchain}) {
    push(@{$list},$self->getChain($self->{defchain}));
    return $list;
  } else {
    return $self->{chain};
  }
}

## method: selectChain(chainid)
## set a chain to be used as default

sub selectChain {
  my $self=shift;
  my $chainid=shift;

  $self->{defchain}=$chainid;
}

## method: removeChain()

sub removeChain {
  my $self=shift;
  my $cid=shift;

  for (my $i=0; $i<=$#{$self->{chain}}; $i++) {
      if ($self->{chain}->[$i]->{id} eq $cid) {
	  splice(@{$self->{chain}},$i,1);
      }
  }
  undef $self->{chainlookup}->{$cid};

  return;
}

## method: removeHetero()
## if there's a heteroatom chain, remove it

sub removeHetero {
  my $self=shift;

  $self->removeChain("+");

  return;
}

## method: numchains([protein_only])
## return the total number of (active) chains

sub numchains {
  my $self=shift;
  my $protein_only=shift;

  $protein_only=0 if (! defined $protein_only);

  return 0 if (! defined $self->activeChains());

  my $nc=0;
  foreach my $c ( @{$self->activeChains()} ) {
      $nc++ if ((! $protein_only) || ($c->{id} ne "+"));
  }

  return $nc;
}

## method: numres()
## return the total number of residues

sub numres {
  my $self=shift;

  my $nr=0;

  return $nr
      if (! defined $self->activeChains());

  foreach my $c ( @{$self->activeChains()} ) {
      $nr+=($#{$c->{res}})+1;
  }

  return $nr;
}

## method: lastResInx(chain)
## return the index of the last residue in the chain

sub lastResInx {
  my $self=shift;
  my $c=shift;

  return $#{$c->{res}};
}

## method: clearInfoLines()
## clear the info hash

sub clearInfoLines {
  my $self=shift;

  $self->{info}={};

  return;
}

## method: addInfoLine(tag, infoline)
## add a line to the info hash, with the key "tag"

sub addInfoLine {
  my $self=shift;
  my $tag=shift;
  my $infoline=shift;

  chomp($infoline);
  $self->{info}->{$tag}=() if (! defined $self->{info}->{$tag});
  push(@{$self->{info}->{$tag}},$infoline);

  return;
}

## method: computeSc()
## Call the "sc" module of CCP4
## Returns the sc value

sub computeSc {
  my $self=shift;

  my $sc_exec = &GenUtil::findExecutable("sc");
  die "cannot find the sc executable" if ( !defined $sc_exec );

  my $sc_script="tmp.".$$.".script";
  my $sc_pdb="tmp.".$$.".pdb";
  my $sc_output="tmp.".$$.".out";

  my $outscript=&GenUtil::getOutputFile($sc_script);
  my $chain1 = $self->{chain}->[0]->{id};
  my $chain2 = $self->{chain}->[1]->{id};
  printf $outscript "  MOLECULE 1\n  CHAIN %s\n  MOLECULE 2\n  CHAIN %s\n", $chain1, $chain2;
  undef $outscript;

  $self->writePDB($sc_pdb,translate=>"GENERIC_NOH",compress=>0,allow_rosetta_water=>0);

  system($sc_exec." XYZIN ".$sc_pdb." < ".$sc_script." > ".$sc_output);

  my $sc = 0.;
  my $sc_in=&GenUtil::getInputFile($sc_output);
  while (<$sc_in>) {
    chomp;
    if ( /Shape complementarity statistic Sc =    ([0-9]\.[0-9]+)/ ) {
      $sc = $1;
    }
  }
  undef $sc_in;

  &GenUtil::remove($sc_script);
  &GenUtil::remove($sc_pdb);
  &GenUtil::remove($sc_output);

  return $sc;
}

## method: countInterfaceContacts([distCut][,chainA][,chainB][,CACBonly][,SConly])
## Count the number of atomic contacts across an interface
## Operates on all pairs of chains if no chain IDs are specified
## Operates on all partners for given chain if one chain ID is specified
## Operates only on one pair if two chain IDs are specified
## Returns an array of hashes with {AAcontacts}, {cidA}, {cidB}

sub countInterfaceContacts {
  my $self=shift;
  my $distCut=shift;
  my $cidA=shift;
  my $cidB=shift;
  my $CACBonly=shift;
  my $SConly=shift;

  $distCut=5 if (! defined $distCut);
  my $SqDistCut=$distCut*$distCut;

  $CACBonly=0 if (! defined $CACBonly);
  $SConly=0 if (! defined $SConly);

  my @ac=@{$self->activeChains()};
  my $AAlist=();

  for (my $i=0; $i<=$#ac; $i++) {
    my $c1=$ac[$i];
    for (my $j=$i+1; $j<=$#ac; $j++) {
      my $c2=$ac[$j];
      if (($#{$c1->{res}}>=0) && ($#{$c2->{res}}>=0) &&
           ((! defined $cidA) || ($cidA eq $c1->{id}) ||
           ($cidA eq $c2->{id})) && ((! defined $cidB) ||
           ($cidB eq $c1->{id}) || ($cidB eq $c2->{id}))) {

	  my $AAcontacts=0;
	  foreach my $a1 ( @{$c1->{atom}} ) {
	    if (($a1->{aux1} >= 0.0) && (! $a1->{hyd}) &&
		((! $SConly) || (! $a1->{bb})) &&
		((! $CACBonly) ||
		 ($a1->{atomname} eq "CA") || ($a1->{atomname} eq "CB"))) {
	      my $x1=$a1->{xcoor};
	      my $y1=$a1->{ycoor};
	      my $z1=$a1->{zcoor};
	      foreach my $a2 ( @{$c2->{atom}} ) {
		  if (($a2->{aux1} >= 0.0) && (! $a2->{hyd}) &&
		      ((! $SConly) || (! $a2->{bb})) &&
		      ((! $CACBonly) ||
		       ($a2->{atomname} eq "CA") || ($a2->{atomname} eq "CB"))) {
		  my $x2=$a2->{xcoor};
		  my $y2=$a2->{ycoor};
		  my $z2=$a2->{zcoor};
		  my $SqDist=($x1-$x2)*($x1-$x2);
		  $SqDist+=($y1-$y2)*($y1-$y2);
		  $SqDist+=($z1-$z2)*($z1-$z2);
		  $AAcontacts++ if ($SqDist < $SqDistCut);

	        }
	      }
	    }
	  }

	  my $h={};
	  $h->{AAcontacts}=$AAcontacts;
	  $h->{cidA}=$c1->{id};
	  $h->{cidB}=$c2->{id};
	  push(@{$AAlist},$h);

      }
    }
  }

  return $AAlist;
}

## method: count_nonUBQcont_Contacts([distCut])
## Works like countInterfaceContacts, but doesn't count contacts in the first 76 residues of either chain
## Returns an integer

sub count_nonUBQcont_Contacts {
  my $self=shift;
  my $distCut=shift;

  $distCut=5 if (! defined $distCut);
  my $SqDistCut=$distCut*$distCut;

  my @ac=@{$self->activeChains()};
  my $AAcontacts=0;

  for (my $i=0; $i<=$#ac; $i++) {
    my $c1=$ac[$i];
    next if ($#{$c1->{res}} < 80);
    my $min_resnum1 = $c1->{res}->[78]->{num};
    for (my $j=$i+1; $j<=$#ac; $j++) {
      my $c2=$ac[$j];
      next if ($#{$c2->{res}} < 80);
      my $min_resnum2 = $c2->{res}->[78]->{num};

      foreach my $a1 ( @{$c1->{atom}} ) {
	next if ( $a1->{resnum} < $min_resnum1 );
	if (($a1->{aux1} >= 0.0) && (! $a1->{hyd})) {
	  my $x1=$a1->{xcoor};
	  my $y1=$a1->{ycoor};
	  my $z1=$a1->{zcoor};
	  foreach my $a2 ( @{$c2->{atom}} ) {
	    next if ( $a2->{resnum} < $min_resnum2 );
	    if (($a2->{aux1} >= 0.0) && (! $a2->{hyd})) {
	      my $x2=$a2->{xcoor};
	      my $y2=$a2->{ycoor};
	      my $z2=$a2->{zcoor};
	      my $SqDist=($x1-$x2)*($x1-$x2);
	      $SqDist+=($y1-$y2)*($y1-$y2);
	      $SqDist+=($z1-$z2)*($z1-$z2);
	      $AAcontacts++ if ($SqDist < $SqDistCut);
	    }
	  }
	}
      }

    }
  }

  return $AAcontacts;
}

## method: measure_UBQ_terminus_sqdist()
## Return the square of the distance between the termini of a UBQ design construct
## Note: measures both, returns the minimum

sub measure_UBQ_terminus_sqdist {
  my $self=shift;

  # Note: the UBQs have different lengths - specify the longest possible
  my $max_UBQ_length=76;

  my $min_sqdist=9999.;

  my $cA=$self->getChain("A");
  my $cB=$self->getChain("B");

  # Get the distance between the non-UBQ N-term of chain A and the C-term of chain B
  my $A_resinx = $max_UBQ_length;
  my $B_resinx = $#{$cB->{res}};
  my $A_atominx = $cA->{res}->[$A_resinx]->{start};
  my $B_atominx = $cB->{res}->[$B_resinx]->{start};
  my $aA = $cA->{atom}->[$A_atominx];
  my $aB = $cB->{atom}->[$B_atominx];
  my $x1=$aA->{xcoor};
  my $y1=$aA->{ycoor};
  my $z1=$aA->{zcoor};
  my $x2=$aB->{xcoor};
  my $y2=$aB->{ycoor};
  my $z2=$aB->{zcoor};
  my $SqDist=($x1-$x2)*($x1-$x2);
  $SqDist+=($y1-$y2)*($y1-$y2);
  $SqDist+=($z1-$z2)*($z1-$z2);
  $min_sqdist = $SqDist if ($SqDist < $min_sqdist);

  # Get the distance between the non-UBQ N-term of chain A and the C-term of chain B
  $A_resinx = $#{$cA->{res}};
  $B_resinx = $max_UBQ_length;
  $A_atominx = $cA->{res}->[$A_resinx]->{start};
  $B_atominx = $cB->{res}->[$B_resinx]->{start};
  $aA = $cA->{atom}->[$A_atominx];
  $aB = $cB->{atom}->[$B_atominx];
  $x1=$aA->{xcoor};
  $y1=$aA->{ycoor};
  $z1=$aA->{zcoor};
  $x2=$aB->{xcoor};
  $y2=$aB->{ycoor};
  $z2=$aB->{zcoor};
  $SqDist=($x1-$x2)*($x1-$x2);
  $SqDist+=($y1-$y2)*($y1-$y2);
  $SqDist+=($z1-$z2)*($z1-$z2);
  $min_sqdist = $SqDist if ($SqDist < $min_sqdist);

  return $min_sqdist;
}


## method: interfaceRes(res,distCut)
## Returns true if any heavy atom of this residue is
## within distCut of a heavy atom on a different (non-hetero) chain

sub interfaceRes {
  my $self=shift;
  my $res=shift;
  my $scOnly=shift;
  my $distCut=shift;

  $scOnly=0 if (! defined $scOnly);

  my $cid=$res->{chain};
  foreach my $c ( @{$self->activeChains()} ) {
      if ($c->{id} eq $cid) {
	  my $start=$res->{start};
	  my $end=$res->{end};
	  for (my $ai=$start; $ai<=$end; $ai++) {
	      my $atom=$c->{atom}->[$ai];
	      if (uc($res->{name}) ne "GLY") {
		  if ((! $atom->{hyd}) && (( ! $scOnly) || ( ! $atom->{bb}))) {
		      return 1 if ($self->interfaceAtom($atom,$distCut));
		  }
	      } else {
		  # Gly is an interface residue if the CA is in the interface
		  if (uc($atom->{atomname}) eq "CA") {
		      return 1 if ($self->interfaceAtom($atom,$distCut));
		  }
	      }
	  }
      }
  }

  return 0;
}

## method: interfaceAtom(atom,distCut)
## Returns true if this atom is within distCut
## of a heavy atom on a different (non-hetero) chain

sub interfaceAtom {
  my $self=shift;
  my $atom=shift;
  my $distCut=shift;

  $distCut=4 if (! defined $distCut);
  my $SqDistCut=$distCut*$distCut;

  my $cid=$atom->{chain};
  my $x1=$atom->{xcoor};
  my $y1=$atom->{ycoor};
  my $z1=$atom->{zcoor};

  foreach my $c ( @{$self->activeChains()} ) {
      if (($c->{id} ne $cid) && ($c->{id} ne "+")) {
	  my $na=$#{$c->{atom}};
	  my $a=$c->{atom};
	  for (my $ai=0; $ai<=$na; $ai++) {
	      my $a2=$a->[$ai];
	      if (! $a2->{hyd}) {
		  my $x2=$a2->{xcoor};
		  my $y2=$a2->{ycoor};
		  my $z2=$a2->{zcoor};
		  my $SqDist=($x1-$x2)*($x1-$x2);
		  $SqDist+=($y1-$y2)*($y1-$y2);
		  $SqDist+=($z1-$z2)*($z1-$z2);
		  return 1 if ($SqDist < $SqDistCut);
	      }
	  }
      }
  }

  return 0;
}


## method: checkCloseHetatm(res,distCut,scOnly)
## Returns true if there's a hetero atom close to this residue

sub checkCloseHetatm {
  my $self=shift;
  my $res=shift;
  my $distCut=shift;
  my $scOnly=shift;

  $scOnly=0 if (! defined $scOnly);
  $distCut=4.0 if (! defined $distCut);

  my $sqDistCut=$distCut*$distCut;

  my $rcid=$res->{chain};
  my $rstart=$res->{start};
  my $rend=$res->{end};
  foreach my $rc ( @{$self->activeChains()} ) {
      if (($rc->{id} eq $rcid) && ($#{$rc->{atom}}>=0)) {
	  foreach my $hc ( @{$self->activeChains()} ) {
	      if (($hc->{id} eq "+") && ($#{$hc->{atom}}>=0)) {
		  for (my $ai=$rstart; $ai<=$rend; $ai++) {
		      my $ratom=$rc->{atom}->[$ai];
		      my $rX=$ratom->{xcoor};
		      my $rY=$ratom->{ycoor};
		      my $rZ=$ratom->{zcoor};
		      if (((( ! $scOnly) || ( ! $ratom->{bb})) ||
			   ((uc($ratom->{resname}) eq "GLY") &&
			    ($ratom->{atomname} eq "CA"))) &&
			    (! $ratom->{hyd})) {
			  my $hend=$#{$hc->{xcoor}};
			  for (my $hai=1; $hai<=$hend; $hai++) {
			      my $hX=$hc->{xcoor}->[$hai];
			      my $hY=$hc->{ycoor}->[$hai];
			      my $hZ=$hc->{zcoor}->[$hai];
			      my $sqdist=($rX-$hX)*($rX-$hX);
			      $sqdist+=($rY-$hY)*($rY-$hY);
			      $sqdist+=($rZ-$hZ)*($rZ-$hZ);
			      return 1 if ($sqdist < $sqDistCut);
			  }
		      }
		  }
	      }
	  }
      }
  }

  return 0;
}


## method: setValidSequence(sequencestring)
## sets the <mark>valid</mark> flag to 1
## for residues that match the given sequence string
## this method recognizes a previous residue subselection

sub setValidSequence {
  my $self=shift;
  my $seqstr=shift;

  my $foundany=0;

  foreach my $c ( @{$self->activeChains()} ) {
    my $ir;
    for ($ir=0; $ir<=$#{$c->{res}}+1-length($seqstr); $ir++) {
      my $tir=0;
      while ($c->{res}->[$tir+$ir]->{valid} && $tir<length($seqstr) && 
	     &_cmpResName($c->{res}->[$tir+$ir]->{name},$Sequence::_seqlong{substr($seqstr,$tir,1)})) {
	$tir++;
      }
      if ($tir<length($seqstr)) {
	$c->{res}->[$ir]->{valid}=0;
      } else {
	$foundany=1;
	$ir+=length($seqstr)-1;
      }
    }
    for ( ;$ir<=$#{$c->{res}}; $ir++) {
      $c->{res}->[$ir]->{valid}=0;
    }
  } 
  return $foundany;
}

## method: setValidChain(chain[,exclude])
## sets the <mark>valid</mark> flag to 1 for
## all residues in the given chain 

sub setValidChain {
  my $self=shift;
  my $chain=shift;
  my $exclmode=shift;
  my $noreset=shift;

  return if (!defined $chain);

  $self->resetValidResidues($exclmode?1:0) unless (defined $noreset && $noreset);

  $exclmode=0 if (!defined $exclmode);

  my $c=$self->getChain($chain);
  if (defined $c) {
    foreach my $r ( @{$c->{res}} ) {
      $r->{valid}=($exclmode)?0:1;
    }    
  }
}

## method: setValidResidues(fraglist)
## sets the <mark>valid</mark> flag to 1
## for residues in the fragment list and to
## 0 for residues outside the list

sub setValidResidues {
  my $self=shift;
  my $fraglist=shift;
  my $exclmode=shift;
  my $noreset=shift;

  return if (!defined $fraglist);

  $exclmode=0 if (!defined $exclmode);

  $self->resetValidResidues($exclmode?1:0) unless (defined $noreset && $noreset);

  foreach my $l ( @{$fraglist} ) {
    my $c=$self->getChain($l->{chain});
    if (defined $c) {
      if (!defined $l->{from}) {
	$l->{from}=$self->firstResNum($c);
	$l->{to}=$self->lastResNum($c);
      } elsif (!defined $l->{to}) {
	$l->{to}=$l->{from};
      }
      for (my $i=$l->{from}; $i<=$l->{to}; $i++) {
	my $r=$self->getResidueInChain($i,$c);
	$r->{valid}=($exclmode?0:1) if (defined $r);
      }    
    }
  }
}

## method: setValidSegment(segname[,exclmode])
## sets the <mark>valid</mark> flag to 1
## for residues with the given segment name 

sub setValidSegment {
  my $self=shift;
  my $segname=shift;
  my $exclmode=shift;

  return if (!defined $segname);
  $exclmode=0 if (!defined $exclmode);

  $self->resetValidResidues($exclmode?1:0);

  foreach my $c ( @{$self->activeChains()} ) {
    foreach my $r ( @{$c->{res}} ) {
      $r->{valid}=($exclmode?0:1) if ($r->{seg} eq $segname);
    }
  }
}


## method: resetValidResidues([value])
## sets the <mark>valid</mark> flag of all
## residues to the given value (default: 1)

sub resetValidResidues {
  my $self=shift;
  my $value=shift;

  $value=1 if (!defined $value);
  
  foreach my $c ( @{$self->activeChains()} ) {
    foreach my $r ( @{$c->{res}} ) {
      $r->{valid}=$value;
    }
  }
}

## method: $list = listFromValid(force)
## returns a list of residues from the residues 
## previously set with <mark>setValidResidues</mark>.

sub listFromValid {
  my $self=shift;

  my $retlist=();
  foreach my $c ( @{$self->activeChains()} ) {
    my @arr;
    foreach my $r ( @{$c->{res}} ) {
      push(@arr,$r->{num}) if ($r->{valid});
    }
    if ($#arr>=0) {
      foreach my $tlist ( @{&GenUtil::fragListFromArray(\@arr,$c->{id})} ) {
	push(@{$retlist},$tlist);
      }
    } else {
      my $rec={};
      $rec->{from}=$self->firstResNum($c);
      $rec->{to}=$self->lastResNum($c);
      $rec->{chain}=$c->{id};
      push(@{$retlist},$rec);
    }
  }
  return $retlist;
}

## method: zapCoordinates([chain][,startres][,stopres])
## sets all coordinates in the current structure
## to zero 

sub zapCoordinates {
  my $self=shift;
  my $chainid=shift;
  my $start=shift;
  my $stop=shift;

  foreach my $c ( @{$self->activeChains($chainid)} ) {
    foreach my $a ( @{$c->{atom}} ) {
      if (((! defined $start) || ($a->{resnum} >= $start)) &&
	  ((! defined $stop) || ($a->{resnum} <= $stop))) {
	$a->{xcoor}=0.0;
	$a->{ycoor}=0.0;
	$a->{zcoor}=0.0;
      }
    }
    $self->_coorCache($c->{id});
  }
}

## method: enantiomer()
## generate the mirror image,
## by changing the sign on all of the z-coors

sub enantiomer {
  my $self=shift;

  foreach my $c ( @{$self->activeChains()} ) {
      foreach my $a ( @{$c->{atom}} ) {
	  $a->{zcoor}*=-1.0;
      }
  }
  $self->_coorCache();
}

## method: fillCoorFromPDB(file)
## updates only coordinates from a PDB file

sub fillCoorFromPDB {
  my $self=shift;
  my $pdb=Molecule->new(shift);

  foreach my $c ( @{$pdb->{chain}} ) {
    my $lc=$self->getChain($c->{id});

    if (defined $lc) {
      my %have;
      foreach my $a ( @{$lc->{atom}} ) {
	my $key=sprintf("%s%d:%s",
			$a->{resname},$a->{resnum},$a->{atomname});
	$have{$key}=$a;
      }

      foreach my $a ( @{$c->{atom}} ) {
	my $key=sprintf("%s%d:%s",
			$a->{resname},$a->{resnum},$a->{atomname});
	if (exists $have{$key}) {
	  $have{$key}->{xcoor}=$a->{xcoor};
	  $have{$key}->{ycoor}=$a->{ycoor};
	  $have{$key}->{zcoor}=$a->{zcoor};
	}
      }
      $self->_coorCache($c->{id});
    }
  }

  undef $pdb;
}

## method: merge(refmol)
## merges the current structure with the structure
## in <mark>refmol</mark> on a per residue basis

sub merge {
  my $self=shift;
  my $ref=shift;

  my $have={};
  my $cid=();

  foreach my $c ( @{$ref->activeChains()} ) {
    my $sc=$self->getChain($c->{id});
    if (!defined $sc) {
      my $crec=$self->_newChain($c->{id});
      my $ainx=1;
      foreach my $r ( @{$c->{res}} ) {
	if ($r->{valid}) {
	  my $rrec={};
	  %{$rrec}=%{$r};
	  push(@{$crec->{res}},$rrec);
	  $rrec->{start}=$#{$crec->{atom}}+1;
	  for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	    my $arec={};
	    %{$arec}=%{$c->{atom}->[$ia]};
	    $arec->{atominx}=$ainx++;
	    push(@{$crec->{atom}},$arec);
	  }
	  $rrec->{end}=$#{$crec->{atom}};
	}
      }
    } else {
      my %have;
      foreach my $r ( @{$sc->{res}} ) {
	my $trec={res=>$r, atom=>$sc->{atom}};
	$have{$r->{num}}=$trec;
      }
      foreach my $r ( @{$c->{res}} ) {
	if ($r->{valid}) {
	  my $trec={res=>$r, atom=>$c->{atom}};
	  $have{$r->{num}}=$trec;
	}
      }
      
      my $newatom=();
      my $newres=();
      my $ainx=1;
      foreach my $k ( sort { $a<=>$b } keys %have ) {
	my $r=$have{$k}->{res};
	my $a=$have{$k}->{atom};
	my $rrec={};
	%{$rrec}=%{$r};
	push(@{$newres},$rrec);
	$rrec->{start}=$#{$newatom}+1;
	for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	  my $arec={};
	  %{$arec}=%{$a->[$ia]};
	  $arec->{atominx}=$ainx++;
	  push(@{$newatom},$arec);
	}
	$rrec->{end}=$#{$newatom};
      }
      $sc->{atom}=$newatom;
      $sc->{res}=$newres;
    }
  }

  # when we merge, keep the OLD header
  # note: commented code below concatenates headers,
  # but this is dangerous if used within a loop!
#  foreach my $k ( keys %{$ref->{info}} ) {
#      my $aref=$ref->{info}->{$k};
#      $k=uc($k);
#      $self->{info}->{$k}=() if (! defined $self->{info}->{$k});
#      my $a=$self->{info}->{$k};
#      foreach my $line ( @{$aref} ) {
#	  push(@{$a},$line);
#      }
#  }

  $self->{segmentlist}=undef;
  $self->_coorCache();
}

## method: $mol = clone([valid])
## generates a copy of the current molecule
## and returns it. If the valid flag is set, only 
## valid residues are copied

sub clone {
  my $self=shift;
  my $valid=shift;

  my $n=$self->new();

  foreach my $s ( @{$self->{ssbond}} ) {
    my $srec={};
    %{$srec}=%{$s};
    push (@{$n->{ssbond}},$srec);
  }

  foreach my $c ( @{$self->activeChains()} ) {
    my $nc=undef;

    foreach my $r (@{$c->{res}}) {
      if ($r->{valid} || !defined $valid || !$valid) {
	$nc=$n->_newChain($c->{id}) if (!defined $nc);

	my $rrec={};
	%{$rrec}=%{$r};
	$rrec->{start}=$#{$nc->{atom}}+1;
	push(@{$nc->{res}},$rrec);
	for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	  my $arec={};
	  %{$arec}=%{$c->{atom}->[$ia]};
	  push(@{$nc->{atom}},$arec);
	}
	$rrec->{end}=$#{$nc->{atom}};
      }
    }
  }

  # Clone the info hash
  $n->{info}={};
  foreach my $k ( keys %{$self->{info}} ) {
      if ($#{$self->{info}->{$k}} >= 0) {
	  $n->{info}->{$k}=();
	  foreach my $infoline ( @{$self->{info}->{$k}} ) {
	      push(@{$n->{info}->{$k}},$infoline);
	  }
      }
  }

  $n->_coorCache();

  return $n;
}

## method: translate(mode)
## translates the current structure according to the mode argument.
## As with <mark>writePDB</mark> the following modes are recognized:
## <mark>CHARMM19</mark>, <mark>CHARMM22</mark>,
## <mark>AMBER</mark>, and <mark>GENERIC</mark>.

sub translate {
  my $self=shift;
  my $translate=shift;
  my $reqHis=shift;

  $translate="GENERIC" if ((!defined $translate) || ($translate eq ""));
  $reqHis="" if (!defined $reqHis);

  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{atom}}>=0) {
      my $cterm=0;
      my $lastres=$c->{res}->[$#{$c->{res}}];
      for (my $i=$lastres->{start}; $i<=$lastres->{end}; $i++) {
	if ($c->{atom}->[$i]->{atomname} eq "OXT") {
	  $cterm=1;
	}
      }

      foreach my $r ( @{$c->{res}} ) {
	if ($translate=~/CHARMM/) {
	  $r->{name}=~s/HOH/TIP3/;
	  if ($r->{name}=~/HSD|HSE|HIS|HSP/) {
	    my $hisName=$r->{name};
	    if ($translate =~ /CHARMM19/) {
	      if ($hisName eq "HSD") {
		$hisName="HIS";
	      } elsif ($hisName eq "HSE") {
		$hisName="HSD";
	      } elsif ($hisName eq "HSP") {
		$hisName="HSC";
	      }
	    } elsif (($hisName eq "HIS") && ($translate =~ /CHARMM22/)) {
	      $hisName="HSD";
	    }
	    $r->{name}=$hisName;
            for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	      $c->{atom}->[$i]->{resname}=$hisName;
	    }
	  }

        } elsif ($translate =~ /AMBER/) {
	  $r->{name}=~s/HSD/HID/;
	  $r->{name}=~s/HSE/HIE/;
	} elsif ($translate =~ /GENERIC/) {
	  $r->{name}=~s/HSD/HIS/;
	  $r->{name}=~s/HSE/HIS/;
	}
      }

      foreach my $a ( @{$c->{atom}} ) {
	if ($translate=~/CHARMM/) {
	  $a->{atomname}="CD" 
	    if ($a->{resname} eq "ILE" && $a->{atomname} eq "CD1");
	  $a->{atomname}="OT1"
	    if ($a->{atomname} eq "O" && $a->{resnum} eq $lastres->{num} && $cterm &&
		$a->{resname} =~/ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HIS|HSP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX/);
	  $a->{atomname}="OT2"
	    if ($a->{atomname} eq "OXT" && $a->{resnum} eq $lastres->{num});
	  $a->{resname}=~s/HOH/TIP3/;
	} 

	if ($translate =~ /AMBER/) {
	  $a->{resname}=~s/HSD/HID/;
	  $a->{resname}=~s/HSE/HIE/;
	  $a->{atomname}="CD1" 
	    if ($a->{resname} eq "ILE" && $a->{atomname} eq "CD");

	  if ($translate !~/CHAMBER/) {
	    $a->{atomname}="H1" if ($a->{atomname} eq "HT1");
	    $a->{atomname}="H2" if ($a->{atomname} eq "HT2");
	    $a->{atomname}="H3" if ($a->{atomname} eq "HT3");
	    $a->{atomname}="H" if ($a->{atomname} eq "HN");
	    $a->{atomname}="HB3" if ($a->{atomname} eq "HB2" && $a->{resname}=~/MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HID|HIE/);
	    $a->{atomname}="HB2" if ($a->{atomname} eq "HB1" && $a->{resname}=~/MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HID|HIE/);
	    $a->{atomname}="HG3" if ($a->{atomname} eq "HG2" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	    $a->{atomname}="HG2" if ($a->{atomname} eq "HG1" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	    $a->{atomname}=(($a->{resnum} eq $c->{atom}->[0]->{resnum}) || ($a->{resnum} eq $lastres->{num} && $cterm))
	      ?"HSG":"HG" if ($a->{atomname} eq "HG1" && $a->{resname}=~/CYS/);
	    $a->{atomname}="HD3" if ($a->{atomname} eq "HD2" && $a->{resname}=~/LYS|ARG|PRO/);
	    $a->{atomname}="HD2" if ($a->{atomname} eq "HD1" && $a->{resname}=~/LYS|ARG|PRO/);
	    $a->{atomname}="HE3" if ($a->{atomname} eq "HE2" && $a->{resname}=~/LYS/);
	    $a->{atomname}="HE2" if ($a->{atomname} eq "HE1" && $a->{resname}=~/LYS/);
	    $a->{atomname}="HA3" if ($a->{atomname} eq "HA2" && $a->{resname}=~/GLY/);
	    $a->{atomname}="HA2" if ($a->{atomname} eq "HA1" && $a->{resname}=~/GLY/);
	    $a->{atomname}="HG"  if ($a->{atomname} eq "HG1" && $a->{resname}=~/SER/);
	    $a->{atomname}="HD11" if ($a->{atomname} eq "HD1" && $a->{resname}=~/ILE/);
	    $a->{atomname}="HD12" if ($a->{atomname} eq "HD2" && $a->{resname}=~/ILE/);
	    $a->{atomname}="HD13" if ($a->{atomname} eq "HD3" && $a->{resname}=~/ILE/);
	    $a->{atomname}="HG13" if ($a->{atomname} eq "HG12" && $a->{resname}=~/ILE/);
	    $a->{atomname}="HG12" if ($a->{atomname} eq "HG11" && $a->{resname}=~/ILE/);
	    $a->{atomname}="H2" if ($a->{atomname} eq "HN1" && $a->{resname}=~/PRO/);
	    $a->{atomname}="H3" if ($a->{atomname} eq "HN2" && $a->{resname}=~/PRO/);
	  }
	} elsif ($translate =~ /IMPACT/) {
	  $a->{resname}="HID" if ($a->{resname} =~/HSD/);
	  $a->{resname}="HIE" if ($a->{resname} =~/HSE/);
	  $a->{atomname}="1H" if ($a->{atomname} eq "HT1");
	  $a->{atomname}="2H" if ($a->{atomname} eq "HT2");
	  $a->{atomname}="3H" if ($a->{atomname} eq "HT3");
	  $a->{atomname}="H" if ($a->{atomname} eq "HN");
	  $a->{atomname}="2H" if ($a->{atomname} eq "HN1");
	  $a->{atomname}="3H" if ($a->{atomname} eq "HN2");
	  $a->{atomname}="3HB" if ($a->{atomname} eq "HB3" && $a->{resname}=~/ALA/);
	  $a->{atomname}="2HB" if ($a->{atomname} eq "HB2" && $a->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $a->{atomname}="1HB" if ($a->{atomname} eq "HB1" && $a->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $a->{atomname}="2HG" if ($a->{atomname} eq "HG2" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $a->{atomname}="1HG" if ($a->{atomname} eq "HG1" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $a->{atomname}="2HA" if ($a->{atomname} eq "HA2" && $a->{resname}=~/GLY/);
	  $a->{atomname}="1HA" if ($a->{atomname} eq "HA1" && $a->{resname}=~/GLY/);
	  $a->{atomname}="3HE" if ($a->{atomname} eq "HE3" && $a->{resname}=~/MET/);
	  $a->{atomname}="2HE" if ($a->{atomname} eq "HE2" && $a->{resname}=~/MET|LYS/);
	  $a->{atomname}="1HE" if ($a->{atomname} eq "HE1" && $a->{resname}=~/MET|LYS/);
	  $a->{atomname}="2HD" if ($a->{atomname} eq "HD2" && $a->{resname}=~/ARG|PRO/);
	  $a->{atomname}="1HD" if ($a->{atomname} eq "HD1" && $a->{resname}=~/ARG|PRO/);
	  $a->{atomname}="1HH1" if ($a->{atomname} eq "HH11" && $a->{resname}=~/ARG/);
	  $a->{atomname}="2HH1" if ($a->{atomname} eq "HH12" && $a->{resname}=~/ARG/);
	  $a->{atomname}="1HH2" if ($a->{atomname} eq "HH21" && $a->{resname}=~/ARG/);
	  $a->{atomname}="2HH2" if ($a->{atomname} eq "HH22" && $a->{resname}=~/ARG/);
	  $a->{atomname}="1HE2" if ($a->{atomname} eq "HE21" && $a->{resname}=~/GLN/);
	  $a->{atomname}="2HE2" if ($a->{atomname} eq "HE22" && $a->{resname}=~/GLN/);
	  $a->{atomname}="1HD1" if ($a->{atomname} eq "HD11" && $a->{resname}=~/LEU/);
	  $a->{atomname}="2HD1" if ($a->{atomname} eq "HD12" && $a->{resname}=~/LEU/);
	  $a->{atomname}="3HD1" if ($a->{atomname} eq "HD13" && $a->{resname}=~/LEU/);
	  $a->{atomname}="1HD2" if ($a->{atomname} eq "HD21" && $a->{resname}=~/LEU|ASN/);
	  $a->{atomname}="2HD2" if ($a->{atomname} eq "HD22" && $a->{resname}=~/LEU|ASN/);
	  $a->{atomname}="3HD2" if ($a->{atomname} eq "HD23" && $a->{resname}=~/LEU/);
	  $a->{atomname}="1HD" if ($a->{atomname} eq "HD1" && $a->{resname}=~/LYS|ARG|PRO/);
	  $a->{atomname}="2HD" if ($a->{atomname} eq "HD2" && $a->{resname}=~/LYS|ARG|PRO/);
	  $a->{atomname}="3HZ" if ($a->{atomname} eq "HZ3" && $a->{resname}=~/LYS/);
	  $a->{atomname}="2HZ" if ($a->{atomname} eq "HZ2" && $a->{resname}=~/LYS/);
	  $a->{atomname}="1HZ" if ($a->{atomname} eq "HZ1" && $a->{resname}=~/LYS/);
	  $a->{atomname}="1HG1" if ($a->{atomname} eq "HG11" && $a->{resname}=~/VAL/);
	  $a->{atomname}="2HG1" if ($a->{atomname} eq "HG12" && $a->{resname}=~/VAL/);
	  $a->{atomname}="3HG1" if ($a->{atomname} eq "HG13" && $a->{resname}=~/VAL/);
	  $a->{atomname}="1HG2" if ($a->{atomname} eq "HG21" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="2HG2" if ($a->{atomname} eq "HG22" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="3HG2" if ($a->{atomname} eq "HG23" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="HG"  if ($a->{atomname} eq "HG1" && $a->{resname}=~/SER|CYS/);
	  $a->{atomname}="1HD1" if ($a->{atomname} eq "HD1" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HD1" if ($a->{atomname} eq "HD2" && $a->{resname}=~/ILE/);
	  $a->{atomname}="3HD1" if ($a->{atomname} eq "HD3" && $a->{resname}=~/ILE/);
	  $a->{atomname}="1HG1" if ($a->{atomname} eq "HG11" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HG1" if ($a->{atomname} eq "HG12" && $a->{resname}=~/ILE/);
	  $a->{atomname}="1HG2" if ($a->{atomname} eq "HG21" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HG2" if ($a->{atomname} eq "HG22" && $a->{resname}=~/ILE/);
	  $a->{atomname}="3HG2" if ($a->{atomname} eq "HG23" && $a->{resname}=~/ILE/);
	} elsif ($translate =~ /ICM/) {
	  $a->{resname}="HIS" if ($a->{resname} =~/HSD/);
	  $a->{resname}="HIS" if ($a->{resname} =~/HID/);
	  $a->{resname}="HIE" if ($a->{resname} =~/HSE/);
	  $a->{atomname}="1H" if ($a->{atomname} eq "HT1");
	  $a->{atomname}="2H" if ($a->{atomname} eq "HT2");
	  $a->{atomname}="3H" if ($a->{atomname} eq "HT3");
	  $a->{atomname}="HN" if ($a->{atomname} eq "HN");
	  $a->{atomname}="2H" if ($a->{atomname} eq "HN1");
	  $a->{atomname}="3H" if ($a->{atomname} eq "HN2");
	  $a->{atomname}="3HB" if ($a->{atomname} eq "HB3" && $a->{resname}=~/ALA/);
	  $a->{atomname}="2HB" if ($a->{atomname} eq "HB2" && $a->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $a->{atomname}="1HB" if ($a->{atomname} eq "HB1" && $a->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $a->{atomname}="2HG" if ($a->{atomname} eq "HG2" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $a->{atomname}="1HG" if ($a->{atomname} eq "HG1" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $a->{atomname}="2HA" if ($a->{atomname} eq "HA2" && $a->{resname}=~/GLY/);
	  $a->{atomname}="1HA" if ($a->{atomname} eq "HA1" && $a->{resname}=~/GLY/);
	  $a->{atomname}="3HE" if ($a->{atomname} eq "HE3" && $a->{resname}=~/MET/);
	  $a->{atomname}="2HE" if ($a->{atomname} eq "HE2" && $a->{resname}=~/MET|LYS/);
	  $a->{atomname}="1HE" if ($a->{atomname} eq "HE1" && $a->{resname}=~/MET|LYS/);
	  $a->{atomname}="2HD" if ($a->{atomname} eq "HD2" && $a->{resname}=~/ARG|PRO/);
	  $a->{atomname}="1HD" if ($a->{atomname} eq "HD1" && $a->{resname}=~/ARG|PRO/);
	  $a->{atomname}="1HH1" if ($a->{atomname} eq "HH11" && $a->{resname}=~/ARG/);
	  $a->{atomname}="2HH1" if ($a->{atomname} eq "HH12" && $a->{resname}=~/ARG/);
	  $a->{atomname}="1HH2" if ($a->{atomname} eq "HH21" && $a->{resname}=~/ARG/);
	  $a->{atomname}="2HH2" if ($a->{atomname} eq "HH22" && $a->{resname}=~/ARG/);
	  $a->{atomname}="1HE2" if ($a->{atomname} eq "HE21" && $a->{resname}=~/GLN/);
	  $a->{atomname}="2HE2" if ($a->{atomname} eq "HE22" && $a->{resname}=~/GLN/);
	  $a->{atomname}="1HD1" if ($a->{atomname} eq "HD11" && $a->{resname}=~/LEU/);
	  $a->{atomname}="2HD1" if ($a->{atomname} eq "HD12" && $a->{resname}=~/LEU/);
	  $a->{atomname}="3HD1" if ($a->{atomname} eq "HD13" && $a->{resname}=~/LEU/);
	  $a->{atomname}="1HD2" if ($a->{atomname} eq "HD21" && $a->{resname}=~/LEU|ASN/);
	  $a->{atomname}="2HD2" if ($a->{atomname} eq "HD22" && $a->{resname}=~/LEU|ASN/);
	  $a->{atomname}="3HD2" if ($a->{atomname} eq "HD23" && $a->{resname}=~/LEU/);
	  $a->{atomname}="1HD" if ($a->{atomname} eq "HD1" && $a->{resname}=~/LYS|ARG|PRO/);
	  $a->{atomname}="2HD" if ($a->{atomname} eq "HD2" && $a->{resname}=~/LYS|ARG|PRO/);
	  $a->{atomname}="3HZ" if ($a->{atomname} eq "HZ3" && $a->{resname}=~/LYS/);
	  $a->{atomname}="2HZ" if ($a->{atomname} eq "HZ2" && $a->{resname}=~/LYS/);
	  $a->{atomname}="1HZ" if ($a->{atomname} eq "HZ1" && $a->{resname}=~/LYS/);
	  $a->{atomname}="1HG1" if ($a->{atomname} eq "HG11" && $a->{resname}=~/VAL/);
	  $a->{atomname}="2HG1" if ($a->{atomname} eq "HG12" && $a->{resname}=~/VAL/);
	  $a->{atomname}="3HG1" if ($a->{atomname} eq "HG13" && $a->{resname}=~/VAL/);
	  $a->{atomname}="1HG2" if ($a->{atomname} eq "HG21" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="2HG2" if ($a->{atomname} eq "HG22" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="3HG2" if ($a->{atomname} eq "HG23" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="HG"  if ($a->{atomname} eq "HG1" && $a->{resname}=~/SER|CYS/);
	  $a->{atomname}="1HD1" if ($a->{atomname} eq "HD1" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HD1" if ($a->{atomname} eq "HD2" && $a->{resname}=~/ILE/);
	  $a->{atomname}="3HD1" if ($a->{atomname} eq "HD3" && $a->{resname}=~/ILE/);
	  $a->{atomname}="1HG1" if ($a->{atomname} eq "HG11" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HG1" if ($a->{atomname} eq "HG12" && $a->{resname}=~/ILE/);
	  $a->{atomname}="1HG2" if ($a->{atomname} eq "HG21" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HG2" if ($a->{atomname} eq "HG22" && $a->{resname}=~/ILE/);
	  $a->{atomname}="3HG2" if ($a->{atomname} eq "HG23" && $a->{resname}=~/ILE/);
	} elsif ($translate =~ /GENERIC/) {
	  $a->{resname}=~s/HSD/HIS/;
	  $a->{resname}=~s/HSE/HIS/;
	  $a->{resname}=~s/HSP/HIS/;
	  $a->{atomname}="CD1" 
	    if ($a->{resname} eq "ILE" && $a->{atomname} eq "CD");
	  $a->{atomname}="O"
	    if ($a->{atomname} eq "OT1");
	  $a->{atomname}="OXT"
	    if ($a->{atomname} eq "OT2");
	}
      }
    }
  }
}

## method: $map = renumber([map|startindex][,chain])
## renumbers the residues either according
## to a map given as an argument or continuously
## beginning from a start index. The method returns
## the translation map which may be applied
## to other structures or used to revert the
## residue numbering to its original form using
## <mark>numberReset</mark>

sub renumber {
  my $self=shift;

  my $start=1;
  my $map;
  my $t=shift;
  if (defined $t) {
    if (ref $t) {
      $map=$t;
    } else {
      $start=$t;
    }
  }

  my $c=$self->activeChains(shift)->[0];

  my $atom=$c->{atom};
  my $res=$c->{res};

  my $nres=$start-1;
  my $lastnum=-999;

  my $havemap=(defined $map);

  for (my $i=0; $i<=$#{$atom}; $i++) {  
    $atom->[$i]->{atominx}=$i+1;

    if ($atom->[$i]->{resnum} != $lastnum) {
      $lastnum=$atom->[$i]->{resnum};
      $nres++;
    }

    my $newnum=$havemap ? $map->{$atom->[$i]->{resnum}} : $nres;
    $atom->[$i]->{resnum}=$newnum;
  }

  my @nn;
  for (my $i=0; $i<=$#{$res}; $i++) {
    my $newnum;
    if ($havemap) {
      $newnum=$map->{$res->[$i]->{num}};
    } else {
      $map->{$res->[$i]->{num}}=$i+$start;
      $newnum=$i+$start;
    }
    $nn[$i]=$newnum;
  }

  foreach my $s ( @{$self->{ssbond}} ) {
  ADD1:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum1}==$res->[$i]->{num} &&
	  $s->{chain1}==$res->[$i]->{chain}) {
	$s->{resnum1}=$nn[$i];
	last ADD1;
      }
    }
  ADD2:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum2}==$res->[$i]->{num} &&
	  $s->{chain2}==$res->[$i]->{chain}) {
	$s->{resnum2}=$nn[$i];
	last ADD2;
      }
    }
  }

  for (my $i=0; $i<=$#{$res}; $i++) {
    $res->[$i]->{num}=$nn[$i];
  }

  $c->{resinx}=undef;

  return $map;
}

## method: renumberAcrossChains()
## renumber sequentially across the chains in order

sub renumberAcrossChains {
  my $self = shift;

  return if ( !defined $self->activeChains() );

  my $nextRes = 1;
  foreach my $c ( @{ $self->activeChains() } ) {
    $self->renumber( $nextRes, $c->{id} );
    $nextRes = $c->{res}->[ $#{ $c->{res} } ]->{num} + 1;
  }

  return;
}

## method: numberReset(map[,chain])
## applies the reverse translation from the one
## given in the map argument

sub numberReset {
  my $self=shift;
  my $map=shift;

  my $c=$self->activeChains(shift)->[0];

  my $reverseMap={};
  
  foreach my $m ( keys %{$map} ) {
    $reverseMap->{$map->{$m}}=$m;
  }

  my $atom=$c->{atom};
  my $res=$c->{res};

  for (my $i=0; $i<=$#{$atom}; $i++) {  
    my $newnum=$reverseMap->{$atom->[$i]->{resnum}};
    $atom->[$i]->{resnum}=$newnum;
  }

  foreach my $s ( @{$self->{ssbond}} ) {
  ADD1:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum1}==$res->[$i]->{num} &&
	  $s->{chain1}==$res->[$i]->{chain}) {
	$s->{resnum1}=$reverseMap->{$res->[$i]->{num}};
	last ADD1;
      }
    }
  ADD2:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum2}==$res->[$i]->{num} &&
	  $s->{chain2}==$res->[$i]->{chain}) {
	$s->{resnum2}=$reverseMap->{$res->[$i]->{num}};
	last ADD2;
      }
    }
  }

  for (my $i=0; $i<=$#{$res}; $i++) {
    my $newnum=$reverseMap->{$res->[$i]->{num}};
    $res->[$i]->{num}=$newnum;
  }

  $c->{resinx}=undef;
}

## method: shiftResNumber(add[,chain])
## shifts all residue numbers by a constant. A chain ID
## may be given to shift residues in that chain

sub shiftResNumber {
  my $self=shift;
  my $add=shift;
  my $cid=shift;

  my $c=$self->activeChains($cid)->[0];

  my $atom=$c->{atom};
  my $res=$c->{res};

  for (my $i=0; $i<=$#{$atom}; $i++) {  
    $atom->[$i]->{resnum}+=$add
  }

  foreach my $s ( @{$self->{ssbond}} ) {
  ADD1:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum1}==$res->[$i]->{num} &&
	  $s->{chain1}==$res->[$i]->{chain}) {
	$s->{resnum1}+=$add;
	last ADD1;
      }
    }
  ADD2:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum2}==$res->[$i]->{num} &&
	  $s->{chain2}==$res->[$i]->{chain}) {
	$s->{resnum2}+=$add;
	last ADD2;
      }
    }
  }

  for (my $i=0; $i<=$#{$res}; $i++) {
    $res->[$i]->{num}+=$add
  }
  
  $c->{resinx}=undef;
}

## method: ($cx,$cy,$cz) = centerOfMass()
## calculates the center of mass for the current structure

sub centerOfMass {
  my $self=shift;
  my $chain=shift;
  
  my $n;
  my ($cx,$cy,$cz)=(0.0,0.0,0.0);

  foreach my $c ( @{$self->activeChains($chain)} ) {
    my $atom=$c->{atom};
    $n+=$#{$atom}+1;
    
    for (my $i=0; $i<$n; $i++) {
      $cx+=$atom->[$i]->{xcoor};
      $cy+=$atom->[$i]->{ycoor};
      $cz+=$atom->[$i]->{zcoor};
    }
    
    $cx/=$n;
    $cy/=$n;
    $cz/=$n;
  }
  
  return ($cx,$cy,$cz);
}  

## method: center()
## centers the current structure by subtracting the center 
## of mass from all coordinates

sub center {
  my $self=shift;
  my $chain=shift;
  
  my ($cx,$cy,$cz)=$self->centerOfMass($chain);

  foreach my $c ( @{$self->activeChains($chain)} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      $atom->[$i]->{xcoor}-=$cx;
      $atom->[$i]->{ycoor}-=$cy;
      $atom->[$i]->{zcoor}-=$cz;
    }
  }
}

## method: move(dx,dy,dz,chain)
## shifts the current structure by the given distances
## in x/y/z direction

sub move {
  my $self=shift;
  my $dx=shift;
  my $dy=shift;
  my $dz=shift;
  my $chain=shift;

  foreach my $c ( @{$self->activeChains($chain)} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      $atom->[$i]->{xcoor}+=$dx;
      $atom->[$i]->{ycoor}+=$dy;
      $atom->[$i]->{zcoor}+=$dz;
    }
  }

  $self->_coorCache();
  return;
}

## method: matrixOperation(A,B,C,D,E,F,G,H,I [,chain])
## apply the matrix operator:
##  A  B  C
##  D  E  F
##  G  H  I

sub matrixOperation {
  my $self=shift;
  my $A=shift;
  my $B=shift;
  my $C=shift;
  my $D=shift;
  my $E=shift;
  my $F=shift;
  my $G=shift;
  my $H=shift;
  my $I=shift;
  my $chain=shift;

  foreach my $c ( @{$self->activeChains($chain)} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      my $x=$atom->[$i]->{xcoor};
      my $y=$atom->[$i]->{ycoor};
      my $z=$atom->[$i]->{zcoor};

      $atom->[$i]->{xcoor}=($A*$x)+($B*$y)+($C*$z);
      $atom->[$i]->{ycoor}=($D*$x)+($E*$y)+($F*$z);
      $atom->[$i]->{zcoor}=($G*$x)+($H*$y)+($I*$z);
    }
  }
}

## method: resetResidueName(resname,resnum[,chain])
## sets the residue name for a given residue number

sub resetResidueName {
  my $self=shift;
  my $resname=shift;
  my $resnum=shift;
  my $chain=shift;

  my $c=$self->activeChains($chain)->[0];
  my $r=$self->getResidueInChain($resnum,$c);
  for (my $ai=$r->{start}; $ai<=$r->{end}; $ai++) {
    $c->{atom}->[$ai]->{resname}=$resname;
  }
  $r->{name}=$resname;
} 

## method: resetResidueNameAll(oldresname,newresname)
## change the residue name for all occurances of a resname

sub resetResidueNameAll {
  my $self=shift;
  my $oldresname=shift;
  my $newresname=shift;

  foreach my $c ( @{$self->activeChains()} ) {
      if ($#{$c->{res}}>=0) {
	  foreach my $r ( @{$c->{res}} ) {
	      $self->resetResidueName($newresname,$r->{num},$c->{id})
		  if ($oldresname eq $r->{name});
	  }
      }
  }

  return;
}


## method: fromSICHO(sequence,chain)
## generates a molecule structure from a SICHO chain file

sub fromSICHO {
  my $self=shift;
  my $seq=shift;
  my $sicho=shift;

  my $nsicho=$#{$sicho->{sidechain}};
  my $nseq=$#{$seq->{sequence}};
  
  my $coff=0;
  $coff=1 if ($nsicho == $nseq+2);

  my $natom=1;

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{ssbond}=();

  my $chainrec=$self->_newChain("");

  for (my $i=0; $i<=$nseq; $i++) {
    my $resname=$seq->{sequence}->[$i]->{residue};
    my $resnum=$seq->{sequence}->[$i]->{index};

    my $resrec={ name => $resname, num => $resnum, start => $natom-1, chain=>"", valid=>1 };
    push (@{$chainrec->{res}},$resrec);
   
    my ($x,$y,$z)=$sicho->fromProjection($i+$coff);
    my $asc={ atominx => $natom++, atomname => "SC", 
	      resname => $resname, resnum  => $resnum, chain=> "",
              xcoor   => $x, ycoor => $y, zcoor => $z };
    push (@{$chainrec->{atom}},$asc);

    if (defined $sicho->{ca}->[$i+$coff]) {
      my ($cax,$cay,$caz)=$sicho->fromProjection($i+$coff,1);
      my $aca={ atominx => $natom++, atomname => "CA", 
		resname => $resname, resnum  => $resnum, chain=> "",
		xcoor => $cax, ycoor => $cay, zcoor => $caz};
      push (@{$chainrec->{atom}},$aca);
    }

    $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
  }

  $self->_coorCache();
}

## method: rebuildFromSICHO(sequence,chain[,fraglistoption,refPDB])
## rebuilds an all-atom structure from a SICHO chain and
## sequence file. For loop modeling loop/fragment residues
## and a protein template in PDB form may be given as
## extra arguments

sub rebuildFromSICHO {
  my $self=shift;
  my $seq=shift;
  my $sicho=shift;
  my $fraglistopt=shift;
  my $refpdb=shift;
  my $fixca=shift;

  my $nsicho=$#{$sicho->{sidechain}};
  my $nseq=$#{$seq->{sequence}};

  my $option="";

  die "empty SICHO chain" if ($nsicho<0);
  die "missing sequence" if ($nseq<0);

  if (defined $fraglistopt) {  
    die "Need reference PDB file" if (!defined $refpdb);
    die "Chain ID in fragment list is not supported for all-atom rebuilding"
      if ($fraglistopt =~ /[A-Za-z]/);
    $option.=" -l ".$fraglistopt;
    $option.=" ".$refpdb;
  }

  if (defined $fixca && $fixca) {
    $option.=" -fixca";
  }

  my $rebuildbin=&GenUtil::findExecutable("rebuild");
  die "cannot find rebuild executable"
    if (!defined $rebuildbin);

  local (*READ,*WRITE);
  my $pid=open2(*READ,*WRITE,"$rebuildbin $option");

  die "cannot open2" if (!defined $pid || $pid<0);

  my $coff=0;
  $coff=1 if ($nsicho == $nseq+2);

  for (my $i=0; $i<=$nseq; $i++) {
    my ($x,$y,$z)=$sicho->fromProjection($i+$coff);
    my $str=sprintf("%d %s %f %f %f",
		    $seq->{sequence}->[$i]->{index},$seq->{sequence}->[$i]->{residue},
		    $x,$y,$z);
    if (defined $sicho->{ca}->[$i+$coff]) {
      my ($cax,$cay,$caz)=$sicho->fromProjection($i+$coff,1);
      $str.=sprintf(" %f %f %f",$cax,$cay,$caz);
    }
    print WRITE $str,"\n";
  }
  close WRITE;

  $self->readPDB(\*READ); close READ;
  waitpid($pid,0);
}


## method: completeResidue()
## completes missing peptide backbone
## and side chain atoms 

sub completeResidue {
  my $self=shift;

  my $defchain=$self->{defchain};

  $self->translate("GENERIC");

  foreach my $c ( @{$self->activeChains()} ) {
    my @sicholist;
    my @backlist;
    my @scwrllist;
    my $sichonogly=0;
    my $scwrlnogly=0;
#    if ($c->{id} ne "+") {
    foreach my $r ( @{$c->{res}} ) {
      if ($r->{valid} && $r->{name} ne "TIP3" ) {
	my %have;
	for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	  $have{$c->{atom}->[$i]->{atomname}}=1;
	}

	my $rn=$r->{name};
	my $side=(&_cside($rn,\%have,qw( ALA CB )) ||
	    &_cside($rn,\%have,qw( CYS CB SG )) ||
	    &_cside($rn,\%have,qw( LEU CB CG CD1 CD2)) ||
	    &_cside($rn,\%have,qw( ILE CB CG2 CD1 )) ||
	    &_cside($rn,\%have,qw( ARG CB CG CD NE CZ NH1 NH2)) ||
	    &_cside($rn,\%have,qw( MET CB CG SD CE )) ||
	    &_cside($rn,\%have,qw( LYS CB CG CD CE NZ )) ||
	    &_cside($rn,\%have,qw( GLY CA )) ||
	    &_cside($rn,\%have,qw( SER CB OG )) ||
	    &_cside($rn,\%have,qw( TRP CB CG CD2 CE2 CE3 CD1 NE1 CZ2 CZ3 CH2 )) ||
	    &_cside($rn,\%have,qw( PRO CD CA CB CG)) ||
	    &_cside($rn,\%have,qw( HIS CB CG ND1 CD2 NE2 CE1)) ||
	    &_cside($rn,\%have,qw( HSE CB CG ND1 CD2 NE2 CE1)) ||
	    &_cside($rn,\%have,qw( HSD CB CG ND1 CD2 NE2 CE1)) ||
	    &_cside($rn,\%have,qw( HSP CB CG ND1 CD2 NE2 CE1)) ||
	    &_cside($rn,\%have,qw( GLU CB CG CD OE1 OE2)) || 
	    &_cside($rn,\%have,qw( GLN CB CG CD OE1 NE2)) || 
	    &_cside($rn,\%have,qw( ASP CB CG OD1 OD2)) || 
	    &_cside($rn,\%have,qw( ASN CB CG OD1 ND2)) || 
	    &_cside($rn,\%have,qw( TYR CB CG CD1 CE1 CD2 CE2 CZ OH)) || 
	    &_cside($rn,\%have,qw( PHE CB CG CD1 CE1 CD2 CE2 CZ)) || 
	    &_cside($rn,\%have,qw( THR CB OG1 CG2)) ||
	    &_cside($rn,\%have,qw( VAL CB CG1 CG2)) );
	my $back=(defined $have{CA} && defined $have{C} &&
		  defined $have{N} && 
		  (defined $have{O} || defined $have{OT1}));
	
	if (!$back) {
	  push(@backlist,$r->{num});
	}
	if (!$side || $rn eq "GLY") {
	  if (defined $have{CB} || $rn eq "GLY") {
	    push(@sicholist,$r->{num});
	    $sichonogly=1 if ($rn ne "GLY");
	  } 
	  if (!defined $have{CB}) {
	    push(@scwrllist,$r->{num});
	    $scwrlnogly=1 if ($rn ne "GLY");
	  } 
	}
      }
    }

#    printf STDERR "backlist %s\n",join(" ",@backlist);
#    printf STDERR "sicholist %s\n",join(" ",@sicholist);
#    printf STDERR "scwrllist %s\n",join(" ",@scwrllist);
    
#    printf STDERR "sichonogly: %d\n",$sichonogly; 
#    printf STDERR "scwrlnogly: %d\n",$scwrlnogly; 

    if ($#backlist>=0 && $#backlist != $#sicholist) {
      my $option="-backonly -fixca ";
      my $refpdb;
      my $fraglist;
      if ($#backlist+1 < $#{$c->{res}}) {
	$self->selectChain($c->{id});
	$refpdb="ref-".int(rand(1000000)).".pdb";
	$self->writePDB($refpdb,ssbond=>0);
	$fraglist=&GenUtil::fragListFromArray(\@backlist);
	$option.="-l ".&GenUtil::fragOptionFromList($fraglist)." ".$refpdb;
      }

      my $rebuildbin=&GenUtil::findExecutable("rebuild");
      die "cannot find rebuild executable"
	if (!defined $rebuildbin);

#      printf STDERR "running 1 $rebuildbin $option\n";

      local (*READ,*WRITE);
      my $pid=open2(*READ,*WRITE,"$rebuildbin $option");

      die "cannot open2" if (!defined $pid || $pid<0);

      foreach my $a ( @{$c->{atom}} ) {
	if ($a->{atomname} eq "CA") {
	  printf WRITE "%d %s %f %f %f %f %f %f\n",
	  $a->{resnum},$a->{resname},$a->{xcoor},$a->{ycoor},$a->{zcoor},
	  $a->{xcoor},$a->{ycoor},$a->{zcoor};
	}
      }

      close WRITE;

      my $tmol=Molecule->new();
      $tmol->readPDB(\*READ); close READ;
      $tmol->setValidResidues($fraglist);
      $tmol->setChain($c->{id});

      $self->merge($tmol);
      &GenUtil::remove($refpdb);
      waitpid($pid,0);
    }

    if ($#sicholist>=0 && $sichonogly) {
      my $option="-fixca ";
      my $refpdb;
      my $fraglist;
      if ($#sicholist+1 < $#{$c->{res}}) {
	$self->selectChain($c->{id});
	$refpdb="ref-".int(rand(1000000)).".pdb";
	$self->writePDB($refpdb,ssbond=>0);
	$fraglist=&GenUtil::fragListFromArray(\@sicholist);
	$option.="-l ".&GenUtil::fragOptionFromList($fraglist)." ".$refpdb;
      }

      my $rebuildbin=&GenUtil::findExecutable("rebuild");
      die "cannot find rebuild executable"
	if (!defined $rebuildbin);

#      printf STDERR "running 2 $rebuildbin $option\n";

      local (*READ,*WRITE);
      my $pid=open2(*READ,*WRITE,"$rebuildbin $option");

      die "cannot open2" if (!defined $pid || $pid<0);

      foreach my $r ( @{$c->{res}} ) {
	my $ca=undef;
	my $cb=undef;
	for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	  if ($c->{atom}->[$ia]->{atomname} eq "CA") {
	    $ca=$c->{atom}->[$ia];
	  } elsif ($c->{atom}->[$ia]->{atomname} eq "CB") {
	    $cb=$c->{atom}->[$ia];
	  }
	}

	if (defined $cb) {
	  printf WRITE "%d %s %f %f %f %f %f %f\n",
	  $ca->{resnum},$ca->{resname},
	  $ca->{xcoor}+$cb->{xcoor},$ca->{ycoor}+$cb->{ycoor},
	  $ca->{zcoor}+$cb->{zcoor},
	  $ca->{xcoor},$ca->{ycoor},$ca->{zcoor};
	} else {
	  printf WRITE "%d %s %f %f %f %f %f %f\n",
	  $ca->{resnum},$ca->{resname},$ca->{xcoor},$ca->{ycoor},$ca->{zcoor},
	  $ca->{xcoor},$ca->{ycoor},$ca->{zcoor};
	}
      }
      close WRITE;

      my $tmol=Molecule->new();
      $tmol->readPDB(\*READ); close READ;
      $tmol->setValidResidues($fraglist);
      $tmol->setChain($c->{id});

      $self->merge($tmol);
      &GenUtil::remove($refpdb);
      waitpid($pid,0);
    } 

    if ($#scwrllist>=0 && $scwrlnogly) {
      my $refseq="seq-".int(rand(1000000));
      my $option="-s $refseq";
      my $fraglist;

      my %haves;
      foreach my $s ( @scwrllist ) {
	$haves{$s}=1;
      }
      my $seqstr="";
      foreach my $r (@{$c->{res}}) {
	if (defined $haves{$r->{num}}) {
	  $seqstr.=$Sequence::_seqabbrev{$r->{name}};
	} else {
	  $seqstr.="x";
	}
      }
      
      open OUT,">$refseq";
      print OUT $seqstr;
      close OUT;
      
      my $scwrlbin=&GenUtil::findExecutable("scwrl");
      die "cannot find scwrl executable"
	if (!defined $scwrlbin);

#     printf STDERR "running $scwrlbin $option\n";

      local (*READ,*WRITE);
      my $pid=open2(*READ,*WRITE,"$scwrlbin $option");

      die "cannot open2" if (!defined $pid || $pid<0);

      $self->selectChain($c->{id});
      $self->writePDB(\*WRITE,ssbond=>0);
      close WRITE;

      my $tmol=Molecule->new();
      $tmol->readPDB(\*READ); close READ;
      $tmol->setValidResidues(&GenUtil::fragListFromArray(\@scwrllist));
      $tmol->setChain($c->{id});

      $self->merge($tmol);
      &GenUtil::remove($refseq);
      waitpid($pid,0);
    }
    #$self->writePDB("lat.$c->{id}.pdb");
  }
  $self->{defchain}=$defchain;
  return $self;
}  

sub _cside {
  my $name=shift;
  my $have=shift;
  my $cname=shift;
  my @atoms=@_;
  return 0 if ($name ne $cname);
  foreach my $a (@atoms) {
    return 0 
      if (!defined $have->{$a});
  }
  return 1;
}

## method: fixCOO()
## corrects the position of OT2/OXT if necessary

sub fixCOO {
  my $self=shift;
  my $nocterm=shift;
 
  foreach my $ch ( @{$self->activeChains()} ) {
    my $lastres=$ch->{res}->[$#{$ch->{res}}];

    my $ca;
    my $c;
    my $ot1;
    my $ot2;
    
    for (my $i=$lastres->{start}; $i<=$lastres->{end}; $i++) {
      my $aname=$ch->{atom}->[$i]->{atomname};

      $ch->{atom}->[$i]->{atomname}="O"
	if ($nocterm && $aname eq "OT1");
      $ch->{atom}->[$i]->{atomname}="OXX"
	if ($nocterm && $aname eq "OXT");

      if ($aname eq "CA") {
	$ca=$ch->{atom}->[$i];
      } elsif ($aname eq "C") {
	$c=$ch->{atom}->[$i];
      } elsif ($aname eq "O" || $aname eq "OT1") {
	$ot1=$ch->{atom}->[$i];
      } elsif ($aname eq "OXT" || $aname eq "OT2") {
	$ot2=$ch->{atom}->[$i];
      }
    }

    return if ($nocterm || !defined $ca || !defined $c || !defined $ot1);

    my $dot=0.0;
    if (defined $ot2) {
      my $dotx=$ot1->{xcoor}-$ot2->{xcoor};
      my $doty=$ot1->{ycoor}-$ot2->{ycoor};
      my $dotz=$ot1->{zcoor}-$ot2->{zcoor};
      $dot=sqrt($dotx*$dotx+$doty*$doty+$dotz*$dotz);
    } else {
      my $arec={};
      %{$arec}=%{$ch->{atom}->[$#{$ch->{atom}}]};
      $arec->{atominx}++;
      $arec->{atomname}="OXT";
      push(@{$ch->{atom}},$arec);
      $ot2=$arec;
    }

    if ($dot<1.5) {
      my $tcx=$c->{xcoor}-$ca->{xcoor};
      my $tcy=$c->{ycoor}-$ca->{ycoor};
      my $tcz=$c->{zcoor}-$ca->{zcoor};

      my $ltc=sqrt($tcx*$tcx+$tcy*$tcy+$tcz*$tcz);
      
      $tcx/=$ltc;
      $tcy/=$ltc;
      $tcz/=$ltc;
      
      my $t1x=$ot1->{xcoor}-$c->{xcoor};
      my $t1y=$ot1->{ycoor}-$c->{ycoor};
      my $t1z=$ot1->{zcoor}-$c->{zcoor};

      my $c1=($tcx*$t1x+$tcy*$t1y+$tcz*$t1z);
      
      $ot2->{xcoor}=$ot1->{xcoor}+2.0*($c1*$tcx-$t1x);
      $ot2->{ycoor}=$ot1->{ycoor}+2.0*($c1*$tcy-$t1y);
      $ot2->{zcoor}=$ot1->{zcoor}+2.0*($c1*$tcz-$t1z);
    }
  }
  $self->_coorCache();
}

## method: solvate([cutoff,[shape]])
## solvates a PDB from a pre-equilibrated water box

sub solvate {
  my $self=shift;
  my $cutoff=shift;
  my $shape=shift;
  my $fraglist=shift;

  $cutoff=9.0 if (!defined $cutoff);

  my $option="-box $ENV{MMTSBDIR}/data/water.pdb -cutoff $cutoff ";
  $option.="-$shape " if (defined $shape);
  
  my $solvatebin=&GenUtil::findExecutable("solvate");
  die "cannot find solvate executable"
    if (!defined $solvatebin);

  local (*READ,*WRITE);
  my $pid=open2(*READ,*WRITE,"$solvatebin $option -");
  
  die "cannot open2" if (!defined $pid || $pid <0);

  my $ss=$self->getSSBonds();
  $self->writePDB(\*WRITE,ssbond=>0);
  close WRITE;

  $self->readPDB(\*READ); close READ;
  $self->setSSBonds($ss);
  waitpid($pid,0);
}

## method: changeResName(list)
## changes the residue names according to the list

sub changeResName {
  my $self=shift;
  my $list=shift;

  if (defined $list && $list ne "") {
    foreach my $n ( split(/_/,$list)) {
      my ($newname,$rnum)=split(/:/,$n);
      my ($chainid,$num)=($rnum=~/([A-Za-z]*)([0-9]+)/);
      $self->resetResidueName($newname,$num,$chainid);
    }
  }
}

## method: fixHistidine(hsdlist,hselist,hsplist)
## changes the residue names to HSD and HSE respectively
## according to the lists of residues given as arguments

sub fixHistidine {
  my $self=shift;
  my $hsd=shift;
  my $hse=shift;
  my $hsp=shift;

  if (defined $hsd && $hsd ne "") {
    foreach my $n ( split(/:/,$hsd)) {
      my ($chainid,$num)=($n=~/([A-Za-z]*)([0-9]+)/);
      $self->resetResidueName("HSD",$num,$chainid);
    }
  }

  if (defined $hse && $hse ne "") {
    foreach my $n ( split(/:/,$hse)) {
      my ($chainid,$num)=($n=~/([A-Za-z]*)([0-9]+)/);
      $self->resetResidueName("HSE",$num,$chainid);
    }
  }

  if (defined $hsp && $hsp ne "") {
    foreach my $n ( split(/:/,$hsp)) {
      my ($chainid,$num)=($n=~/([A-Za-z]*)([0-9]+)/);
      $self->resetResidueName("HSP",$num,$chainid);
    }
  }
}


## method: setChain(chainid)
## sets the chain ID for the current structure

sub setChain {
  my $self=shift;
  my $chainid=shift;

  my $c=$self->activeChains()->[0];

  $self->{chainlookup}->{$c->{id}}=undef;
  $self->{defchain}=undef;
  $self->{segmentlist}=undef;

  $c->{id}=$chainid;

  foreach my $a ( @{$c->{atom}} ) {
    $a->{chain}=$chainid;
  }

  foreach my $r ( @{$c->{res}} ) {
    foreach my $s ( @{$self->{ssbond}} ) {
      if ($s->{resnum1}==$r->{num} &&
	  $s->{chain1}==$r->{chain}) {
	$s->{chain1}=$chainid;
      }
      if ($s->{resnum2}==$r->{num} &&
	  $s->{chain2}==$r->{chain}) {
	$s->{chain2}=$chainid;
      }
    }
  }

  foreach my $r ( @{$c->{res}} ) {
    $r->{chain}=$chainid;
  }

  $self->{chainlookup}->{uc $chainid}=$c;
}

## method: match(refMolecule)
## matches residues between two structures
## with different residue numbering and/or
## missing residues. This method changes
## the residue numbering of the current structure
## to match the reference structure given
## as the argument

sub match {
  my $self=shift;
  my $ref=shift;

  my $sc=$self->activeChains(shift)->[0];

  die "cannot find chain"
    if (!defined $sc);

  my $rc=$ref->activeChains($sc->{id})->[0];
  $rc=$ref->activeChains()->[0] if (!defined $rc);
  
  die "cannot match chain"
    if (!defined $rc);

  my $nself=$#{$sc->{res}}+1;
  my $nref=$#{$rc->{res}}+1; 

  my $maxshift=0;
  my $maxfirst=0;
  my $maxlast=0;
  my $maxmatch=0;
  
  for (my $shift=-$nref; $shift<=$nself+$nref; $shift++) {
    my $nmatch=0;
    my $first=undef;
    my $last=undef;

    for (my $j=0; $j<$nself; $j++) {
      my $sr=$sc->{res}->[$j];
      my $rr=($j>=$shift)?$rc->{res}->[$j-$shift]:undef;
      if (defined $rr && &_cmpResName($sr->{name},$rr->{name})) {
	$first=$j if (!defined $first);
	$last=$j;
	$nmatch++;
      } 
    } 
    if ($nmatch>$maxmatch) {
      $maxmatch=$nmatch;
      $maxshift=$shift;
      $maxfirst=$first;
      $maxlast=$last;
    }
  }    

  my $delta=$rc->{res}->[$maxfirst-$maxshift]->{num}-$sc->{res}->[$maxfirst]->{num};
  $self->shiftResNumber($delta,$sc->{id});

  return ($sc->{res}->[$maxfirst]->{num},$sc->{res}->[$maxlast]->{num});
}

## method: $index = getResidue(resnum[,chain])
## returns the residue index for a residue number. 
## A chain ID may be given as additional argument for
## multi-domain structures

sub getResidue {
  my $self=shift;
  my $inx=shift;
  my $chain=shift;

  my $c=$self->activeChains($chain)->[0];
  return $self->getResidueInChain($inx,$c);
}  

## method: $index = getResidueInChain(resnum[,chain])
## returns the residue index for a residue number. 
## A reference to a chain structure may be given as 
## additional argument for multi-domain structures

sub getResidueInChain {
  my $self=shift;
  my $inx=shift;
  my $c=shift;

  return undef if (!defined $c);

  if (!defined $c->{resinx}) {
    $c->{resinx}={};
    for (my $ir=0; $ir<=$#{$c->{res}}; $ir++) {
      my $key=$c->{res}->[$ir]->{num};
      $c->{resinx}->{$key}=$ir;
    }
  }

  return undef if (!defined $c->{resinx}->{$inx});

  return $c->{res}->[$c->{resinx}->{$inx}];
}  

## method: $distance = minDistance(ir,jr)
## returns the minimum distance between heavy atoms of residue
## structures ir and jr. "sc_only" flag requires a sidechain-sidechain
## contact, while "i_sc" flag requires a contact between a sidechain
## atom of i and any heavy atom of residue j.

sub minDistance {
  my $self=shift;
  my $ir=shift;
  my $jr=shift;
  my $sc_only=shift;
  my $i_sc=shift;

  $sc_only=0 if (! defined $sc_only);

  my $ic=$self->getChain($ir->{chain});
  my $jc=$self->getChain($jr->{chain});

  my $ia=$ic->{atom};
  my $ja=$jc->{atom};

  my $ix=$ic->{xcoor};
  my $iy=$ic->{ycoor};
  my $iz=$ic->{zcoor};

  my $jx=$jc->{xcoor};
  my $jy=$jc->{ycoor};
  my $jz=$jc->{zcoor};

  my $mind=1E99;
  for (my $i=$ir->{start}; $i<=$ir->{end}; $i++) {
    my $ai=$ia->[$i];
    if ((!$ai->{hyd}) && ((! $ai->{bb}) || ((! $sc_only) && (! $i_sc)))) {
      for (my $j=$jr->{start}; $j<=$jr->{end}; $j++) {
        my $aj=$ja->[$j];
	if ((!$aj->{hyd}) && ((! $sc_only) || (! $aj->{bb}))) {
	  my $dx=$ix->[$i]-$jx->[$j];
	  my $dy=$iy->[$i]-$jy->[$j];
	  my $dz=$iz->[$i]-$jz->[$j];
	  my $d=($dx*$dx+$dy*$dy+$dz*$dz);
	  $mind=$d if ($d<$mind);
	}
      }
    }
  }
  return sqrt($mind);
}

## method: number = firstResNum([chain])
## returns the number of the first residue. If
## a chain reference is given as argument it returns the
## number of the first residue of that chain

sub firstResNum {
  my $self=shift;
  my $chain=shift;

  $chain=$self->activeChains()->[0]
    if (!defined $chain);

  return undef if (!defined $chain);

  return $chain->{res}->[0]->{num};
}

## method: number = lastResNum([chain])
## returns the number of the last residue. If
## a chain reference is given as argument it returns the
## number of the last residue of that chain.

sub lastResNum {
  my $self=shift;
  my $chain=shift;

  $chain=$self->activeChains()->[0]
    if (!defined $chain);

  return undef if (!defined $chain);

  return $chain->{res}->[$#{$chain->{res}}]->{num};
}

## method: ret = empty()
## returns 1 or 0 whether the current molecule
## is empty

sub empty {
  my $self=shift;
  return (!defined $self->{chain} || $#{$self->{chain}}<0);
}

sub _coorCache {
  my $self=shift;
  my $chainid=shift;

  my $chains=$self->activeChains($chainid);
  if (defined $chains) {
    foreach my $c ( @{$self->activeChains($chainid)} ) {
      my $a=$c->{atom};
      for (my $i=0; $i<=$#{$a}; $i++) {
	$c->{xcoor}->[$i]=$a->[$i]->{xcoor};
	$c->{ycoor}->[$i]=$a->[$i]->{ycoor};
	$c->{zcoor}->[$i]=$a->[$i]->{zcoor};
      }
    }
  }
}

sub _pdbLine {
  my $pdbrec=shift;
  my $chmode=shift;
  my $longaux2=shift;
  $chmode=0 if (!defined $chmode);

  my $chainid=$pdbrec->{chain};
  $chainid=" " if (!defined $chainid || $chainid eq "" || $chainid eq "+");

  my $resnumstr;
  if ($pdbrec->{resnum}>999 && $chmode) {
    if ($pdbrec->{resnum}>9999) {
      $resnumstr=sprintf("%5d",$pdbrec->{resnum});
    } else {
      $resnumstr=sprintf(" %4d ",$pdbrec->{resnum});
    }
  } else {
    $resnumstr=sprintf("%4d  ",$pdbrec->{resnum});
  }

  my $aux1=1.0;
  my $aux2=0.0;

  $aux1=$pdbrec->{aux1} if (defined $pdbrec->{aux1});
  $aux2=$pdbrec->{aux2} if (defined $pdbrec->{aux2});

  my $atomnum;
  if ($pdbrec->{chain} eq "+") {
    $atomnum=sprintf("HETATM%5d",$pdbrec->{atominx});
  } else {
    $atomnum=sprintf("ATOM %6d",$pdbrec->{atominx});
  }
  if (length($pdbrec->{atomname})>3 || ($pdbrec->{atomname}=~/[0-9]H.*/)) {
    if (defined $longaux2 && $longaux2) {
      return sprintf "%11s %-4s %-4s%1s%-6s  %8.3f%8.3f%8.3f %5.2f %6.3f     %-4s",
	$atomnum,$pdbrec->{atomname}, 
	  $pdbrec->{resname},$chainid,$resnumstr, 
	    $pdbrec->{xcoor}, $pdbrec->{ycoor}, $pdbrec->{zcoor}, $aux1, $aux2, $pdbrec->{seg};
    } else {
      return sprintf "%11s %-4s %-4s%1s%-6s  %8.3f%8.3f%8.3f %5.2f %5.2f      %-4s",
	$atomnum,$pdbrec->{atomname}, 
	  $pdbrec->{resname},$chainid,$resnumstr, 
	    $pdbrec->{xcoor}, $pdbrec->{ycoor}, $pdbrec->{zcoor}, $aux1, $aux2, $pdbrec->{seg};
    }
  } else {
    if (defined $longaux2 && $longaux2) {
      return sprintf "%11s  %-3s %-4s%1s%-6s  %8.3f%8.3f%8.3f %5.2f %6.3f     %-4s",
	$atomnum,$pdbrec->{atomname}, 
	  $pdbrec->{resname},$chainid,$resnumstr, 
	    $pdbrec->{xcoor}, $pdbrec->{ycoor}, $pdbrec->{zcoor}, $aux1, $aux2, $pdbrec->{seg};
    } else {
      return sprintf "%11s  %-3s %-4s%1s%-6s  %8.3f%8.3f%8.3f %5.2f %5.2f      %-4s",
	$atomnum,$pdbrec->{atomname}, 
	  $pdbrec->{resname},$chainid,$resnumstr, 
	    $pdbrec->{xcoor}, $pdbrec->{ycoor}, $pdbrec->{zcoor}, $aux1, $aux2, $pdbrec->{seg};
    }      
   } 
}

sub _newChain {
  my $self=shift;
  my $cid=shift;

#  printf "new Chain >%s<\n",$cid;

  my $chainrec={};
  $chainrec->{id}=($cid eq " " || $cid eq "")?"":uc $cid;
  $chainrec->{atom}=();
  $chainrec->{res}=();
  $chainrec->{xcoor}=();
  $chainrec->{ycoor}=();
  $chainrec->{zcoor}=();
  push(@{$self->{chain}},$chainrec);

  $self->{chainlookup}->{$chainrec->{id}}=$chainrec
    if ($chainrec->{id} ne "");

  return $chainrec;
}

sub _cmpResName {
  my $name1=shift;
  my $name2=shift;
  return ($name1 eq $name2 || 
	  ($name1 =~ /HIS|HSD|HSE|HSP/ && $name2 =~ /HIS|HSD|HSE|HSP/));
}

# Try renaming His residues by looking for proton(s)
sub _correctHisType {
  my $self=shift;

  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{res}}>=0) {
      my $ra=$c->{res};
      foreach my $r ( @{$ra} ) {
	if ($r->{name}=~/HSD|HSE|HIS|HSP/) {

	  my $hisType="HIS";
	  for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	    my $aname=$c->{atom}->[$i]->{atomname};
	    if ($aname eq "HD1") {
	      if ($hisType eq "HSE") {
		$hisType="HSP";
	      } else {
		$hisType="HSD";
	      }
	    } elsif ($aname eq "HE2") {
	      if ($hisType eq "HSD") {
		$hisType="HSP";
	      } else {
		$hisType="HSE";
	      }
	    }
	  }

	  for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	    $c->{atom}->[$i]->{resname}=$hisType;
	    }
	  $r->{name}=$hisType;
        }
      }
    }
  }

  return;
}


## method: setBfactorsFromAux3()
## set Bfactors to Aux3 field

sub setBfactorsFromAux3 {
  my $self = shift;

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{atom} } >= 0 ) {
      foreach my $a ( @{ $c->{atom} } ) {
        $a->{aux2}=$a->{aux3};
      }
    }
  }

  return;
}

## method: setAux3fromBfactors()
## set Aux3 field to Bfactors

sub setAux3fromBfactors {
  my $self = shift;

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{atom} } >= 0 ) {
      foreach my $a ( @{ $c->{atom} } ) {
        $a->{aux3}=$a->{aux2};
      }
    }
  }

  return;
}

## method: operateOnAllAtoms(operator,[noH],[scOnly])
## given an operator, apply it to every atom
## note: the operator should expect an atom to be passed to it

sub operateOnAllAtoms {
  my $self = shift;
  my $operator = shift;
  my $noH = shift;
  my $scOnly = shift;

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{atom} } >= 0 ) {
      foreach my $a ( @{ $c->{atom} } ) {
        if ( ( ( !$noH ) || ( !$a->{hyd} ) )
            && ( ( !$scOnly ) || ( !$a->{bb} ) ) ) {
          &$operator($a);
        }
      }
    }
  }

  return;
}

## method: zapAux3()
## set the aux3 field to zero for every atom
## note: demonstrates usage of operateOnAllAtoms

sub zapAux3 {
  my $self = shift;

  sub zap {
    my $a = shift;
    $a->{aux3}=0.0;
    return;
  }

  $self->operateOnAllAtoms(\&{zap});

  return;
}

## method: operateOnAllAtomsVsRef(refMol,operator,[noH],[scOnly])
## given an operator, pass it every corresponding pair of atoms
## in two structures, ie. operator(self_atom,ref_atom)

sub operateOnAllAtomsVsRef {
  my $self = shift;
  my $refMol = shift;
  my $operator = shift;
  my $noH = shift;
  my $scOnly = shift;

  foreach my $c ( @{ $self->activeChains() } ) {
    if ( $#{ $c->{res} } >= 0 ) {
      my $refC = $refMol->getChain( $c->{id} );
      foreach my $r ( @{ $c->{res} } ) {
        my $start    = $r->{start};
        my $end      = $r->{end};
        my $refres   = $refMol->getResidueInChain( $r->{num}, $c );
        my $refstart = $refres->{start};
        my $refend   = $refres->{end};
        for ( my $ai = $start ; $ai <= $end ; $ai++ ) {
          my $atom = $c->{atom}->[$ai];
          if ( ( ( !$noH ) || ( !$atom->{hyd} ) ) &&
             ( ( !$scOnly ) || ( !$atom->{bb} ) ) ) {
            my $atomname = $atom->{atomname};
            for ( my $rai = $refstart ; $rai <= $refend ; $rai++ ) {
              if ( $refC->{atom}->[$rai]->{atomname} eq $atomname ) {
                &$operator($atom,$refC->{atom}->[$rai]);
              }
            }
          }
        }
      }
    }
  }

  return;
}

## method: sutstituteAla(cid,resnum)
## put an Ala at the desired residue, using the existing CB position
## sidechain hydrogens will not be built,
## CB will be put at HA2 if WT Gly
sub substituteAla {
  my $self=shift;
  my $resnum=shift;
  my $cid=shift;
  my $new_resname=shift;

  $new_resname="ALA" if (! defined $new_resname) ;

  my $c=$self->activeChains($cid)->[0];

  my $shift=0;
  my $found=0;
  for (my $ir=0; $ir<=$#{$c->{res}}; $ir++) {
      my $r=$c->{res}->[$ir];

      if ($r->{num} == $resnum) {
	  if ($r->{name} ne "GLY") {
	      for (my $ai=$r->{start}; $ai<=$r->{end}; $ai++) {
		  my $a=$c->{atom}->[$ai];
		  $c->{atom}->[$ai]->{resname}=$new_resname;
		  if ((! $a->{bb}) && ($a->{atomname} ne "CB")) {
		      splice(@{$c->{atom}},$ai,1);
		      $ai--;
		      $r->{end}--;
		      $shift--;
		  }
	      }

	  } else {
	      my $foundHA1=0;
	      for (my $ai=$r->{start}; $ai<=$r->{end}; $ai++) {
		my $a=$c->{atom}->[$ai];
		$c->{atom}->[$ai]->{resname}=$new_resname;
		if (($a->{atomname} eq "HA1") || ($a->{atomname} eq "1HA")) {
		    $a->{atomname}="CB";
		    $a->{hyd}=0;
		    $foundHA1=1;
		} elsif (($a->{atomname} eq "HA2") || ($a->{atomname} eq "2HA")) {
		    $a->{atomname}="HN";
		}
	      }
	      die "Cannot substituteAla for Gly ".$resnum.": 1HA/HA1 not present"
		  if ($foundHA1==0);
	  }
	  $r->{name}=$new_resname;
	  $found=1;

      } elsif ($found) {
	  # Update atom lookup indices for subsequent residues
	  $r->{start}+=$shift;
	  $r->{end}+=$shift;
      }
  }
  $self->_coorCache();

  return;
}

sub getSSBonds {
  my $self=shift;

  my $ns=();
  if ($#{$self->{ssbond}}>=0) {
    foreach my $s ( @{$self->{ssbond}} ) {
      my $nrec={};
      %{$nrec}=%{$s};
      push(@{$ns},$nrec);
    }
  }
  return $ns;
}

sub setSSBonds {
  my $self=shift;
  my $ss=shift;

#  $self->{ssbond}=();
  if (defined $ss && $#{$ss}>=0) {
    foreach my $s ( @{$ss} ) { 
      my $nrec={};
      %{$nrec}=%{$s};
      push(@{$self->{ssbond}},$nrec);
    }
  }
}

sub findSSBonds {
  my $self=shift;

  my $distCut=3.0;
  my $sqDistCut=$distCut*$distCut;
  $self->{ssbond}=();

  # First collect all the Cys S atoms
  my $Slist=();

  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{atom}}>=0) {
      my $aa=$c->{atom};
      foreach my $a ( @{$aa} ) {
	  push(@{$Slist},$a) if ($a->{atomname} eq "SG");
      }
    }
  }

  # Make nearby S atoms SS bonded
  for (my $i=0; $i<$#{$Slist}; $i++) {
    my $ia=$Slist->[$i];
    my $ix=$ia->{xcoor};
    my $iy=$ia->{ycoor};
    my $iz=$ia->{zcoor};

    for (my $j=$i+1; $j<=$#{$Slist}; $j++) {
       my $ja=$Slist->[$j];
       my $jx=$ja->{xcoor};
       my $jy=$ja->{ycoor};
       my $jz=$ja->{zcoor};
       my $dx=$ix-$jx;
       my $dy=$iy-$jy;
       my $dz=$iz-$jz;
       my $d=($dx*$dx+$dy*$dy+$dz*$dz);

       if ($d<$sqDistCut) {
	 my $trec={};
	 $trec->{chain1}=$ia->{chain};
	 $trec->{resnum1}=$ia->{resnum};
	 $trec->{chain2}=$ja->{chain};
	 $trec->{resnum2}=$ja->{resnum};
	 push(@{$self->{ssbond}},$trec);
       }
    }
  }


  # Ensure none have multiple partners
  # (required because of the naive method used for pairing)
  my $usedRes={};
  foreach my $SS ( @{$self->{ssbond}} ) {
      my $Sa=$SS->{chain1}.$SS->{resnum1};
      if (! defined $usedRes->{$Sa}) {
	  $usedRes->{$Sa}=1;
      } else {
	  die "Unable to perform SS pairing";
      }
      $Sa=$SS->{chain2}.$SS->{resnum2};
      if (! defined $usedRes->{$Sa}) {
	  $usedRes->{$Sa}=1;
      } else {
	  die "Unable to perform SS pairing";
      }
  }

  # Return the number of SS bonds found
  return ($#{$self->{ssbond}}+1);
}


## method: mutlist = generateMutlist(WTmol)
## Build a list of residues which differ

sub generateMutlist {
  my $self=shift;
  my $WTmol=shift;

  my $mutlist=();
  foreach my $Nc ( @{$WTmol->activeChains()} ) {
      if (($Nc->{id} ne "+") && ($#{$Nc->{res}}>=0)) {
	  my $Mc=$self->getChain($Nc->{id});
	  if (! defined $Mc) {
	      printf STDERR "Chain mismatch, chain .%s. not present",$Nc->{id};
	      printf STDERR " in mutant PDBfile\n";
	      exit 1;
	  }

	  foreach my $Nres ( @{$Nc->{res}} ) {
	      my $Mres=$self->getResidueInChain($Nres->{num},$Mc);
	      if (! defined $Mres) {
		  printf STDERR "Residue mismatch, chain .%s., residue .%s.\n",
  		  $Nc->{id},$Nres->{num};
		  exit 1;
	      }

	      # Compare one-letter codes, to avoid His tautomer problems
	      if ($Sequence::_seqabbrev{$Nres->{name}} ne
		  $Sequence::_seqabbrev{$Mres->{name}}) {
		  my $mut={};
		  $mut->{cid}=$Nc->{id};
		  $mut->{resnum}=$Nres->{num};
		  $mut->{natres}=$Nres->{name};
		  $mut->{natoneres}=$Sequence::_seqabbrev{$Nres->{name}};
		  $mut->{mutres}=$Mres->{name};
		  $mut->{mutoneres}=$Sequence::_seqabbrev{$Mres->{name}};
		  push(@{$mutlist},$mut);
	      }
	  }
      }
  }

  return $mutlist;
}


## function: mutstr = generateMutStr($mutlist)
## Convert a mutlist to a mutstring,
##  limiting to a given chain if desired

sub generateMutStr {
  my $mutlist=shift;
  my $cid=shift;

  my @mutarr=();
  foreach my $mut ( @{$mutlist} ) {
      push(@mutarr,($mut->{natoneres}.$mut->{resnum}.$mut->{mutoneres}))
	  if ((! defined $cid) || ($mut->{cid} eq $cid));
  }

  return join(",",@mutarr);
}


## function: mutlist = mutStrToMutlist($mutstr)
## Convert a mutstring to a mutlist

sub mutStrToMutlist {
  my $mutStr=shift;

  my $mutlist=();
  my @mutarr=split(/,/,$mutStr);
  foreach my $m ( @mutarr ) {
      my $mut={};
      my $WTone=substr($m,0,1);
      $mut->{natoneres}=$WTone;
      $mut->{natres}=$Sequence::_seqlong{$WTone};
      my $mutOne=substr($m,-1,1);
      $mut->{mutoneres}=$mutOne;
      $mut->{mutres}=$Sequence::_seqlong{$mutOne};
      $mut->{resnum}=substr($m,1,-1);
      push(@{$mutlist},$mut);
  }

  return $mutlist;
}


## method: rasmolDisplay(seleStr,spacefill,title)
## Show the featured selection in rasmol

sub rasmolDisplay {
  my $self=shift;
  my $seleStr=shift;
  my $spacefill=shift;
  my $title=shift;

  my $rasmolExec=&GenUtil::findExecutable("rasmol");

  my $outscript="tmp.mut.".$$.".script";
  my $outpdb="tmp.mut.".$$.".pdb";
  my $outF=&GenUtil::getOutputFile($outscript);
  printf $outF "cartoon\n";
  printf $outF "color chain\n";
  printf $outF "wireframe off\n";
  printf $outF "select %s\n",$seleStr;
  printf $outF "center selected\n";
  printf $outF "zoom 350\n";
  printf $outF "define sele (%s)\n",$seleStr;
  printf $outF "define reg (within (6.0,%s))\n",$seleStr;
  printf $outF "select within (7.0,(%s)) and (water or *.?w??) and not wss\n",$seleStr;
  printf $outF "spacefill 120\n";
  printf $outF "color yellow\n";
  printf $outF "select within (7.0,(%s)) and not (water or *.?w??) and not wss\n",$seleStr;
  if ($spacefill) {
      printf $outF "spacefill\n";
  } else {
      printf $outF "wireframe 100\n";
  }
  printf $outF "color cpk\n";
  if (defined $title) {
      printf $outF "title\"%s\"\n",$title;
  }
  undef $outF;
  $self->writePDB($outpdb);
  system($rasmolExec." -script ".$outscript." ".$outpdb);
  &GenUtil::remove($outscript);
  &GenUtil::remove($outpdb);

  return;
}


1;
