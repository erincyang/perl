# Analyze package
# analysis functions
#
# http://mmtsb.scripps.edu/doc/Analyze.pm.html
# 2000-2002, Michael Feig, Brooks group, TSRI

package Analyze;

require 5.004;

use strict;
no warnings "redefine";

use FileHandle;
use IPC::Open2;
use Sys::Hostname;

use GenUtil;
use Molecule;
use Sequence;

## data: refmol
## Molecule object containing reference structure

## data: contactReferenceList
## reference residue contact list

## data: qReferenceList
## reference pair distance list for q scores

## constructor: new(refmol[,fraglist])
## creates a new Analyze object. A reference molecule
## is required as an argument. A fragment residue list 
## may be given in addition to restrict analysis functions
## only to those residues.

sub new {
  my $self={};

  $self->{refmol}=shift;
  $self->{fraglist}=shift;
  $self->{contactReferenceList}=undef;

  bless $self;
  return $self;
}

## method: ($number,$fraction) = gdt(cmpmol[,cutoff])
## determines the number of residues that can be matched
## between two structures with an RMSD within a given 
## distance cutoff. Arguments are the Molecule object 
## for the structure to be compared against the 
## reference structure and the value of the distance cutoff.

sub gdt {
  my $self=shift;
  my $cmp=shift;
  my $cutoff=shift;

  $cutoff=1.0 if (!defined $cutoff);
  return split(' ',$self->_rungdt($cmp,"-cutoff $cutoff"));
}

## method: $val = gdtts(cmpmol)
## calculates an approximate GDT_TS score to measure
## similarity between two structures based on a combination
## of the fractions of matching residues within a distance cutoffs
## of 1, 2, 4, and 8 A. A Molecule object for the structure
## to be compared against the reference structure is required
## as argument.

sub gdtts {
  my $self=shift;
  my $cmp=shift;
  
  return $self->_rungdt($cmp,"-ts");
}

sub _rungdt {
  my $self=shift;
  my $cmp=shift;
  my $options=shift;

  my $ref=$self->{refmol};
  
  my $gdtexe=&GenUtil::findExecutable("gdt");
  die "cannot find gdt executable"
    if (!defined $gdtexe);

  my $fname=hostname."-gdt$$";
  if (!defined $ref->{defchain}) {
    $ref->selectChain("");
    $ref->writePDB($fname,ssbond=>0);
    $ref->{defchain}=undef;
  } else {
    $ref->writePDB($fname,ssbond=>0);
  }

  local (*READ,*WRITE);
  my $pid=open2(*READ,*WRITE,"$gdtexe $options $fname");

  if (!defined $cmp->{defchain}) {
    $cmp->selectChain("");
    $cmp->writePDB(\*WRITE,ssbond=>0);
    $cmp->{defchain}=undef;
  } else {
    $cmp->writePDB(\*WRITE,ssbond=>0);
  }
  close WRITE;

  my $response=<READ>;
  chomp $response;
  close READ;

  waitpid($pid,0);

  &GenUtil::remove($fname);

  return $response;
}

## method: lsqfit(cmpmol,selmode)
## performs a least squares fit of the structure in 
## <mark>cmpmol</mark> vs. the reference structure. 
## <mark>selmode</mark> is used to select which atoms
## are used for the fit (<mark>ca</mark>, <mark>cb</mark>,
## <mark>cab</mark>, or <mark>all</mark>).

sub lsqfit {
  my $self=shift;
  my $cmp=shift;
  my $selmode=shift;
  my $warn=shift;
  my $resnumonly=shift;

  $warn=1 if (!defined $warn);

  die "need reference structure" 
    if (!defined $self->{refmol});

  undef $self->{_cmpatom};
  
  my $refatom=$self->_getRefAtom($self->{refmol},$resnumonly);
  my ($cmpatom,$keylist)=$self->_getCmpAtom($cmp,$resnumonly);
  &_translateRotate($refatom,$cmpatom,$keylist,$selmode,$warn);

  $cmp->_coorCache();
}

## method: $rmsd = rmsd(cmpmol[,warn[,keep]])
## calculates root mean square deviations of the structure in
## the argument <mark>cmpmol</mark> with the reference
## structure. 
## The method returns a data structure with the following fields
## for RMSD values for atom subsets:
## <mark>CA</mark>, <mark>CB</mark>, <mark>CAB</mark>,
## <mark>C</mark>, <mark>N</mark>, <mark>O</mark>,
## <mark>side</mark>, <mark>back</mark>, <mark>all</mark>.
## Please note that this function does not perform a
## least squares fit prior to calculating the RMSD values.

sub rmsd {
  my $self=shift;
  my $cmp=shift;
  my $warn=shift;
  my $keep=shift;
  my $resnumonly=shift;

  $warn=1 if (!defined $warn);

  undef $self->{_cmpatom} unless (defined $keep && $keep);

  die "need reference structure" 
    if (!defined $self->{refmol});

  my $refatom=$self->_getRefAtom($self->{refmol},$resnumonly);
  my ($cmpatom,$keylist)=$self->_getCmpAtom($cmp,$resnumonly);
  
  my $rms={};
  my %nrms;

  my $alt1rms=0.0;
  my $alt2rms=0.0;

#    printf STDERR "key: $key, akey: $akey, defined: %d\n",(defined $refatom->{$key})?1:0;

  for (my $ik=0; $ik<=$#{$keylist}; $ik++) {
    my $key=$keylist->[$ik]->{main};
    my $akey=$keylist->[$ik]->{alt};

    if (!defined $refatom->{$key} || !defined $refatom->{$akey}) {
      printf STDERR "cannot find matching atom for %s\n",$key
	if ($key !~ /OXT/ && $warn);
    } else {
      my $aname=$cmpatom->{$key}->{atomname};
      my $rname=$cmpatom->{$key}->{resname};
      my $rnum=$cmpatom->{$key}->{resnum};
      my $chain=$cmpatom->{$key}->{chain};
      my $ctag=(defined $chain)?"_$chain":"_";

      my $dx=$cmpatom->{$key}->{xcoor}-$refatom->{$key}->{xcoor};
      my $dy=$cmpatom->{$key}->{ycoor}-$refatom->{$key}->{ycoor};
      my $dz=$cmpatom->{$key}->{zcoor}-$refatom->{$key}->{zcoor};
      my $s=$dx*$dx+$dy*$dy+$dz*$dz;

      if ($aname eq "CA" || $aname eq "C" || $aname eq "O" ||
	  $aname eq "N" || $aname eq "CB") {

	if ($aname eq "CA") {
	  $rms->{"CA$ctag"}+=$s;   $nrms{"CA$ctag"}++;
	  $rms->{"CAB$ctag"}+=$s;  $nrms{"CAB$ctag"}++;
	}
	      
	if ($aname eq "C") {
	  $rms->{"C$ctag"}+=$s;    $nrms{"C$ctag"}++;
	}
	      
	if ($aname eq "O") {
	  $rms->{"O$ctag"}+=$s;    $nrms{"O$ctag"}++;
	}
	      
	if ($aname eq "N") {
	  $rms->{"N$ctag"}+=$s;    $nrms{"N$ctag"}++;
	}

	if ($aname eq "CB") {
	  $rms->{"CB$ctag"}+=$s;   $nrms{"CB$ctag"}++;
	  $rms->{"CAB$ctag"}+=$s;  $nrms{"CAB$ctag"}++;
	  
	  $rms->{"side$ctag"}+=$s;  $nrms{"side$ctag"}++;
	  $rms->{"side:$rname$ctag"}+=$s; $nrms{"side:$rname$ctag"}++;
	  $rms->{"side:$rnum$ctag"}+=$s; $nrms{"side:$rnum$ctag"}++;
	} else {
	  $rms->{"back$ctag"}+=$s; $nrms{"back$ctag"}++;
	  $rms->{"back:$rnum$ctag"}+=$s; $nrms{"back:$rnum$ctag"}++;
	}
	
	$rms->{"all$ctag"}+=$s;  $nrms{"all$ctag"}++;
	$rms->{"all:$rnum$ctag"}+=$s;  $nrms{"all:$rnum$ctag"}++;
      } else {
	my $s1=$s;
	my $s2=$s;
	      
	if ($akey ne $key) {
	  $dx=$cmpatom->{$akey}->{xcoor}-$refatom->{$akey}->{xcoor};
	  $dy=$cmpatom->{$akey}->{ycoor}-$refatom->{$akey}->{ycoor};
	  $dz=$cmpatom->{$akey}->{zcoor}-$refatom->{$akey}->{zcoor};
	  $s2=$dx*$dx+$dy*$dy+$dz*$dz;
	}
	
	$alt1rms+=$s1;
	$alt2rms+=$s2;
	      
	$nrms{"side$ctag"}++;
	$nrms{"side:$rname$ctag"}++;
	$nrms{"side:$rnum$ctag"}++;
	$nrms{"all$ctag"}++;
	$nrms{"all:$rnum$ctag"}++;
      }

      if ($ik==$#{$keylist} || 
	  $cmpatom->{$keylist->[$ik+1]->{main}}->{resnum} != $rnum) {
	if ($alt1rms<=$alt2rms) {
	  $rms->{"side$ctag"}+=$alt1rms;
	  $rms->{"side:$rname$ctag"}+=$alt1rms;
	  $rms->{"side:$rnum$ctag"}+=$alt1rms;
	  $rms->{"all$ctag"}+=$alt1rms;
	  $rms->{"all:$rnum$ctag"}+=$alt1rms;
	} else {
	  $rms->{"side$ctag"}+=$alt2rms;
	  $rms->{"side:$rname$ctag"}+=$alt2rms;
	  $rms->{"side:$rnum$ctag"}+=$alt2rms;
	  $rms->{"all$ctag"}+=$alt2rms;
	  $rms->{"all:$rnum$ctag"}+=$alt2rms;
	}
	$alt1rms=$alt2rms=0.0;
      }	
    }
  }

  my @kn=keys %nrms;
  my @kk=sort map(/all(_.*)/,@kn);

  foreach my $k ( @kk ) {
    foreach my $t ( map(/(.*[a-zA-Z]+)$k/,@kn) ) {
      if (defined $nrms{"$t$k"}) {
	if (!defined $nrms{$t}) {
	  $rms->{$t}=$rms->{"$t$k"};
	  $nrms{$t}=$nrms{"$t$k"};
	} else {
	  $rms->{$t}+=$rms->{"$t$k"};
	  $nrms{$t}+=$nrms{"$t$k"};
	}
      }
    }
  }

  foreach my $k (keys %nrms) {
    if (defined $nrms{$k} && $nrms{$k}>0) {
      $rms->{$k}=sqrt($rms->{$k}/$nrms{$k});
    }
  }

  return $rms;
}


## method: $list = pairwiseList([mol])
## determines a C-alpha based pairwise residue contact list 
## for the calculation of q scores

sub pairwiseList {
  my $self=shift;
  my $mol=shift;

  $mol=$self->{refmol} unless (defined $mol);
  die "no molecule defined" if (!defined $mol);

  my $clist={};

  foreach my $c ( @{$mol->activeChains()} ) {
    my $cid=$c->{id};
    my $a=$c->{atom};
    for (my $ir=0; $ir<$#{$c->{res}}; $ir++) {
      my $r1=$c->{res}->[$ir];
      my $rn1=$r1->{num};
      my $ca1;

      
    CA1:
      for (my $i=$r1->{start}; $i<=$r1->{end}; $i++) {
	if ($a->[$i]->{atomname} eq "CA") {
	  $ca1=$a->[$i];
	  last CA1;
	}
      }
      for (my $jr=$ir+1; $jr<=$#{$c->{res}}; $jr++) {
	my $r2=$c->{res}->[$jr];
	my $rn2=$r2->{num};
	my $ca2;
      CA2:
	for (my $i=$r2->{start}; $i<=$r2->{end}; $i++) {
	  if ($a->[$i]->{atomname} eq "CA") {
	    $ca2=$a->[$i];
	    last CA2;
	  }
	}
	
	my $dx=$ca1->{xcoor}-$ca2->{xcoor};
	my $dy=$ca1->{ycoor}-$ca2->{ycoor};
	my $dz=$ca1->{zcoor}-$ca2->{zcoor};
	my $d=sqrt($dx*$dx+$dy*$dy+$dz*$dz);

	$clist->{"$cid:$rn1:$rn2"}=$d;
	$clist->{"$cid:$rn2:$rn1"}=$d;
      }
    }
  }
  return $clist;
}


## method: $qscore = qscore(cmpmol[,keep])
## calculates Q scores for the structure in
## the argument <mark>cmpmol</mark> with respect to the reference
## structure. 
## The method returns a data structure with the following fields:
## <mark>all</mark>, <mark>short</mark>, <mark>medium</mark>,
## <mark>long</mark>
## Please note that this function does not perform a
## least squares fit prior to calculating the Q scores

sub qscore {
  my $self=shift;
  my $cmp=shift;
  my $keep=shift;

  undef $self->{qReferenceList} unless (defined $keep && $keep);

  die "need reference structure" 
    if (!defined $self->{refmol});

  my $ref=$self->{refmol};

  die "need contact list from reference structure" 
    if (!defined $ref && !defined $self->{qReferenceList});

  $self->{qReferenceList}=$self->pairwiseList($ref) 
    if (!defined $self->{qReferenceList});

  my $rlist=$self->{qReferenceList};

  my $qsc={};
  $qsc->{all}=0.0;
  $qsc->{short}=0.0;
  $qsc->{medium}=0.0;
  $qsc->{long}=0.0;

  my $nq=0;
  my $nqshort=0;
  my $nqmed=0;
  my $nqlong=0;

  foreach my $c ( @{$cmp->activeChains()} ) {
    my $cid=$c->{id};
    my $a=$c->{atom};
    for (my $ir=0; $ir<$#{$c->{res}}-1; $ir++) {
      my $r1=$c->{res}->[$ir];
      my $rn1=$r1->{num};
      my $ca1;
    CA1:
      for (my $i=$r1->{start}; $i<=$r1->{end}; $i++) {
	if ($a->[$i]->{atomname} eq "CA") {
	  $ca1=$a->[$i];
	  last CA1;
	}
      }
      for (my $jr=$ir+2; $jr<=$#{$c->{res}}; $jr++) {
	my $r2=$c->{res}->[$jr];
	my $rn2=$r2->{num};
	my $ca2;
      CA2:
	for (my $i=$r2->{start}; $i<=$r2->{end}; $i++) {
	  if ($a->[$i]->{atomname} eq "CA") {
	    $ca2=$a->[$i];
	    last CA2;
	  }
	}
	
	my $dx=$ca1->{xcoor}-$ca2->{xcoor};
	my $dy=$ca1->{ycoor}-$ca2->{ycoor};
	my $dz=$ca1->{zcoor}-$ca2->{zcoor};
	my $d=sqrt($dx*$dx+$dy*$dy+$dz*$dz);

	die "cannot find matching reference atom for $cid:$rn1:$rn2\n"
	  if (!exists($rlist->{"$cid:$rn1:$rn2"}));
	
	my $di=$rn2-$rn1;
	my $wid=1.0/(1.0*exp(log($di)*0.15));
	my $del=($d-$rlist->{"$cid:$rn1:$rn2"})*$wid;

	my $q=exp(-$del*$del*0.5);
	$qsc->{all}+=$q;
	$nq++;
	if ($di<5) {
	  $qsc->{short}+=$q;
	  $nqshort++;
	} elsif ($di<13) {
	  $qsc->{medium}+=$q;
	  $nqmed++;
	} else {
	  $qsc->{long}+=$q;
	  $nqlong++;
	}
      }
    }
  }

  $qsc->{all}/=$nq;
  $qsc->{short}/=$nqshort;
  $qsc->{medium}/=$nqmed;
  $qsc->{long}/=$nqlong;

  return $qsc;
}


sub _getKey {
  my $atom=shift;
  my $newname=shift;
  my $chainflag=shift;
  my $noresname=shift;

  my $resname;
  ($resname=$atom->{resname})=~s/HSE|HSD|HSP|HIE|HID/HIS/;

  $resname="" if (defined $noresname && $noresname);

  if (!defined $chainflag || $chainflag) {
    return sprintf("%s%1d%s:%s",$resname,$atom->{resnum},$atom->{chain},
		   (defined $newname)?$newname:$atom->{atomname});
  } else {
    return sprintf("%s%1d:%s",$resname,$atom->{resnum},
		   (defined $newname)?$newname:$atom->{atomname});
  }
}

sub _getRefAtom {
  my $self=shift;
  my $mol=shift;
  my $resnumonly=shift;

  if (!defined $self->{_refatom}) {
    $self->{_refatom}={};

    $self->{_multiref}=($#{$mol->activeChains()}>0)?1:0;

    foreach my $c ( @{$mol->activeChains()} ) {
      foreach my $a ( @{$c->{atom}} ) {
	$self->{_refatom}->{&_getKey($a,undef,undef,$resnumonly)}=$a;
	$self->{_refatom}->{&_getKey($a,undef,0,$resnumonly)}=$a
	    if (!$self->{_multiref});
      }
    }
  }

  return $self->{_refatom};
}

sub _getCmpAtom {
  my $self=shift;
  my $mol=shift;
  my $resnumonly=shift;

  if (!defined $self->{_cmpatom} || !defined $self->{_keylist}) {

    $self->{_cmpatom}={};
    $self->{_keylist}=();

    foreach my $c ( @{$mol->activeChains()} ) {
      foreach my $r ( @{$c->{res}} ) {
	for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	  my $a=$c->{atom}->[$ia];
	  
	  my $key=&_getKey($a,undef,undef,$resnumonly);

	  $self->{_cmpatom}->{$key}=$a;
	  
	  if ($r->{valid}) {
	    my $rname=$a->{resname};
	    my $rnum=$a->{resnum};
	    my $aname=$a->{atomname};
	    
	    if ($rname!~/TIP3|HOH/ && $aname!~/^[0-9]*H/) {
	      my $keyrec={};
	      $keyrec->{alt}=$keyrec->{main}=$key;
	      
	      if ($rname eq "ARG" && $aname eq "NH1") {
		$keyrec->{alt}=&_getKey($a,"NH2",undef,$resnumonly);
	      } elsif ($rname eq "ARG" && $aname eq "NH2") {
		$keyrec->{alt}=&_getKey($a,"NH1",undef,$resnumonly);
	      } elsif ($rname eq "ASP" && $aname eq "OD1") {
		$keyrec->{alt}=&_getKey($a,"OD2",undef,$resnumonly);
	      } elsif ($rname eq "ASP" && $aname eq "OD2") {
		$keyrec->{alt}=&_getKey($a,"OD1",undef,$resnumonly);
	      } elsif ($rname eq "GLU" && $aname eq "OE1") {
		$keyrec->{alt}=&_getKey($a,"OE2",undef,$resnumonly);
	      } elsif ($rname eq "GLU" && $aname eq "OE2") {
		$keyrec->{alt}=&_getKey($a,"OE1",undef,$resnumonly);
	      } elsif ($rname eq "LEU" && $aname eq "CD1") {
		$keyrec->{alt}=&_getKey($a,"CD2",undef,$resnumonly);
	      } elsif ($rname eq "LEU" && $aname eq "CD2") {
		$keyrec->{alt}=&_getKey($a,"CD1",undef,$resnumonly);
	      } elsif ($rname eq "VAL" && $aname eq "CG1") {
		$keyrec->{alt}=&_getKey($a,"CG2",undef,$resnumonly);
	      } elsif ($rname eq "VAL" && $aname eq "CG2") {
		$keyrec->{alt}=&_getKey($a,"CG1",undef,$resnumonly);
	      } elsif ($rname eq "PHE" && $aname eq "CD1") {
		$keyrec->{alt}=&_getKey($a,"CD2",undef,$resnumonly);
	      } elsif ($rname eq "PHE" && $aname eq "CD2") {
		$keyrec->{alt}=&_getKey($a,"CD1",undef,$resnumonly);
	      } elsif ($rname eq "PHE" && $aname eq "CE1") {
		$keyrec->{alt}=&_getKey($a,"CE2",undef,$resnumonly);
	      } elsif ($rname eq "PHE" && $aname eq "CE2") {
		$keyrec->{alt}=&_getKey($a,"CE1",undef,$resnumonly);
	      } elsif ($rname eq "TYR" && $aname eq "CD1") {
		$keyrec->{alt}=&_getKey($a,"CD2",undef,$resnumonly);
	      } elsif ($rname eq "TYR" && $aname eq "CD2") {
		$keyrec->{alt}=&_getKey($a,"CD1",undef,$resnumonly);
	      } elsif ($rname eq "TYR" && $aname eq "CE1") {
		$keyrec->{alt}=&_getKey($a,"CE2",undef,$resnumonly);
	      } elsif ($rname eq "TYR" && $aname eq "CE2") {
		$keyrec->{alt}=&_getKey($a,"CE1",undef,$resnumonly);
	      }
	      push (@{$self->{_keylist}},$keyrec);
	    }
	  }
	}
      }
    }
  }
  return ($self->{_cmpatom},$self->{_keylist});
}

sub _translateRotate {
  my $refatom=shift;
  my $cmpatom=shift;
  my $keylist=shift;
  my $selmode=shift;
  my $warn=shift;

  $warn=1 if (!defined $warn);

  die "no residues selected" if ($#{$keylist}<0);

  $selmode="cab" if (!defined $selmode || $selmode eq "");
  
  $selmode=lc $selmode;

  my $search;
  $search=":(CA|CB)\$" if ($selmode eq "cab");
  $search=":CA\$"      if ($selmode eq "ca");
  $search=":CB\$"      if ($selmode eq "cb");
  $search=":O\$" if ($selmode eq "backbone_hbonders");
  $search=":.*\$"      if ($selmode eq "all" || $selmode eq "heavy");

  my ($cxref,$cyref,$czref);
  my ($cxcmp,$cycmp,$czcmp);
  $cxref=$cyref=$czref=$cxcmp=$cycmp=$czcmp=0.0;

  my $nref=0;
  my @dokey;

  for (my $ik=0; $ik<=$#{$keylist}; $ik++) {
    $dokey[$ik]=0;

    my $key=$keylist->[$ik]->{main};

    if ($key=~/$search/) {
      if (!defined $refatom->{$key} || !defined $cmpatom->{$key}) {
	printf STDERR "cannot match atom $key\n" if ($warn);
      } else {
	$dokey[$ik]=1;

	$cxref+=$refatom->{$key}->{xcoor};
	$cyref+=$refatom->{$key}->{ycoor};
	$czref+=$refatom->{$key}->{zcoor};

	$cxcmp+=$cmpatom->{$key}->{xcoor};
	$cycmp+=$cmpatom->{$key}->{ycoor};
	$czcmp+=$cmpatom->{$key}->{zcoor};

	$nref++;
      }
    }
  }

  $cxref/=$nref;
  $cyref/=$nref;
  $czref/=$nref;
  
  $cxcmp/=$nref;
  $cycmp/=$nref;
  $czcmp/=$nref;

  my $r=();
  for (my $i=1; $i<=3; $i++) {
    $r->[$i]=();
    for (my $j=1; $j<=3; $j++) {
      $r->[$i]->[$j]=0.0;
    }
  }

  for (my $ik=0; $ik<=$#{$keylist}; $ik++) {
    if ($dokey[$ik]) {
      my $key=$keylist->[$ik]->{main};
      my $cx=$cmpatom->{$key}->{xcoor}-$cxcmp;
      my $cy=$cmpatom->{$key}->{ycoor}-$cycmp;
      my $cz=$cmpatom->{$key}->{zcoor}-$czcmp;
      my $rx=$refatom->{$key}->{xcoor}-$cxref;
      my $ry=$refatom->{$key}->{ycoor}-$cyref;
      my $rz=$refatom->{$key}->{zcoor}-$czref;

      $r->[1]->[1]+=$cx*$rx;
      $r->[2]->[1]+=$cx*$ry;
      $r->[3]->[1]+=$cx*$rz;

      $r->[1]->[2]+=$cy*$rx;
      $r->[2]->[2]+=$cy*$ry;
      $r->[3]->[2]+=$cy*$rz;

      $r->[1]->[3]+=$cz*$rx;
      $r->[2]->[3]+=$cz*$ry;
      $r->[3]->[3]+=$cz*$rz;
    }
  }

  my $u=_frotu($r);

  foreach my $k ( keys %{$cmpatom} ) {
    my $cmp=$cmpatom->{$k};
    my $x=$cmp->{xcoor}-$cxcmp;
    my $y=$cmp->{ycoor}-$cycmp;
    my $z=$cmp->{zcoor}-$czcmp;
    
    my $tx=$u->[1]->[1]*$x+$u->[1]->[2]*$y+$u->[1]->[3]*$z;
    my $ty=$u->[2]->[1]*$x+$u->[2]->[2]*$y+$u->[2]->[3]*$z;
    my $tz=$u->[3]->[1]*$x+$u->[3]->[2]*$y+$u->[3]->[3]*$z;
    
    $cmp->{xcoor}=$tx+$cxref;
    $cmp->{ycoor}=$ty+$cyref;
    $cmp->{zcoor}=$tz+$czref;
  }  
}    

sub _frotu {
  my $r=shift;

  my ($i,$j,$k);

  my $det=0.0;
  for ($i=1; $i<=3; $i++) {
    my $i1=$i+1;
    $i1-=3 if ($i1>3);
    my $i2=$i+2;
    $i2-=3 if ($i2>3);
    $det+=$r->[$i]->[1]*($r->[$i1]->[2]*$r->[$i2]->[3]-$r->[$i2]->[2]*$r->[$i1]->[3]);
  }

  my $ipt=0;
  my $w;

  for ($i=1; $i<=3; $i++) {
    for ($j=$i; $j<=3; $j++) {
      $ipt++;
      $w->[$ipt]=0.0;
      for ($k=1; $k<=3; $k++) {
	$w->[$ipt]+=$r->[$j]->[$k]*$r->[$i]->[$k];
      }
    }
  }

  my $trace=$w->[1]+$w->[4]+$w->[6];

  my $u=();

  if ($trace<3.0E-6) {
    for ($i=1; $i<=3; $i++) {
      $u->[$i]=();
      for ($j=1; $j<=3; $j++) {
	$u->[$i]->[$j]=0.0;
      }
      $u->[$i]->[$i]=1.0;
    }
    return $u;
  }

  my ($vec,$ev)=_diagq(3,3,$w);

  my $a=();
  $ipt=1;
  for ($i=1; $i<=3; $i++) {
    for ($j=1; $j<=3; $j++) {
      $a->[$j]->[$i]=$vec->[$ipt];
      $ipt++;
    }
  }

  for ($i=1; $i<=3; $i++) {
    $ev->[$i]=sqrt(abs($ev->[$i]));
    $ev->[$i]=1.0E-6 if ($ev->[$i]<1.0E-6);
  }

  $ev->[1]=-$ev->[1] if ($det<0.0);

  my $b=();
  for ($j=1; $j<=3; $j++) {
    $b->[$j]=();
  }

  for ($j=1; $j<=3; $j++) {
    my $evs=$ev->[$j];
    for ($i=1; $i<=3; $i++) {
      $b->[$i]->[$j]=0.0;
      for ($k=1; $k<=3; $k++) {
	$b->[$i]->[$j]+=$r->[$k]->[$i]*$a->[$k]->[$j]/$evs;
      }
    }
  }

  $det=0.0;
  for ($i=1; $i<=3; $i++) {
    my $i1=$i+1;
    $i1-=3 if ($i1>3);
    my $i2=$i+2;
    $i2-=3 if ($i2>3);
    $det+=$a->[$i]->[1]*($a->[$i1]->[2]*$a->[$i2]->[3]-$a->[$i2]->[2]*$a->[$i1]->[3]);
  }

  for ($j=1; $j<=3; $j++) {
    if (abs($ev->[$j]) <= 1.0E-6) {
      my $jp=$j+1;
      my $jq=$j+2;
      $jp-=3 if ($jp>3);
      $jq-=3 if ($jq>3);
      for ($k=1; $k<=3; $k++) {
	my $kp=$k+1;
	my $kq=$k+2;
	$kp-=3 if ($kp>3);
	$kq-=3 if ($kq>3);
	$b->[$k]->[$j]=$b->[$kp]->[$jp]*$b->[$kq]->[$jq]-$b->[$kp]->[$jq]*$b->[$kq]->[$jp];
	$b->[$k]->[$j]=-$b->[$k]->[$j] if ($det<0.0);
      }
    }

    my $c=0.0;
    for ($k=1; $k<=3; $k++) {
      $c+=$b->[$k]->[$j]*$b->[$k]->[$j];
    }

    $c=($c>1.0E-10)?1.0/sqrt($c):0.0;

    for ($k=1; $k<=3; $k++) {
      $b->[$k]->[$j]*=$c;
    }
  }

  for ($j=1; $j<=3; $j++) {
    $u->[$i]=();
    for ($i=1; $i<=3; $i++) {
      $u->[$i]->[$j]=0.0;
      for ($k=1; $k<=3; $k++) {
	$u->[$i]->[$j]+=$a->[$i]->[$k]*$b->[$j]->[$k];
      }
    }
  }

  for ($j=1; $j<=3; $j++) {
    my $c=0.0;
    for ($k=1; $k<=3; $k++) {
      $c+=$u->[$k]->[$j]*$u->[$k]->[$j];
    }
    $c=($c>1.0E-10)?1.0/sqrt($c):0.0;

    for ($k=1; $k<=3; $k++) {
      $u->[$k]->[$j]*=$c;
    }
  }

  $det=0.0;
  
  for ($i=1; $i<=3; $i++) {
    my $i1=$i+1;
    $i1-=3 if ($i1>3);
    my $i2=$i+2;
    $i2-=3 if ($i2>3);
    $det+=$u->[$i]->[1]*($u->[$i1]->[2]*$u->[$i2]->[3]-$u->[$i2]->[2]*$u->[$i1]->[3]);
  }
  
  printf STDERR "non-unitary rotation matrix, determinant: %f\n",$det
    if (abs($det-1.0)>1.0E-4);
      
  return $u;
}

sub _diagq {
  my $nx=shift;
  my $nfrqx=shift;
  my $dd=shift;

  my $vec=();
  my $ev=();

  my @a;
  my @b;
  my @p;
  my @w;
  my @ta;
  my @tb;
  my @y;
  
  my $nadd=0;

  my $eta=2.22045E-16;
  my $theta=4.4923E+307;

  my $n=$nx;
  my $nev=$nfrqx;
  my $nevadd=$nev+$nadd; 

  my $del1=$eta/100.0;
  my $delta=$eta*$eta*100.0;
  my $small=$eta*$eta/100.0;
  my $delbig=$theta*$delta/1000.0;
  my $theta1=1000.0/$theta;
  my $toler=100.0*$eta;
  my $rpower=8388608.0;
  my $rpow1=$rpower*0.50;
  my $rand1=$rpower-3.0;
  my $dunity=1.0;

  my $factor=0.0;
  my $ntot=int(($n*($n+1))/2);

  my ($i,$j,$k,$l,$m);

  for ($i=1; $i<=$ntot; $i++) {
    $factor=abs($dd->[$i]) if ($factor<abs($dd->[$i]));
  }

  if ($factor<$theta1) {
    printf STDERR "zero matrix passed in _diagq\n";

    for ($i=1; $i<=$nev; $i++) {
      $ev->[$i]=0.0;
      my $ipt=($i-1)*$n;
      for ($j=1; $j<$n; $j++) {
	$ipt++;
	$vec->[$ipt]=($i+$nadd == $j)?1.0:0.0;
      }
    }
    return ($vec,$ev);
  }

  my $ij=0;
  my $anorm=0.0;
  
  for ($i=1; $i<=$n; $i++) {
    for ($j=$i; $j<=$n; $j++) {
      $ij++;
      my $u=($dd->[$ij]/$factor)*($dd->[$ij]/$factor);
      $u*=0.5 if ($i == $j);
      $anorm+=$u;
    }
  }
  $anorm=sqrt($anorm+$anorm)*$factor;
  my $anormr=$dunity/$anorm;

  for ($i=1; $i<=$ntot; $i++) {
    $dd->[$i]*=$anormr;
  }

  my $nn=$n-1;
  my $mi=0;
  my $mi1=$n-1;

  for ($i=1; $i<=$nn; $i++) {
    my $sum1=0.0;
    $b[$i]=0.0;
    my $ji=$i+1;
    my $ipt=$mi+$i;
    $a[$i]=$dd->[$ipt];
    $ipt++;
    my $bx=$dd->[$ipt];
    my $ji2=$ji+1;
    if ($ji == $n) {
      $b[$i]=$bx;
      $dd->[$mi+$ji]=0.0;
      $mi+=$mi1;
      $mi1--;
    } else {
      for ($j=$ji2; $j<=$n; $j++) {
	$ipt++;
	$sum1+=$dd->[$ipt]*$dd->[$ipt];
      }
      
      if ($sum1<=$small) {
	$b[$i]=$bx;
	$dd->[$mi+$ji]=0.0;
	$mi+=$mi1;
	$mi1--;
      } else {
	my $s=sqrt($sum1+$bx*$bx);
	my $sgn=($bx>=0.0)?abs($dunity):-abs($dunity);
	my $temp=abs($bx);
	$w[$ji]=sqrt(0.5*($dunity+($temp/$s)));
	$ipt=$mi+$ji;
	$dd->[$ipt]=$w[$ji];
	my $ii=$i+2;
	if ($ii<=$n) {
	  $temp=$sgn/(2.0*$w[$ji]*$s);
	  for ($j=$ii; $j<=$n; $j++) {
	    $ipt++;
	    $w[$j]=$temp*$dd->[$ipt];
	    $dd->[$ipt]=$w[$j];
	  }
	}

	$b[$i]=-$sgn*$s;
	
	for ($j=$ji; $j<=$n; $j++) {
	  $p[$j]=0.0;
	}

	my $ml=$mi+$mi1;
	my $ml1=$mi1-1;


	for ($l=$ji; $l<=$n; $l++) {
	  $ipt=$ml+$l;
	  for ($m=$l; $m<=$n; $m++) {
	    $bx=$dd->[$ipt];
	    $p[$l]+=$bx*$w[$m];
	    $p[$m]+=$bx*$w[$l] if ($l!=$m);
	    $ipt++;
	  }
	  $ml+=$ml1;
	  $ml1--;
	}
	my $xkap=0.0;
	
	for ($k=$ji; $k<=$n; $k++) {
	  $xkap+=$w[$k]*$p[$k];
	}

	for ($l=$ji; $l<=$n; $l++) {
	  $p[$l]-=$xkap*$w[$l];
	}
	
	my $mj=$mi+$mi1;
	my $mj1=$mi1-1;

	for ($j=$ji; $j<=$n; $j++) {
	  for ($k=$j; $k<=$n; $k++) {
	    my $expr=($p[$j]*$w[$k])+($p[$k]*$w[$j]);
	    $dd->[$mj+$k]-=$expr+$expr;
	  }
	  $mj+=$mj1;
	  $mj1--;
	}
	
	$mi+=$mi1;
	$mi1--;
      }
    }
  }

  $a[$n]=$dd->[$mi+$n];
  $b[$n]=0.0;

  my $alimit=1.0;
  for ($i=1; $i<=$n; $i++) {
    $w[$i]=$b[$i];
    $b[$i]*=$b[$i];
  }
  
  for ($i=1; $i<=$nevadd; $i++) {
    $ev->[$i]=$alimit;
  }
  my $rootl=-$alimit;

  for ($i=1; $i<=$nevadd; $i++) {  
    my $rootx=$alimit;
    for ($j=$i; $j<=$nevadd; $j++) {
      $rootx=$ev->[$j] if ($ev->[$j]<$rootx);
    }
    $ev->[$i]=$rootx;

    my $trial=($rootl+$ev->[$i])*0.5;

    while(abs($trial-$rootl) > 1.0E-15 && abs($trial-$ev->[$i]) > 1.0E-15) {
      my $nomtch=$n;
      $j=1;
      
      do {
	my $f0=$a[$j]-$trial;
	
	while ($j<=$n && abs($f0)>=$theta1) {
	  $nomtch-- if ($f0>=0.0);
	  $j++;
	  $f0=$a[$j]-$trial-$b[$j-1]/$f0;
	}  
	if ($j<=$n) {
	  $j+=2;
	  $nomtch--;
	}
      }	while ($j<=$n);

      if ($nomtch>=$i) {
	$ev->[$i]=$trial;
	my $nom=($nevadd<=$nomtch)?$nevadd:$nomtch;
	$ev->[$nom]=$trial;
      } else {
	$rootl=$trial;
      }

      $trial=($rootl+$ev->[$i])*0.5;      
    }
  }

  for ($i=1; $i<=$nev; $i++) {
    $ev->[$i]=$ev->[$i+$nadd];
  }

  my $ia=0;
  
  for ($i=1; $i<=$nev; $i++) {
    my $aroot=$ev->[$i];
    for ($j=1; $j<=$n; $j++) {
      $y[$j]=1.0;
    }
    $ia=-1 if ($i==1 || abs($ev->[$i-1]-$aroot)>=$toler);
    $ia++;

    my $elim1=$a[1]-$aroot;
    my $elim2=$w[1];

    for ($j=1; $j<=$nn; $j++) {
      my $temp;
      if (abs($elim1)<=abs($w[$j])) {
	$ta[$j]=$w[$j];
	$tb[$j]=$a[$j+1]-$aroot;
	$p[$j]=$w[$j+1];
	$temp=(abs($w[$j])>$theta1)?$elim1/$w[$j]:1.0;
	$elim1=$elim2-$temp*$tb[$j];
	$elim2=-$temp*$w[$j+1];
      } else {
	$ta[$j]=$elim1;
	$tb[$j]=$elim2;
	$p[$j]=0.0;
	$temp=$w[$j]/$elim1;
	$elim1=$a[$j+1]-$aroot-$temp*$elim2;
	$elim2=$w[$j+1];
      }
      $b[$j]=$temp;
    }

    $ta[$n]=$elim1;
    $tb[$n]=0.0;
    $p[$n]=0.0;
    $p[$nn]=0.0;
    my $iter=1;
    
    if ($ia!=0) { 
      for ($j=1; $j<=$n; $j++) {
	my $rand1=(4099.0*$rand1 % $rpower);
	$y[$j]=$rand1/$rpow1-1.0;
      }
    }

    do {
      $l=$n+1;
      
      for ($j=1; $j<=$n; $j++) {
	$l--;
	do {
	  if (($n-$l-1)<0) {
	    $elim1=$y[$l];
	  } elsif (($n-$l-1)==0) {
	    $elim1=$y[$l]-$y[$l+1]*$tb[$l];
	  } else {
	    $elim1=$y[$l]-$y[$l+1]*$tb[$l]-$y[$l+2]*$p[$l];
	  }
	  if ($elim1>$delbig || $elim1<-$delbig) {
	    for ($k=1; $k<=$n; $k++) {
	      $y[$k]=$y[$k]/$delbig;
	    }
	  } 	  
	} while ($elim1>$delbig || $elim1<-$delbig);
	
	my $temp=$ta[$l];
	$temp=$delta if (abs($temp)<$delta);
	$y[$l]=$elim1/$temp;
      }
      
      if ($iter==1) {
	$elim1=$y[1];
	for ($j=1; $j<=$nn; $j++) {
	  if (abs($ta[$j]-$w[$j])<1E-15) {
	    $y[$j]=$y[$j+1];
	    $elim1=$elim1-$y[$j+1]*$b[$j];
	  } else {
	    $y[$j]=$elim1;
	    $elim1=$y[$j+1]-$elim1*$b[$j];
	  }
	}
	$y[$n]=$elim1;
      } 
      $iter++;
    } while ($iter<=2);
    
    my $ipt;
    if ($ia != 0) {
      for (my $j1=1; $j1<=$ia; $j1++) {
	$k=$i-$j1;
	my $temp=0.0;
	$ipt=($k-1)*$n;
	for ($j=1; $j<=$n; $j++) {
	  $ipt++;
	  $temp+=$y[$j]*$vec->[$ipt];
	}
	$ipt=($k-1)*$n;
	for ($j=1; $j<=$n; $j++) {
	  $ipt++;
	  $y[$j]-=$temp*$vec->[$ipt];
	}
      }
    }
    
    $elim1=0.0;
    for ($j=1; $j<=$n; $j++) {
      $elim1=abs($y[$j]) if (abs($y[$j])>$elim1);
    }
    my $temp=0.0;
    for ($j=1; $j<=$n; $j++) {
      $elim2=$y[$j]/$elim1;
      $temp+=$elim2*$elim2;
    }
    $temp=$dunity/(sqrt($temp)*$elim1);
    for ($j=1; $j<=$n; $j++) {
      $y[$j]*=$temp;
      $y[$j]=0.0 if (abs($y[$j])<$del1);
    }
    $ipt=($i-1)*$n;
    for ($j=1; $j<=$n; $j++) {
      $ipt++;
      $vec->[$ipt]=$y[$j];
    }
  }

  my $ipt;
  my $kk;
  for ($i=1; $i<=$nev; $i++) {
    $ipt=($i-1)*$n;
    for ($j=1; $j<=$n; $j++) {
      $ipt++;
      $y[$j]=$vec->[$ipt];
    }
    $l=$n-2;
    my $mk=($n*($n-1))/2-3;
    my $mk1=3;
    
    my $t;
    for ($j=1; $j<=$l; $j++) {
      $t=0.0;
      $k=$n-$j-1;
      $m=$k+1;
      for ($kk=$m; $kk<=$n; $kk++) {
	$t+=$dd->[$mk+$kk]*$y[$kk];
      }
      for ($kk=$m; $kk<=$n; $kk++) {
	my $epr=$t*$dd->[$mk+$kk];
	$y[$kk]-=$epr+$epr;
      }
      $mk-=$mk1;
      $mk1++;
    }

    $t=0.0;
    for ($j=1; $j<=$n; $j++) {
      $t+=$y[$j]*$y[$j];
    }
    
    my $xnorm=sqrt($t);
    my $xnorm1=$dunity/$xnorm;
    for ($j=1; $j<=$n; $j++) {
      $y[$j]*=$xnorm1;
    }
    
    $ipt=($i-1)*$n;
    for ($j=1; $j<=$n; $j++) {
      $ipt++;
      $vec->[$ipt]=$y[$j];
    }
  }

  for ($i=1; $i<=$n; $i++) {
    $ev->[$i]=$ev->[$i]*$anorm;
  }

  return ($vec,$ev);
}

## method: $list = contactList([mol[,fraglist]])
## determines a residue contact list based
## on a minimum distance between heavy atoms of
## two residues of less than 4.2 A. By default
## the contact list is calculated from the reference 
## structure, but a different structure may be
## given as an argument.

sub contactList {
  my $self=shift;
  my $mol=shift;

  $mol=$self->{refmol} unless (defined $mol);
  die "no molecule defined" if (!defined $mol);

  my $clist=();

  foreach my $c ( @{$mol->activeChains()} ) {
    for (my $ir=0; $ir<=$#{$c->{res}}; $ir++) {
      my $r1=$c->{res}->[$ir];
      for (my $jr=$ir+1; $jr<=$#{$c->{res}}; $jr++) {
	my $r2=$c->{res}->[$jr];
	if (($r1->{valid} || $r2->{valid}) && 
	    ($r2->{num}>=$r1->{num}+5 || $r1->{chain} ne $r2->{chain})) {
	  my $mind=$mol->minDistance($r1,$r2);
	  if ($mind<4.2) {
	    my $rec={};
	    $rec->{res1}=$r1->{num};
	    $rec->{res2}=$r2->{num};
	    $rec->{chain1}=$r1->{chain};
	    $rec->{chain2}=$r2->{chain};
	    $rec->{d}=$mind;
	    push (@{$clist},$rec);
	  }
	}
      }
    }    
  }
  return $clist;
}

## method: ($ncont,$fraction,$rho,$list) = contacts(cmpmol)
## compares the residue contact list in <mark>cmpmol</mark>
## with a reference contact list that has to be either loaded
## previously, or will be calculated from the reference
## structure if not available. It returns the number of
## matching contacts, the fraction with respect to all
## contacts in the reference structure, the quantity rho
## measuring the contact fraction in a more continuous manner,
## and the contact list itself.

sub contacts {
  my $self=shift;
  my $cmp=shift;

  my $ref=$self->{refmol};

  die "need contact list from reference structure" 
    if (!defined $ref && !defined $self->{contactReferenceList});

  $self->{contactReferenceList}=$self->contactList($ref) 
    if (!defined $self->{contactReferenceList});
  
  my $ncont=0;
  my $scont=0.0;

  my $contactList=();
  my $ntot=0;

  foreach my $cl ( @{$self->{contactReferenceList}} ) {
    my $r1=$cmp->getResidue($cl->{res1},$cl->{chain1});
    my $r2=$cmp->getResidue($cl->{res2},$cl->{chain2});

    if (($r1->{valid} || $r2->{valid})) {
      my $d=$cl->{d};
      
      my $mind=$cmp->minDistance($r1,$r2);
      $ncont++ if ($mind<=4.2);
      $scont+=1.0/(1.0+exp(20.0*($mind-4.45)));
      
      $ntot++;

      $cl->{compd}=$mind;
      push(@{$contactList},$cl);
    }
  }

  if ($ntot<=0) {
    return (0,0.0,0.0,$contactList);
  } else {
    return ($ncont,$ncont/$ntot,$scont/$ntot,$contactList);
  }
}

## method: readContacts(file)
## reads a residue contact list from a file to 
## be used for comparison with <mark>contacts</mark>

sub readContacts {
  my $self=shift;
  my $cfile=&GenUtil::getInputFile(shift);

  $self->{contactReferenceList}=();

  while (<$cfile>) {
    if (!/^\#/) {
      chomp;
      s/^ +//g;

      my $rec={};
      my ($v1,$v2);
      ($v1, $v2, $rec->{d})=split(/ +/);
      ($rec->{res1},$rec->{chain1})=split(/:/,$v1);
      ($rec->{res2},$rec->{chain2})=split(/:/,$v2);
      push(@{$self->{contactReferenceList}},$rec);
    }
  }

  undef $cfile;
}

## method: writeContacts(file)
## writes the current reference contact list to a file

sub writeContacts {
  my $self=shift;

  die "no contact list defined"
    unless (defined $self->{contactReferenceList});

  my $cfile=&GenUtil::getOutputFile(shift);

  foreach my $c ( @{$self->{contactReferenceList}}) {
    my $tc1=(defined $c->{chain1} && $c->{chain1} ne "")?":$c->{chain1}":"";
    my $tc2=(defined $c->{chain2} && $c->{chain2} ne "")?":$c->{chain2}":"";
    printf $cfile "%d%s %d%s %f\n",
    $c->{res1},$tc1,$c->{res2},$tc2,$c->{d};
  }

  undef $cfile;
}

## function: $value = radiusOfGyration(mol[,caonly])
## calculates the radius of gyration for a molecule

sub radiusOfGyration {
  my $mol=shift;
  my $caonly=shift;

  my ($cx,$cy,$cz);
  $cx=$cy=$cz=0.0;

  my $val=0.0;
  my $ns=0;

  foreach my $c ( @{$mol->activeChains()} ) {
    foreach my $a ( @{$c->{atom}} ) {
      if ($a->{atomname}!~/^[0-9]*H/ && (!defined $caonly || $caonly==0 || $a->{atomname} eq "CA")) {
	$cx+=$a->{xcoor};
	$cy+=$a->{ycoor};
	$cz+=$a->{zcoor};
	$ns++;
      }
    }
  }

  $cx/=$ns;
  $cy/=$ns;
  $cz/=$ns;

  foreach my $c ( @{$mol->activeChains()} ) {
    foreach my $a ( @{$c->{atom}} ) {
      if ($a->{atomname}!~/^[0-9]*H/ && (!defined $caonly || $caonly==0 || $a->{atomname} eq "CA")) {
	my $dx=$a->{xcoor}-$cx;
	my $dy=$a->{ycoor}-$cy;
	my $dz=$a->{zcoor}-$cz;

	$val+=$dx*$dx+$dy*$dy+$dz*$dz;
      }
    }
  }

  $val/=$ns;
  return sqrt($val);
}

## function: chi1(mol)
## calculates chi1 side chain dihedral angles for
## all residues in the given molecule

sub chi1 {
  my $mol=shift;

  my %inx1;
  my %inx2;
  my %inx3;
  my %inx4;

  foreach my $c ( @{$mol->activeChains()} ) {
    foreach my $a ( @{$c->{atom}} ) {
      my $key="$a->{chain}$a->{resnum}";
      $inx1{$key}=$a if ($a->{atomname} eq "N");
      $inx2{$key}=$a if ($a->{atomname} eq "CA");
      $inx3{$key}=$a if ($a->{atomname} eq "CB");
      $inx4{$key}=$a 
	if (($a->{atomname} eq "CG" && 
	     $a->{resname}=~/(GLU|PRO|LYS|GLN|ARG|LEU|PHE|TYR|TRP|ASN|ASP|MET|HIS)/) ||
	    ($a->{atomname} eq "OG" && $a->{resname}=~/(SER)/) ||
	    ($a->{atomname} eq "OG1" && $a->{resname}=~/(THR)/) ||
	    ($a->{atomname} eq "CG1" && $a->{resname}=~/(VAL|ILE)/) ||
	    ($a->{atomname} eq "SG" && $a->{resname}=~/(CYS)/));
    }
  }
    
  foreach my $c ( @{$mol->activeChains()} ) {
    foreach my $r ( @{$c->{res}} ) {
      if ($r->{valid}) {
	my $key="$r->{chain}$r->{num}";
	if (defined $inx1{$key} && defined $inx2{$key} && 
	    defined $inx3{$key} && defined $inx4{$key}) {
	  $r->{chi1}=
	    &GenUtil::dihedral($inx1{$key},$inx2{$key},$inx3{$key},$inx4{$key});
	} else {
	  $r->{chi1}=0.0;
	}
      }
    }
  }
}

## function: phipsi(mol)
## calculates phi and psi backbone dihedral angles
## for all residues in the molecule given as argument

sub phipsi {
  my $mol=shift;

  my %cinx;
  my %cainx;
  my %ninx;

  foreach my $c ( @{$mol->activeChains()} ) {
    foreach my $a ( @{$c->{atom}} ) {
      my $key="$a->{chain}$a->{resnum}";
      $cinx{$key}=$a  if ($a->{atomname} eq "C");
      $cainx{$key}=$a if ($a->{atomname} eq "CA");
      $ninx{$key}=$a  if ($a->{atomname} eq "N");
    }
  }

  foreach my $c ( @{$mol->activeChains()} ) {
    foreach my $r ( @{$c->{res}} ) {
      if ($r->{valid}) {
	my $key="$r->{chain}$r->{num}";
	my $keyp=sprintf("%s%d",$r->{chain},$r->{num}+1);
	my $keym=sprintf("%s%d",$r->{chain},$r->{num}-1);

	if (defined $ninx{$key} && defined $cainx{$key} && 
	    defined $cinx{$key}) {
	  if (defined $cinx{$keym}) {
	    $r->{phi}=
	      &GenUtil::dihedral($cinx{$keym},$ninx{$key},$cainx{$key},$cinx{$key});
	  } 
	  if (defined $ninx{$keyp}) {
	    $r->{psi}=
	      &GenUtil::dihedral($ninx{$key},$cainx{$key},$cinx{$key},$ninx{$keyp});
	  }
	  if (defined $ninx{$keyp} && defined $cainx{$keyp}) {
	    $r->{omega}=
	      &GenUtil::dihedral($cainx{$key},$cinx{$key},$ninx{$keyp},$cainx{$keyp});
	  }
	}
      }
    }
  }
}

## method: ($rmsd,$percentage) = chi1RMSD(cmpmol)
## calculates the root mean square deviation and percentage
## of correct chi1 dihedral values
## in <mark>cmpmol</mark> compared with the chi1 dihedral values
## in the reference structure. 

sub chi1RMSD {
  my $self=shift;
  my $cmp=shift;

  if (!defined $self->{_havechi1}) {
    &chi1($self->{refmol});
    $self->{_havechi1}=1;
  }

  &chi1($cmp);
  my ($chi1,$correctchi1)=$self->_angleRMSD($self->{refmol},$cmp,"chi1");
  return ($chi1,$correctchi1);
}

## method: ($phirmsd,$psirmsd,$phipercent,$psipercent) = 
## method:  phipsiRMSD(cmpmol[,fraglist])
## calculates the root mean square deviation and percentage
## of correct phi and psi backbone dihedral values
## in <mark>cmpmol</mark> compared with the phi/psi dihedral values
## in the reference structure. 

sub phipsiRMSD {
  my $self=shift;
  my $cmp=shift;
  my $fraglist=shift;

  if (!defined $self->{_havepsi} || !defined $self->{_havephi}) {
    &phipsi($self->{refmol});
    $self->{_havepsi}=$self->{_havephi}=1;
  }

  &phipsi($cmp);
  my ($phi,$correctphi)=$self->_angleRMSD($self->{refmol},$cmp,"phi");
  my ($psi,$correctpsi)=$self->_angleRMSD($self->{refmol},$cmp,"psi");
  return ($phi,$psi,$correctphi,$correctpsi);
}


### _angleRMSD ######

sub _angleRMSD {
  my $self=shift;
  my $ref=shift;
  my $cmp=shift;
  my $angle=shift;

  my $sum=0.0;
  my $nsum=0;
  my $correct=0;

  foreach my $c ( @{$cmp->activeChains()} ) {
    foreach my $r ( @{$c->{res}} ) {
      my $rr=$ref->getResidue($r->{num},$r->{chain});
      if ($r->{valid} && defined $rr && defined $r->{$angle} && $rr->{$angle}) {
	my $da=$r->{$angle}-$rr->{$angle};
	$da+=360.0 if ($da<-180.0);
	$da-=360.0 if ($da>180.0);
	$sum+=$da*$da;
	$nsum++;
	$correct++ if ($da>-40 && $da<40);
      }
    }
  }
  return ($nsum>0)?(sqrt($sum/$nsum),100.0*$correct/$nsum):(0,0);
}

## function: $hist = radialDistribution(mol,select,resolution)
## calculates the radial distribution histogram of the
## selected subset of atoms from the center of the molecule
## given in <mark>mol</mark>. The selection criteria
## may be a residue name, or one of <mark>POL</mark> 
## (polar residues), <mark>HYD</mark> (hydrophobic residues),
## <mark>CHG</mark> (charged residues), <mark>POS</mark>
## (positively charged residues), <mark>NEG</mark>
## (negatively charged residues). The resolution arguments determines
## the width of the histogram bins.

sub radialDistribution {
  my $mol=shift;
  my $sel=shift;
  my $resolution=shift;

  $sel="*" if (!defined $sel);
  $resolution=1.0 if (!defined $resolution);
  
  my ($cx,$cy,$cz)=$mol->centerOfMass();
  
  my $hist=();
  for (my $i=0; $i<=int(100.0/$resolution); $i++) {
    $hist->[$i]=0;
  }

  my $nmatch=0;

  foreach my $c ( @{$mol->activeChains()} ) {
    foreach my $a ( @{$c->{atom}} ) {
      if (&_matchSelection($sel,$a)) {
#      printf STDERR "%s:%s\n",$a->[$ia]->{resname},$a->[$ia]->{atomname};
	my $x=$a->{xcoor}-$cx;
	my $y=$a->{ycoor}-$cy;
	my $z=$a->{zcoor}-$cz;
      
	my $d=sqrt($x*$x+$y*$y+$z*$z);
      
	if ($d<100.0) {
	  my $inx=int($d/$resolution);
	  $hist->[$inx]++;
	}
	$nmatch++;
      }
    }
  }

  my $tfac=4.0*3.141592654/3.0*$nmatch/100.0;
  for (my $i=0; $i<=int(100.0/$resolution); $i++) {
    my $r1=$i*$resolution;
    my $r2=$r1+$resolution;
    $hist->[$i]/=$tfac*($r2*$r2*$r2-$r1*$r1*$r1);
  }
  
  return $hist;
}

## function: $value = radialRMSD(mol,select)
## calculates the radial RMSD of a subset of atoms 
## from the center of mass in the molecule given by
## <mark>mol</mark>. For an explanation of the
## selection argument see <mark>radialDistribution</mark>
## above.

sub radialRMSD {
  my $mol=shift;
  my $sel=shift;

  $sel="*" if (!defined $sel);
  
  my ($cx,$cy,$cz)=$mol->centerOfMass();

  my $sum=0.0;
  my $nmatch=0.0;
      
  foreach my $c ( @{$mol->activeChains()} ) {
    foreach my $a ( @{$c->{atom}} ) {
      if (&_matchSelection($sel,$a)) {
	my $x=$a->{xcoor}-$cx;
	my $y=$a->{ycoor}-$cy;
	my $z=$a->{zcoor}-$cz;
      
	my $d=sqrt($x*$x+$y*$y+$z*$z);
	$sum+=$d*$d;
	$nmatch++;
      }
    }
  }

  return ($nmatch>0)?sqrt($sum/$nmatch):0.0;
}



### _matchSelection ######

sub _matchSelection {
  my $sel=shift;
  my $atom=shift;

  foreach my $s ( split(',',$sel) ) {
    if ($s eq "none") {
    } else {
      my ($r,$a)=split(/:/,$s);
      return 1       
	if (&_matchResidue($r,$atom) && &_matchAtom($a,$atom));
    }
  }
  return 0;
}

### _matchResidue ######

sub _matchResidue {
  my $res=shift;
  my $atom=shift;
  
  return 1 if (!defined $res || $res eq "" || $res eq "*" || 
	       $res eq "ALL" || $res eq "all");

  my $cres=uc $atom->{resname};
  $cres=~s/(HSD|HSE|HSP|HID|HIE)/HIS/;
  
  foreach my $r ( split(/\+/,$res) ) { 
    $r=uc $r;
    $r=~s/(HSD|HSE|HSP|HID|HIE)/HIS/;
    
    return 1 if ($r eq $cres);
    return 1 if ($Sequence::_seqlong{$r} eq $cres);
    return 1 if ($r eq "POL" && $cres=~/(SER|THR|TYR|CYS|ASN|GLN|HIS|TRP)/);
    return 1 if ($r eq "HYD" && $cres=~/(ALA|VAL|PHE|ILE|PRO|MET|LEU)/);
    return 1 if ($r eq "CHG" && $cres=~/(ASP|GLU|LYS|ARG)/);
    return 1 if ($r eq "POS" && $cres=~/(LYS|ARG)/);
    return 1 if ($r eq "NEG" && $cres=~/(ASP|GLU)/);
  }
  return 0;
}

### _matchAtom ######

sub _matchAtom {
  my $at=shift;
  my $atom=shift;
  
  return 1 if (!defined $at || $at eq "" || $at eq "*" || $at eq "all" || $at eq "ALL");
  
  my $cat=uc $atom->{atomname};
  
  foreach my $a ( split(/\+/,$at) ) {
    $a=uc $a;
    
    return 1 if ($a eq $cat);
    return 1 if ($a eq "NOH" && $cat!~/^[0-9]*H.*/);
    return 1 if ($a eq "ALLC" && $cat=~/^[0-9]*C.*/);
    return 1 if ($a eq "ALLN" && $cat=~/^[0-9]*N.*/);
    return 1 if ($a eq "ALLO" && $cat=~/^[0-9]*O.*/);
    return 1 if ($a eq "ALLH" && $cat=~/^[0-9]*H.*/);
  }

  return 0;
}
