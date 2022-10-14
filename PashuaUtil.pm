# Utilities for useful calls to Pashua (for gui windows)
# 
# 2006, John Karanicolas, Baker group, UW

package PashuaUtil;

require 5.004;

use strict;

use Pashua;
use GenUtil;


## function: exitWithPopup(message)
## pops up the desired message and exits

sub exitWithPopup {
  my $message=shift;

  my $p=PashuaUtil->new();
  $p->setTitle("Exiting...");
  $p->addText($message);
  $p->run();
  exit(1);

  return;
}


## function: [ outfileHandle, appending ] = appendWithWarning(outFname)
## given a filename, check whether the output file exists. If
## not, return an output file handle. If so, check whether the
## user would rather append to the output file, overwrite it, or quit.

sub appendWithWarning {
  my $outFname=shift;

  my $found;
  my $outfile;
  if (&GenUtil::zexists($outFname)) {
      my $p=PashuaUtil->new();
      $p->setTitle("File exists");
      my $rstr="The file \"".$outFname."\" already exists. How should we proceed?";
      my $rt=$p->addRadio($rstr,1,"Append","Overwrite","Quit");
      my $pout = $p->run();
      my $sel=$pout->{$rt};
      if ($sel eq "Quit") {
	  exit(0);
      } elsif ($sel eq "Append") {
	  $outfile=&GenUtil::getAppendFile($outFname);
	  $found=1;
	  return ( $outfile, $found );
      }
  }

  $outfile=&GenUtil::getOutputFile($outFname);
  $found=0;

  return ( $outfile, $found );
}


#######################################################################
#                                                                     #
# Code below is for the "PashuaUtil" object (used by the utils above) #
#                                                                     #
#######################################################################


sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $str = shift;

  my $self = {};
  bless($self,$class);

  $self->{title}="";
  $self->{pieces}=();

  $self->{str}="";
  $self->{str}=$str if (defined $str);

  $self->{ord}={};
  $self->{ord}->{startOrd}=ord("a");
  $self->{ord}->{lastOrd}=ord("z");
  $self->{ord}->{ordA}=ord("a");
  $self->{ord}->{ordB}=ord("a");

  return $self;
}


sub getTag {
  my $self=shift;

  my $tag=chr($self->{ord}->{ordA}).chr($self->{ord}->{ordB});

  $self->{ord}->{ordB}++;
  if ($self->{ord}->{ordA} > $self->{ord}->{lastOrd}) {
      $self->{ord}->{ordA}++;
      $self->{ord}->{ordB}=$self->{ord}->{startOrd};
  }

  return $tag;
}


sub run {
  my $self=shift;

  my $fstr="*.title = ".$self->{title}."\n\n";
  $fstr .= "*.transparency = 1.0\n\n";
  $fstr .= "*.appearance = metal\n\n";

  foreach my $p ( @{$self->{pieces}} ) {
      my $tag=shift(@{$p});
      foreach my $k ( @{$p} ) {
	  $fstr .= $tag.".".$k->{key}." = ".$k->{val}."\n";
      }
      $fstr.="\n";
  }

#  printf STDERR "jk Using Pashua string below:\n\n\n%s\n",$fstr.$self->{str};

  my %result = &Pashua::run($fstr.$self->{str});
  $self->{result} = \%result;

  return $self->{result};
}


sub setTitle {
  my $self=shift;
  my $title=shift;

  if ( defined $title ) {
      $self->{title} = $title;
  } else {
      $self->{title}="";
  }

  return;
}


sub addText {
  my $self=shift;
  my $text=shift;

  my $k=(); my $p={};
  my $tag=$self->getTag(); push(@{$k},$tag);
  $p={}; $p->{key}="type"; $p->{val}="text"; push(@{$k},$p);
  $p={}; $p->{key}="default"; $p->{val}=$text; push(@{$k},$p);

  push(@{$self->{pieces}},$k);

  return $tag;
}


sub addInputTextField {
  my $self=shift;
  my $label=shift;
  my $defaultText=shift;

  $defaultText="" if (! defined $defaultText);

  my $k=(); my $p={};
  my $tag=$self->getTag(); push(@{$k},$tag);
  $p={}; $p->{key}="type"; $p->{val}="textfield"; push(@{$k},$p);
  $p={}; $p->{key}="label"; $p->{val}=$label; push(@{$k},$p);
  $p={}; $p->{key}="default"; $p->{val}=$defaultText; push(@{$k},$p);

  push(@{$self->{pieces}},$k);

  return $tag;
}


sub addRadio {
  my $self=shift;
  my $label=shift;
  my $x=shift;
  my $y=shift;
  my $defnum=shift;
  my @opt=@_;

  my $k=(); my $p={};
  my $tag=$self->getTag(); push(@{$k},$tag);
  $p={}; $p->{key}="type"; $p->{val}="radiobutton"; push(@{$k},$p);
  $p={}; $p->{key}="label"; $p->{val}=$label; push(@{$k},$p);
  $p={}; $p->{key}="x"; $p->{val}=$x; push(@{$k},$p) if (defined $x);
  $p={}; $p->{key}="y"; $p->{val}=$y; push(@{$k},$p) if (defined $y);

  my $onum=1;
  my $defname;
  while ( $#opt >= 0) {
      my $optname=shift(@opt);
      my $p={}; $p->{key}="option"; $p->{val}=$optname; push(@{$k},$p);
      $defname=$optname if ($onum == $defnum);
      $onum++;
  }
  $p={}; $p->{key}="default"; $p->{val}=$defname; push(@{$k},$p);

  push(@{$self->{pieces}},$k);

  return $tag;
}


sub addCancelButton {
  my $self=shift;
  my $label=shift;
  my @opt=@_;

  my $k=(); my $p={};
  my $tag=$self->getTag(); push(@{$k},$tag);
  $p={}; $p->{key}="type"; $p->{val}="cancelbutton"; push(@{$k},$p);
  if ( defined $label ) {
    $p={}; $p->{key}="label"; $p->{val}=$label; push(@{$k},$p);
  }

  push(@{$self->{pieces}},$k);

  return $tag;
}



1;

