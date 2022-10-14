# Citation package
# for Literature lookups, EndNote interfacing, etc.
#
# 2006, John Karanicolas, Baker group, UW

package Citation;

require 5.004;

use strict;
no warnings "redefine";

use GenUtil;
use UWutil;
use Literature;
use LWP::Simple;


use vars qw ( $LitDir @headerFields %_abbrevSub %_abbrevJ );

BEGIN {
    $LitDir = "/Users/johnk/Literature/";

    @headerFields=( "Author", "Year", "Title", "Journal", "Volume", "Pages", "Abstract", "URL", "Electronic Resource Number", "Research Notes");

    # Replace each of these words in a journal name if it's not the only word
    %_abbrevSub = ( "Biophysical" => "Biophys", Nature => "Nat", Protein => "Prot" );

    # Put a period after each of these in a journal name
    %_abbrevJ = ( Acad => 1, Am => 1, Amer => 1, Biophys => 1, Biol => 1,
		  Biotechnol => 1, Chem => 1, Curr => 1, Eng => 1,
		  Engl => 1, Eur => 1, Euro => 1, Genet => 1,
		  Immunol => 1, Lett => 1,
		  Mol => 1, Med => 1, Nat => 1, N => 1, Natl => 1,
		  Opin => 1, Phys => 1, Proc => 1, Prot => 1,
		  Rev => 1, Sci => 1, Soc => 1, Struct => 1, J => 1 );

}


## data:   PDF
## data:   HTML
## data:   TI[], AU[], DA[], etc.
##
## information about this citation


sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $PMID = shift;
  my $PDF = shift;

  my $self = {};
  bless($self,$class);

  if (defined $PDF) {
      $self->{PDF}=$PDF;
      $self->{HTML}=$PDF.".html";
  }

  $self->fillFromPMID($PMID) 
    if (defined $PMID);

  return $self;
}


# Associate a PDF with the current citation
sub associatePDF {
  my $self=shift;
  my $PDF=shift;

  $self->{PDF}=$PDF;
  if (defined $PDF) {
      $self->{HTML}=$PDF.".html";
  }

  return;
}


# fill in all fields of the current citation, given Medline output
# recNum refers to which number citation in the Medline output to read
# returns 1 if a citation was read, 0 if not
# (eg. if recNum is larger than the number of citations in the Medline output)
sub fillFromMedlineOutput {
  my $self=shift;
  my $mOutput=shift;
  my $recNum=shift;

  $recNum=1 if (! defined $recNum);

  my @mArr=split("\n",$mOutput);
  my $done=0;
  while ( ! $done ) {
      return 0 if ($#mArr < 0);
      if ($mArr[0] ne "") {
	  if ( $recNum == 1 ) {

	      # Read the current record, collect info until we get to ""
	      my $lastTag="";
	      while ($mArr[0] ne "") {
		  my $currLine=shift(@mArr);

		  # Get the type of info
		  ( my $currTag = substr($currLine,0,4) ) =~ s/ //g;
		  $currTag=$lastTag if ($currTag eq "");
		  $lastTag=$currTag;

		  # Store the info
		  $self->{$currTag}=() if (! defined $self->{$currTag});
		  my $currInfo = substr($currLine,6);
		  push(@{$self->{$currTag}},$currInfo);

#		  printf "%s: %s\n",$currTag,$currInfo;

	      }
	      return 1;

	  } elsif ( $recNum > 1 ) {

	      # Skip this record
	      my $gotLast=0;
	      while (! $gotLast) {
		  my $curr=shift(@mArr);
		  $gotLast=1 if ($curr eq "");
	      }
	      $recNum--;
	  }

      } else {
	  $done=1;
      }
  }

  return 0;
}


# fill in all fields, given a PMID
sub fillFromPMID {
  my $self=shift;
  my $PMID=shift;

  my $mOutput=&Literature::medlineFetch(&Literature::medlinePost($PMID));
  $self->fillFromMedlineOutput($mOutput);

  return;
}


# test whether the authors from this citation appear in the html file
# returns a boolean
sub htmlAuthorMatch {
  my $self=shift;

  return 0 if (! defined $self->{HTML});

  return 1 if (! &GenUtil::zexists($self->{HTML}));

  # Set the status of all to "not found"
  my $status={};
  for (my $aNum=1; $aNum <= $#{$self->{AU}}; $aNum++) {
    my $surname=lc($self->authorSurname($aNum));
    $surname=~s/ //g;
    $surname=~s/\-//g;
    $surname=~s/\.//g;
    $status->{$surname}=0;
  }

  # Mark names as "found" once we find them
  my $inF=&GenUtil::getInputFile($self->{HTML});
  while (<$inF>) {
    my $inline=$_;
    chomp($inline);
    # Replace breaks with spaces
    $inline=~s/<br>/ /g;
    # Remove all other tags
    $inline=~s/<.*?>//g;
    for (my $aNum=1; $aNum <= $#{$self->{AU}}; $aNum++) {
      my $surname=lc($self->authorSurname($aNum));
      $surname=~s/ //g;
      $surname=~s/\-//g;
      $surname=~s/\.//g;
      $inline=~s/ //g;
      $inline=~s/\-//g;
      $inline=~s/\.//g;
      if ( ! $status->{$surname} ) {
	$status->{$surname}=1 if (lc($inline) =~ $surname);
      }
    }
  }
  undef $inF;

  # Return 0 if any names are still "not found"
  for (my $aNum=1; $aNum <= $#{$self->{AU}}; $aNum++) {
    my $surname=lc($self->authorSurname($aNum));
    $surname=~s/ //g;
    $surname=~s/\-//g;
    $surname=~s/\.//g;
    return 0 if (! $status->{$surname});
  }

  return 1;
}


# test whether the authors from this citation appear in the Trapeze html file
# returns a boolean
sub htmlAuthorMatch_Trapeze {
  my $self=shift;

  return 0 if (! defined $self->{HTML});

  return 1 if (! &GenUtil::zexists($self->{HTML}));

  # Read in the first page of the html file
  my @lines=();
  my $inF=&GenUtil::getInputFile($self->{HTML});
  my $first_page=1;
  while (<$inF>) {
      $first_page=0 if ( /PAGE BREAK/ );
      my $inline=$_;
      if ($first_page) {
	  chomp($inline);
	  push(@lines,$inline);
      }
  }
  undef $inF;

  for (my $aNum=1; $aNum <= $#{$self->{AU}}; $aNum++) {

      my $surname=lc($self->authorSurname($aNum));
      my $found=0;
      foreach my $testline ( @lines ) {
	  if (! $found) {

	      # Clean up issues with non-alphanumeric characters
	      # JK expand this section as we encounter errors (ie. train the script)
	      # Note: there's a similar stretch of text in "Literature.pm",
	      #       in titleFromHTML()

	      $testline =~ s/\&\#252\;/u/g;
	      $testline =~ s/\&\#[0-9]+\;//g;
	      $testline =~ s/\&amp\;//g;
	      $testline =~ s/\&nbsp\;//g;
	      $testline =~ s/ //g;
	      $found=1 if (lc($testline) =~ $surname);
	  }
      }

      return 0 if (! $found);
  }

  return 1;
}


# test whether the title from this citation appears in the html file
# returns a boolean
sub htmlTitleMatch {
  my $self=shift;

  return 0 if (! defined $self->{HTML});

  return 1 if (! &GenUtil::zexists($self->{HTML}));

  # JK FIX CODE BELOW, THEN RETURN RESULTS FROM TESTING (INSTEAD OF JUST "1")
  return 1;

  # Read in the first page of the html file
  my $inF=&GenUtil::getInputFile($self->{HTML});
  my $first_page=1;
  my $longline="";
  while (<$inF>) {
      $first_page=0 if ( /PAGE BREAK/ );
      my $inline=$_;
      if ($first_page) {
	  chomp($inline);
	  $longline.=" ".$inline;
      }
  }
  undef $inF;

  # Clean up issues with non-alphanumeric characters
  # JK expand this section as we encounter errors (ie. train the script)
  # Note: there's a similar stretch of text in "Literature.pm",
  #       in titleFromHTML()

  $longline =~ s/\&\#252\;/u/g;
  $longline =~ s/\&\#[0-9]+\;//g;
  $longline =~ s/\&amp\;//g;
  $longline =~ s/\&nbsp\;//g;
  $longline =~ s/ //g;

  # Remove tags
  while ($longline=~/^(.*?)(\<.*?\>)(.*?)$/) {
      $longline=$1.$3;
  }

  my $title=lc($self->title());
  my @t=split(' ',$title);
  foreach my $tw ( @t ) {
      $tw="(".$tw.")";
  }
  my $pattern=join(".*?",@t);

  return 1 if ($longline =~ /$pattern/);

  return 0;
}


# return the PMID
sub PMID {
  my $self=shift;

  return $self->{PMID}->[0];
}


# return the year of publication
sub year {
  my $self=shift;

  return "" if (! defined $self->{DP});
  return substr($self->{DP}->[0],0,4);
}


# return the title of this citation
sub title {
  my $self=shift;

  return "" if (! defined $self->{TI});
  return join(" ",@{$self->{TI}});
}


# return the last name of the requested author
sub authorSurname {
  my $self=shift;
  my $aNum=shift;

  $aNum=1 if (! defined $aNum);
  $aNum--;
  return "" if (! defined $self->{AU});
  if ($aNum <= $#{$self->{AU}}) {
      my $au=$self->{AU}->[$aNum];
      my @p=split(' ',$au);
      return $p[0];
  }

  return "";
}


# move the current PDF into the "LitDir" directory, renaming it
sub storePDF {
  my $self=shift;

  my $currname=$self->{PDF};
  my $newname=$self->{PMID}->[0].".pdf";
  printf STDERR "Storing %s as %s\n",$currname,$newname;
  &GenUtil::rename($currname,$LitDir.$newname);
  $self->{PDF}=$newname;

  return;
}


# check whether this PDF is already in our list of stored entries
sub alreadyStored {
  my $self=shift;

#  my $inname=$LitDir."titles.txt";
#  my $inF=&GenUtil::getInputFile($inname);
#  while (<$inF>) {
#      my $inline=$_;
#      chomp($inline);
#      return 1 if ($inline eq $self->title());
#  }
#  undef $inF;

  return 0;
}


# remove the HTML file corresponding to the current PDF
sub removeHtml {
  my $self=shift;

  &GenUtil::remove($self->{HTML})
      if (&GenUtil::zexists($self->{HTML}));

  return;
}


# write EndNote header info
# Note: this is a function, NOT a method!
# Note: take a filehandle, NOT a filename
sub writeEndNoteHeader {
  my $outF=shift;

  printf $outF "*Journal Article\n%s\n",join("\t",@headerFields);

  return;
}


# write info about the current citation, for import into EndNote
# Note: take a filehandle, NOT a filename
sub writeEndNoteInfo {
  my $self=shift;
  my $outF=shift;

  my @infoArr=();

  foreach my $field ( @headerFields ) {
      my $val="";

      if ($field eq "URL") {
	$val="None";
	$val="file://jk:".$self->{PDF} if (defined $self->{PDF});

      } elsif ($field eq "Electronic Resource Number") {
	$val="http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed&cmd=retrieve&list_uids=".
	  $self->PMID();

      } elsif ($field eq "Research Notes") {
	  $val="Imported from Literature folder";

      } elsif ($field eq "Year") {
	  $val=$self->year();

      } elsif ($field eq "Title") {
	  $val=&_catLines($self->{TI}) if (defined $self->{TI});
	  $val=~s/\.$//;

      } elsif ($field eq "Abstract") {
	  $val=&_catLines($self->{AB}) if (defined $self->{AB});

      } elsif ($field eq "Volume") {
	  $val=&_catLines($self->{VI}) if (defined $self->{VI});

      } elsif ($field eq "Pages") {
	  $val=&_catLines($self->{PG}) if (defined $self->{PG});

	  if ($val =~ "-") {
	      # Expand page numbers from "1556-9" format to "1556-1559" format
	      my @p=split("-",$val);
	      my $startPage=shift(@p);
	      my $endPage=shift(@p);
	      if ($endPage < $startPage) {
		  $endPage=substr($startPage,0,-length($endPage)).$endPage;
	      }
	      $val=$startPage."-".$endPage;
	  }

      } elsif ($field eq "Author") {
	  $val=&_catLines($self->{AU},";") if (defined $self->{AU});

	  # Put periods between author initials
	  my @aulist=split(";",$val);
	  for (my $ai=0; $ai<=$#aulist; $ai++) {
	      my $au=$aulist[$ai];
	      my @p=split(' ',$au);
	      my $initials=pop(@p);
	      $p[$#p]=$p[$#p].",";
	      my @ilist=split("",$initials);
	      for (my $i=0; $i<=$#ilist; $i++) {
		  $ilist[$i] .= "." if ($ilist[$i] ne "-");
#		  $ilist[$i] .= " " if (($i<$#ilist) && ($ilist[$i+1] ne "-"));
	      }
	      push(@p,join("",@ilist));
	      $aulist[$ai]=join(" ",@p);
	  }
	  $val=join(";",@aulist);

      } elsif ($field eq "Journal") {
	  $val=&_catLines($self->{TA}) if (defined $self->{TA});

	  # Abbreviate journal names not already abbreviated
	  # JK expand this section as we encounter errors (ie. train the script)
	  my @p=split(' ',$val);

	  # Replace "U S A" with "USA" (for PNAS)
	  for (my $i=0; $i<=($#p-2); $i++) {
	      if (($p[$i] eq "U") && ($p[$i+1] eq "S") && ($p[$i+2] eq "A")) {
		  $p[$i]="USA";
		  splice(@p,$i+1,2);
	      }
	  }

	  foreach my $w ( @p ) {

	      # Abbreviate any recognized words, IF it's not a one-word title
	      if (($#p > 0) && (defined $_abbrevSub{$w})) {
		  $w = $_abbrevSub{$w};
	      }

	      # Put a period after each abbreviation
	      if ( defined $_abbrevJ{$w} ) {
		  $w.=".";
	      }

	  }
	  $val=join(" ",@p);
      }

      # Make sure there are no tabs in "val", then push it
      $val=~s/\t/ /;
      push(@infoArr,$val);
  }

  printf $outF "%s\n",join("\t",@infoArr);

  return;
}


# put back together info spread across multiple lines
sub _catLines {
  my $inArr=shift;
  my $delim=shift;

  $delim=" " if (! defined $delim);

  return join($delim,@{$inArr});
}


sub storedCitationFromTitle {
    my $PDFname=shift;
    my $title=shift;

    printf STDERR "Searching for title: %s\n",$title;

    my @c=split(/\//,$PDFname);
    my $finalName=pop(@c);

    my $search_outhash=&Literature::titleSearch($title);
    my $cArr=&Literature::citationArrFromQuery($search_outhash);

    if ($#{$cArr} >= 0) {

	my $cMatch;
	my $numMatch=0;
	foreach my $c ( @{$cArr} ) {
	    $c->associatePDF($PDFname);
	    # If the HTML file doesn't exist, we can't check authorship
	    # Call this a match if there's only one hit with the correct title,
	    # since it must have been manually entered
	    if ($c->htmlAuthorMatch() && $c->htmlTitleMatch()) {
		$cMatch=$c;
		$numMatch++;
	    }
	}

	if ( $numMatch == 0 ) {
	    printf STDERR "No matching citation found for hits from: %s\n",$finalName;
	} elsif ( $numMatch > 1 ) {
	    printf STDERR "Multiple matching citations found for hits from: %s\n",$finalName;
	} else {
	    # Check whether this title is on the list of titles we've already stored
	    if ( ! $cMatch->alreadyStored() ) {
		return $cMatch;
	    } else {
		printf STDERR "Found existing entry matching %s\n",$finalName;
	    }
	}

    } else {
	printf STDERR "Could not find a Medline hit for file %s\n",$finalName;
    }

    return;
}


1;

