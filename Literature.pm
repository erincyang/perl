# Utilities for Literature lookups, EndNote interfacing, etc.
# 
# 2006, John Karanicolas, Baker group, UW

package Literature;

require 5.004;

use strict;
no warnings "redefine";

use GenUtil;
use UWutil;
use Citation;
use LWP::Simple;
use Text::Aspell;
use vars qw ( $lastContact $utils $email $tool $db $retmax
	      $EndNote_app $Trapeze_app $pdftohtml_bin );


BEGIN {
    $lastContact = 0;
    $EndNote_app="/Applications/EndNote\\ 9/EndNote\\ 9.0.app/";
    $Trapeze_app="/Applications/Trapeze.app/";
    $pdftohtml_bin="/Users/johnk/bin/pdftohtml";
}


## function: HTMLfromPDF(inPDF)
## call pdftohtml to convert a PDF to HTML

sub HTMLfromPDF {
  my $inPDF=shift;
  my $outHTML=shift;

  $outHTML=$inPDF.".html" if (! defined $outHTML);
  return if (&GenUtil::zexists($outHTML));

  my $genHTML= substr( $inPDF, 0, -4) . ".html";

  my $convCmd=$pdftohtml_bin." -c -i -noframes ".$inPDF;
  system($convCmd);

  &GenUtil::rename($genHTML,$outHTML)
    if ( &GenUtil::zexists($genHTML) );

  return;
}


## function: HTMLfromPDF_Trapeze(inPDF)
## call Trapeze to convert a PDF to HTML

sub HTMLfromPDF_Trapeze {
  my $inPDF=shift;
  my $outHTML=shift;

  $outHTML=$inPDF.".html" if (! defined $outHTML);
  return if (&GenUtil::zexists($outHTML));

  # Rename to a short name, since Trapeze truncates long names
  my $newPDFname="a.pdf";
  &GenUtil::rename($inPDF,$newPDFname);
  my $HTMLname=$newPDFname.".html";

  my $convCmd="open -a ".$Trapeze_app." ".$newPDFname;
  system($convCmd);

  # Wait a maximum of 10 seconds - after this, Trapeze probably failed
  my $oldSize=-1;
  my $newSize=0;
  my $waitTime=0;
  while ( ( $waitTime < 10 ) && ( ( $newSize == 0 ) || ( $oldSize != $newSize ) ) ) {
      sleep(2);
      $waitTime += 2;
      if ( &GenUtil::zexists($HTMLname) ) {
	  $oldSize=$newSize;
	  $newSize = (-s $HTMLname);
      }
  }
  sleep(1);

  &GenUtil::rename($newPDFname,$inPDF);
  &GenUtil::rename($HTMLname,$outHTML)
      if ( &GenUtil::zexists($HTMLname) );

  return;
}


## function: title = titleFromHTML(infile)
## given an HTML file generated by pdftohtml, extract the title

sub titleFromHTML {
  my $infile=shift;

  return "" if ( ! &GenUtil::zexists($infile));

  my $title="";

  my $currsize=0;
  my $bestsize=0;

  # The font size lookup hash
  my $fs={};

  my $inF=&GenUtil::getInputFile($infile);
  while (<$inF>) {
    my $inline=$_;
    chomp($inline);

    # If this line defines a new font code, store this code
    if ( $inline =~  /\.ft([0-9]+).font-size.([0-9]+)px;/ ) {
      # Sample line with a font code:
      #   .ft5{font-size:5px;font-family:Times;color:#282526;}
      $fs->{$1}=$2;
    }

    # If this line contains a font code, set the new size
    if ( $inline =~  /<span class=\"ft([0-9]+?)\">/ ) {
      $currsize = $fs->{$1};
    }

    if ( $currsize > $bestsize ) {
      # If this is the largest potential title yet, clear the previous title
      # Replace breaks with spaces
      $inline=~s/<br>/ /;
      # Remove all other tags
      $inline=~s/<.*?>//g;
      if (_validTitle($inline)) {
	$bestsize = $currsize;
	$title=$inline;
      }
    } elsif ( $currsize == $bestsize ) {
      # Append to the title
      # Replace breaks with spaces
      $inline=~s/<br>/ /;
      # Remove all other tags
      $inline=~s/<.*?>//g;
      if (_validTitle($inline)) {
	$title.=" ".$inline;
      }
    }

  }
  undef $inF;

  printf STDERR "no valid title found in %s\n",$infile
    if ($title eq "");

  return $title;
}


## function: title = titleFromTrapezeHTML(infile)
## given an HTML file generated by Trapeze, extract the title

sub titleFromTrapezeHTML {
  my $infile=shift;

  return "" if ( ! &GenUtil::zexists($infile));

  # Read in the first page of the file
  my @lines=();
  my $inF=&GenUtil::getInputFile($infile);
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

  my $HTMLline=join(" ",@lines);
  my @sizeArr=split(/size\=/,$HTMLline);

  my $prevsize=0;
  my $contents={};
  my $largestSize=-1;

  # Parse out content between size tags
  while ($HTMLline =~ /^(.*?)(size\=[0-9]+)(.*?)(size\=[0-9]+)(.*)$/) {

    # Save the remainder for next time
    $HTMLline = $1.$4.$5;
    # Keep the current size tag and contents
    my $newSizeTag=$2;
    my $newContents=$3;
    $newSizeTag=~s/size\=//;

    # Remove inner tags
    while ($newContents=~/^(.*?)(\<.*?\>)(.*?)$/) {
      $newContents=$1.$3;
    }
    $newContents=~/^(.*?\>)(.*?)(\<.*?)$/;
    my $textLine=$2;

    # Remove "&nbsp;" insertions
    $textLine =~ s/\&nbsp\;//g;

    $contents->{$newSizeTag}=() if (! defined $contents->{$newSizeTag});

    # If the size tag matches the last size tag, append to the previous entry
    if ( $prevsize == $newSizeTag) {
      my $lastinx=$#{$contents->{$newSizeTag}};
      $contents->{$newSizeTag}->[$lastinx] .= " ".$textLine;
    } else {
      push(@{$contents->{$newSizeTag}},$textLine);
      $prevsize=$newSizeTag;
      $largestSize=$newSizeTag if ($newSizeTag > $largestSize);
    }

  }

  if ( $largestSize == -1) {
    printf STDERR "no size tag detected in %s\n",$infile;
    return "";
  }

  my $titleFound=0;
  for (my $currSize=$largestSize; $currSize>=0; $currSize--) {
    if (! $titleFound) {
      while ($#{$contents->{$currSize}} >= 0) {
	my $title = shift(@{$contents->{$currSize}});

	# Remove sets of 4 or more particular known characters
	$title =~ s/\.{4,}//;
	$title =~ s/\>{4,}//;
	$title =~ s/\<{4,}//;

	# Remove leading spaces, trim sets of consecutive spaces
	$title =~ s/^\ *//;
	$title =~ s/\ {2,}/\ /;

	if (_validTitle($title)) {
	  # Clean up issues with non-alphanumeric characters
	  # JK expand this section as we encounter errors (ie. train the script)
	  # Note: there's a similar stretch of text in "Citation.pm",
	  #       in htmlAuthorMatch()
	  $title =~ s/\&\#252\;/u/g;
	  $title =~ s/\&\#[0-9]+\;/\-/;
	  $title =~ s/\&amp\;//;
	  return $title;
	}
      }
    }
  }

  printf STDERR "no valid title found in %s\n",$infile;
  return "";
}


## function: title = askForTitle(PDFname)
## given a PDF filename, open the PDF then use Pashua
## to ask the user to enter the title, which is returned

sub askForTitle {
  my $PDFname=shift;

  (my $oname=$PDFname) =~ s/ /\\ /g;
  system("/usr/bin/open ".$oname);

  my @c=split(/\//,$PDFname);
  my $finalName=pop(@c);

  my $p=PashuaUtil->new();
  $p->setTitle("Unknown title");
  my $rstr="Could not extract the title from \"".$finalName.
      "\". Please enter it manually below:";
  my $rt=$p->addInputTextField($rstr);
  my $cb=$p->addCancelButton();
  my $pout = $p->run();

  return "" if ($pout->{$cb});
  return $pout->{$rt};
}


## function: _validTitle()
## exclude a few types of lines from being valid titles

sub _validTitle {
    my $title=shift;

    (my $nospace_title = $title) =~ s/ //g;

    # Eliminate anything less than 12 characters (excluding whitespace)
    return 0 if (length($nospace_title) < 10);

    # Eliminate specific non-title text that is often in a large font
    # JK expand this section as we encounter errors (ie. train the script)
    return 0 if (lc($nospace_title) =~ /\.\.\.\.\./);
    return 0 if (lc($nospace_title) =~ /\>\>\>\>\>/);
    return 0 if (lc($nospace_title) =~ /\<\<\<\<\</);
    return 0 if (lc($nospace_title) =~ "publishinggroup");
    return 0 if (lc($nospace_title) =~ "reviewarticles");
    return 0 if (lc($nospace_title) eq "newenglandjournalofmedicine");
    return 0 if (lc($nospace_title) eq "commentary");
    return 0 if (lc($nospace_title) eq "letterstonature");
    return 0 if (lc($nospace_title) eq "researcharticles");
    return 0 if (lc($nospace_title) eq "bioinformaticsapplicationsnote");
    return 0 if (lc($nospace_title) eq "sciencescope");
    return 0 if (lc($nospace_title) eq "thejournalofbiologicalchemistry");
    return 0 if (lc($nospace_title) eq "bmcmolecularbiology");
    return 0 if (lc($nospace_title) eq "minireview");
    return 0 if (lc($nospace_title) eq "perspective");
    return 0 if (lc($nospace_title) eq "advancesinimmunology");
    return 0 if ((lc($nospace_title) =~ "news") && (lc($nospace_title) =~ "views"));

    return 1;
}


## function: setupMedlineGenParams()
## setup general parts of Medline "esearch" and "efetch" calls
## Note: code adapted from the NCBI Entrez "EUtils example" script

sub setupMedlineGenParams {

    $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
#    $utils = "http://www.ncbi.nlm.nih.gov/entrez";
#    $utils = "http://eutils.ncbi.nlm.nih.gov";
    $email = &UWutil::reportEmailAddress();
    $tool = "jklit";
    $db="pubmed";
    $retmax = 99;

    return;
}


## function: $search_result = medlineTitleSearch(query)
## call Medline to evaluate a query, return a "search_outhash", which is a
## hash with keys: "count" (int), "queryKey" (int), and "webEnv" (string).
## Note: code adapted from the NCBI Entrez "EUtils example" script

sub medlineTitleSearch {
  my $query=shift;

  my @w=split(' ',$query);
  while ($#w > 8) {
      pop(@w);
  }
  $query=join(" ",@w);
  printf STDERR "Query is: %s\n",$query;

  &setupMedlineGenParams();

  $query = &URLencode($query);

  my $esearch = "$utils/esearch.fcgi?db=$db&".
      "retmax=$retmax&field=titl&usehistory=y&term=".$query;

#  my $esearch = "$utils/esearch.fcgi?db=$db&".
#      "email=$email&tool=$tool&retmax=$retmax&usehistory=y&term=".$query;

  # Run the search
   my $esearch_result = &contactNCBI($esearch);

  # Parse the output, return search_outhash
  $esearch_result =~ 
      m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;
  my $search_outhash={};
  $search_outhash->{count} = $1;
  $search_outhash->{queryKey} = $2;
  $search_outhash->{webEnv} = $3;

  return $search_outhash;
}


## function: $search_result = medlinePost(PMID)
## call Medline to lookup a given PMID, return a "post_outhash", which is a
## hash with keys: "queryKey" (int), and "webEnv" (string). This is the same
## as a "search_outhash" except that there's no "count" field
## Note: code adapted from the NCBI Entrez "EUtils example" script

sub medlinePost {
  my $PMID=shift;

#  printf STDERR "PMID is %s\n",$PMID;

  &setupMedlineGenParams();

  my $epost = "$utils/epost.fcgi?" .
      "db=$db&usehistory=y&id=".$PMID;

  # Run the post
   my $epost_result = &contactNCBI($epost);

  # Parse the output, return post_outhash
  $epost_result =~ 
      m|<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;
  my $post_outhash={};
  $post_outhash->{queryKey} = $1;
  $post_outhash->{webEnv} = $2;

  return $post_outhash;
}


## function: $mOutput = medlineFetch(search_outhash)
## retrieve citations matching a previous search result
## Note: code adapted from the NCBI Entrez "EUtils example" script

sub medlineFetch {
  my $search_outhash=shift;

  my $queryKey=$search_outhash->{queryKey};
  my $webEnv=$search_outhash->{webEnv};

  &setupMedlineGenParams();

  my $efetch = "$utils/efetch.fcgi?" .
      "db=$db&retstart=0&retmax=$retmax&rettype=medline&" .
      "retmode=text&query_key=$queryKey&WebEnv=$webEnv";

  # Run the search
  return &contactNCBI($efetch);
}


## function: $search_outhash = titleSearch(searchString)
## given a title, search medline

sub titleSearch {
  my $searchString=shift;

  # Run the search
  my $search_outhash = &medlineTitleSearch($searchString);

  if ($search_outhash->{count} > 0) {
      return $search_outhash;

  } else {
      # Broaden the search slightly, try again
      # Note: we can account for errors in "Trapeze" conversion here...
      my $newSearchString=$searchString;

      # Try removing non-alphanumerics
      $newSearchString =~ s/\,|\.|\*|\!|\?//g;

      if (($newSearchString eq $searchString)) {

	  # If there's a word which isn't in the dictionary, try to build a word
	  # by concatenating with previous or following words
	  my $speller = Text::Aspell->new();
	  my @p=split(' ',$searchString);

	  for (my $i=0; $i<=$#p; $i++) {
	      if (! $speller->check($p[$i])) {
		  my $testword=$p[$i];
		  for (my $j=$i+1; $j<=$#p; $j++) {
		      $testword.=$p[$j];
		      if ($speller->check($testword)) {
			  $p[$i]=$testword;
			  splice(@p,$i+1,$j-$i);
			  $i--;
			  $j=$#p+1;
		      }
		  }
	      }
	  }

	  for (my $i=0; $i<=$#p; $i++) {
	      if (! $speller->check($p[$i])) {
		  my $testword=$p[$i];
		  for (my $j=$i-1; $j>=0; $j--) {
		      $testword=$p[$j].$testword;
		      if ($speller->check($testword)) {
			  $p[$j]=$testword;
			  splice(@p,$j+1,$i-$j);
			  $i--;
			  $j=0;
		      }
		  }
	      }
	  }

	  $newSearchString=join(" ",@p);
      }

      # If fixing based on spell checking didn't work,
      # try sending only the words that were in the dictionary
      if (($newSearchString eq $searchString)) {
	  my $speller = Text::Aspell->new();
	  my @p=split(' ',$searchString);
	  my @w=();
	  while ($#p >= 0) {
	      my $testword=shift(@p);
	      push(@w,$testword) if $speller->check($testword);
	  }
	  $newSearchString=join(" ",@w);
      }

      if (($newSearchString eq $searchString)) {

	  # If there's any words which CAN be put together, even if they are both
	  # complete words of their own, put them together
	  my $speller = Text::Aspell->new();
	  my @p=split(' ',$searchString);

	  for (my $i=0; $i<=$#p; $i++) {
	      my $testword=$p[$i];
	      for (my $j=$i+1; $j<=$#p; $j++) {
		  $testword.=$p[$j];
		  if ($speller->check($testword)) {
		      $p[$i]=$testword;
		      splice(@p,$i+1,$j-$i);
		      $i--;
		      $j=$#p+1;
		  }
	      }
	  }

	  $newSearchString=join(" ",@p);
      }

      # Try sending only words of more than 4 letters
      if (($newSearchString eq $searchString)) {
	  my @p=split(' ',$searchString);
	  my @w=();
	  while ($#p >= 0) {
	      my $testword=shift(@p);
	      push(@w,$testword) if (length($testword) > 3);
	  }
	  $newSearchString=join(" ",@w);
      }

      # If fixing based on spell checking didn't work,
      # try pulling off one word at a time until we have a match,
      # never sending less than 4 words
      if (($newSearchString eq $searchString)) {
	  my @p=split(' ',$searchString);
	  if ($#p > 3) {
	      pop(@p);
	      $newSearchString=join(" ",@p);
	  }
      }

#      printf STDERR "search string is %s\n",$newSearchString;

      # Make sure not to do the same search twice in a row, this would be recursive
      if (($newSearchString ne $searchString)) {
	  return &titleSearch($newSearchString);
      }
  }

  return;
}


## function: $citationArray = citationArrFromQuery($search_outhash)
## given a previous search, return an array of citations for each hit

sub citationArrFromQuery {
  my $search_outhash=shift;

  my $cArr=();
  return $cArr if ($search_outhash->{count} > $retmax);
  return $cArr if ($search_outhash->{count} <= 0);

  my $mOutput = &medlineFetch($search_outhash);

  my $gotCitation=1;
  my $recNum=1;
  while ($gotCitation) {
      my $c=Citation->new();
      $gotCitation = $c->fillFromMedlineOutput($mOutput,$recNum);
      push(@{$cArr},$c) if ($gotCitation);
      $recNum++;
  }

  return $cArr;
}


# Contact NCBI, making sure not to do it with too high a frequency
sub contactNCBI {
    my $url=shift;

    my $t=time();
    my $elapsed=$t-$lastContact;
    if ($elapsed < 4) {
	sleep(4-$elapsed);
    }

    return get($url);
}

# Apply URL encoding - subroutine taken from the web
sub URLencode {
    my $theURL = $_[0];
   $theURL =~ s/([\W])/"%" . uc(sprintf("%2.2x",ord($1)))/eg;
   return $theURL;
}

# Remove URL encoding - subroutine taken from the web
sub URLdecode {
    my $theURL = $_[0];
    $theURL =~ tr/+/ /;
    $theURL =~ s/%([a-fA-F0-9]{2,2})/chr(hex($1))/eg;
    $theURL =~ s/<!--(.|\n)*-->//g;
    return $theURL;
}


1;

