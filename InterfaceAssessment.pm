# InterfaceAssessment package
# for evaluating interfaces
#
# 2005, John Karanicolas, Baker lab

package InterfaceAssessment;

require 5.004;

use strict;

use GenUtil;
use PDButil;
use Molecule;
use RosettaMolecule;
use CHARMM;

use vars qw ( $debug $par $datadir $outPDBdir );

BEGIN {
    $debug=0;
    $datadir="datadir/";
    $outPDBdir="outPDBs/";
}


## constructor: new([inPDBlist])
## creates a new InterfaceAssessment object using the given inPDBlist

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $inPDBlist = shift;
  my $doClean = shift;

  my $self = {};
  $self->{inPDBlist}=$inPDBlist;
  $self->{intDat}={};
  $self->{resDat}={};
  $self->{files}={};

  bless($self,$class);

  printf STDERR "Debug mode is ON\n" if ($debug);

  &UWutil::setCondorCpus(0);

  $self->setupPDBset($doClean);

  return $self;
}


## method: setupPDBset()
## creates maps from each input PDB name to the clean complex and free PDBs

sub setupPDBset {

  my $self=shift;
  my $doClean = shift;

  $doClean=1 if (! defined $doClean);

  printf STDERR "Setting up PDBs\n";

  my $cleandir="clean_pdbs/";
  &GenUtil::makeDir($cleandir);

  my $mapName=$cleandir."PDBmap.out";

  # Reuse any existing cleaned files
  if (-e $mapName) {
      my $inmap=&GenUtil::getInputFile($mapName);
      while (<$inmap>) {
	  my $inline=$_;
	  chomp($inline);
	  my @p=split(' ',$inline);
	  $self->{cmap}->{$p[0]}=$p[1];
	  if ($#p > 1) {
	      $self->{fmap}->{$p[0]}=$p[2];
	  }
      }
      undef $inmap;
  }

  # Create entries for any that don't exist
  my $not_present=();
  foreach my $f ( @{$self->{inPDBlist}} ) {
      if ((! defined $self->{cmap}->{$f})) {
	  push(@{$not_present},$f);
      }
  }

  if ($#{$not_present} >= 0) {
      if ($doClean) {
	  $self->cleanPDBset($mapName,$cleandir,$not_present,0);
      } else {
	  # If we're not cleaning, just point to the original PDB
	  my $outmap=&GenUtil::getAppendFile($mapName);
	  $outmap->autoflush(1);
	  my $fnum=1;
	  foreach my $inPDB ( @{$not_present} ) {
	      printf $outmap "%s %s\n",$inPDB,$inPDB;
	      $self->{cmap}->{$inPDB}=$inPDB;
	  }
	  undef $outmap;
      }
  }

  return;
}


## method: cleanPDBset()
## creates maps from each input PDB name to the clean complex and free PDBs

sub cleanPDBset {

    my $self=shift;
    my $mapName=shift;
    my $cleandir=shift;
    my $not_cleaned=shift;

    printf STDERR "Creating clean PDBs\n";

    if ($#{$not_cleaned} >= 0) {
	# Create clean files for any that don't exist
	my $outmap=&GenUtil::getAppendFile($mapName);
	$outmap->autoflush(1);
	my $fnum=1;
	foreach my $inPDB ( @{$not_cleaned} ) {
	    while (&GenUtil::zexists($cleandir.$fnum.".pdb")) {
		$fnum++;
	    }
	    my $outCname=$cleandir.$fnum.".pdb";
	    $self->cleanPDB($inPDB,$outCname);
	    printf $outmap "%s %s\n",$inPDB,$outCname;
	    $self->{cmap}->{$inPDB}=$outCname;
	}
	undef $outmap;
    }

    return;
}


## method: cleanPDB()
## clean a single PDB

sub cleanPDB {

    my $self=shift;
    my $inPDB=shift;
    my $outCname=shift;

    printf STDERR "Cleaning %s\n",$inPDB;
    my $mol=RosettaMolecule->new($inPDB);

    # If there are more than two chains, look for a single TER
    # and split the chains here
    if ($mol->numchains()>2) {

	# Look for a TER, split them here
	my $foundTER=0;
	my $gotTwo=0;
	my $infile=&GenUtil::getInputFile($inPDB);
	my $tmpfile="tmp.pdb";
	my $outfile=&GenUtil::getOutputFile($tmpfile);

	while (<$infile>) {
	    if (! $gotTwo) {
		my $inline=$_;
		chomp($inline);
		if ($inline =~ "TER") {
		    $gotTwo=1 if ($foundTER);
		    $foundTER=1;
		    printf $outfile "%s\n",$inline;
		} elsif ($inline =~ /^END /) {
		    $gotTwo=1;
		} elsif ($inline =~ /^ATOM /) {
		    if ($foundTER) {
			# Second chain - call it "B"
			substr($inline,21,1,"B");
		    } else {
			# First chain - call it "A"
			substr($inline,21,1,"A");
		    }
		    printf $outfile "%s\n",$inline;
		}
	    }
	}

	printf $outfile "END\n";
	undef $infile;
	undef $outfile;
	undef $mol;
	$mol=RosettaMolecule->new($tmpfile);
    }

    $mol->writePDB($outCname,"coor");

    return;
}


## method: collectAllInterfaceData(par,[wtStruct])
## call all the required interface collection functions

sub collectAllInterfaceData {

    my $self=shift;
    $par=shift;
    my $wtStruct=shift;

    $self->{files}->{wtStruct}=$wtStruct;

    if (defined $wtStruct) {
	die "Cannot find WT struct" if (! &GenUtil::exists($wtStruct));
	my $resfile=&PDButil::prefixAndPath($wtStruct).".resfile";
	$resfile =~ s/\.native// if (! &GenUtil::zexists($resfile));
	if (&GenUtil::zexists($resfile)) {
	    $self->{files}->{resfile}=$resfile;
	    printf STDERR "Using resfile %s\n",$resfile;
	}
	printf STDERR "Collecting list of mutations\n";
	$self->collectInterfaceDataType("Mutlist.out",0,\&{_call_mutlist});
    }

    if (&_test_par("RosMode") eq "fast") {
	printf STDERR "Collecting Rosetta data\n";
	$self->collectInterfaceDataType("fastRos.out",0,\&{_call_Ros});
    }

    if (&_test_par("RosMode") eq "full") {
	printf STDERR "Collecting Rosetta data with water\n";
	$self->collectInterfaceDataType("fullRos.out",0,\&{_call_Wat});
    }

    if (&_test_par("useCha")) {
	printf STDERR "Collecting charmm energies\n";
	$self->collectInterfaceDataType("Cha.out",0,\&_call_Cha);
    }

    if ((&_test_par("useMon")) && (defined $wtStruct)) {
	printf STDERR "Collecting data on monomer stability\n";
	$self->collectInterfaceDataType("Monomers.out",0,\&{_call_Mon});
    }

    if (&GenUtil::exists($datadir."Experimental.out")) {
        printf STDERR "Checking for experimental data\n";
        $self->collectInterfaceDataType("Experimental.out",0,\&{_fail_for_expt_data});
    }

    if ((&_test_par("byRes")) && (defined $wtStruct)) {

	foreach my $inPDB ( @{$self->{inPDBlist}} ) {
	    my $mutStr=$self->{intDat}->{$inPDB}->{"Mutlist"};
	    my @mutlist=split(",",$mutStr);
	    foreach my $m ( @mutlist ) {
		$self->{resDat}->{$inPDB}->{$m}={};
		$self->{resDat}->{$inPDB}->{$m}->{resnum}=substr($m,1,-1);
	    }
	}

	printf STDERR "Collecting free stability data\n";
	$self->collectInterfaceDataType("DEScontext.out",1,\&{_call_free_dG});

	printf STDERR "Collecting data for each mutation\n";
	$self->collectInterfaceDataType("ddG.out",1,\&{_call_ddG});

	printf STDERR "Collecting WT stability diffs\n";
	$self->collectInterfaceDataType("WTcontext.out",1,\&{_call_WT_dG});

	printf STDERR "Collecting min/max data\n";
	$self->collectInterfaceDataType("ResMinMax.out",0,\&{_call_ResMinMax});

	printf STDERR "Collecting Rosetta consensus sequence\n";
	$self->collectInterfaceDataType("Consensus.out",1,\&{_call_consensus_seq});

    }

    return;
}


## method: _check_for_expt_data()
## check for (manually added) experimental data
sub _check_for_expt_data {
    my $self=shift;
    die "Cannot compute experimental data";
    return;
}


## method: _byres_setup
## transfer mutation list for usability

sub _byres_setup {

    my $self=shift;

    foreach my $inPDB ( @{$self->{inPDBlist}} ) {
	my $mutStr=$self->{intDat}->{$inPDB}->{"Mutlist"};
	my @mutlist=split(",",$mutStr);
	foreach my $m ( @mutlist ) {
	    $self->{resDat}->{$inPDB}->{$m}={};
	}
    }

    return;
}


## method: collectInterfaceDataType(datFname,defaultField_fxn,compute_fxn,struct)
## collect data of a generic type

sub collectInterfaceDataType {

    my $self=shift;
    my $datFname=shift;
    my $by_res=shift;
    my $compute_fxn=shift;
    my $struct=shift;

    my $dat={};
    &GenUtil::makeDir($datadir);
    my $fullDatName=$datadir.$datFname;

    # Reuse data where possible
    my $missing_inputs=();
    my $print_header=0;
    if (-s $fullDatName) {
	my $first_field=$self->_read_existing($fullDatName,$by_res);
	$missing_inputs=$self->_find_missing_inputs($first_field,$by_res);
    } else {
	foreach my $inPDB ( @{$self->{inPDBlist}} ) {
	    if ($by_res) {
		my $mutStr=$self->{intDat}->{$inPDB}->{"Mutlist"};
		my @mutlist=split(",",$mutStr);
		foreach my $mut ( @mutlist ) {
		    push(@{$missing_inputs},($inPDB."::".$mut));
		}
	    } else {
		push(@{$missing_inputs},$inPDB);
	    }
	}
	$print_header=1;
    }

    # Compute missing data, if required
    if ($#{$missing_inputs} >= 0) {
	my $outdat=&GenUtil::getAppendFile($fullDatName);
	$outdat->autoflush(1);
	&$compute_fxn($self,$missing_inputs,$outdat,$print_header,$struct);
	undef $outdat;
    }

    return;
}


## function: _test_par(parname)
## tests whether a parameter is set and true

sub _test_par {
    my $key=shift;
    return 0 if (! defined $par->{$key});
    return $par->{$key};
}


## method: first_field = _read_existing(fname)
## read data from an existing outfile whenever possible

sub _read_existing {

    my $self=shift;
    my $inFname=shift;
    my $by_res=shift;

    $by_res=0 if (! defined $by_res);

    my $indat=&GenUtil::getInputFile($inFname);

    # Parse the first line separately
    my $inline=<$indat>;
    chomp($inline);
    my @fields=split(' ',$inline);
    shift(@fields);
    shift(@fields) if ($by_res);

    while (<$indat>) {
	my $inline=$_;
	chomp($inline);
	my @p=split(' ',$inline);
	my $fname=shift(@p);
	if (defined $self->{cmap}->{$fname}) {
	    if (! $by_res) {
		for (my $i=0; $i<=$#fields; $i++) {
		    $self->{intDat}->{$fname}->{$fields[$i]}=shift(@p);
		}
	    } else {
		my $mut=shift(@p);
		for (my $i=0; $i<=$#fields; $i++) {
		    $self->{resDat}->{$fname}->{$mut}->{$fields[$i]}=shift(@p);
		}
	    }
	}
    }
    undef $indat;

    return $fields[0];
}


## method: _find_missing_inputs(first_field)
## list any structures for which we need to generate data
sub _find_missing_inputs {
    my $self=shift;
    my $first_field=shift;
    my $by_res=shift;

    $by_res=0 if (! defined $by_res);

    my $missing=();
    foreach my $inPDB ( @{$self->{inPDBlist}} ) {
	if (! $by_res) {
	    if (! defined $self->{intDat}->{$inPDB}->{$first_field}) {
		push(@{$missing},$inPDB);
	    }
	} else {
	    my $mutStr=$self->{intDat}->{$inPDB}->{"Mutlist"};
	    my @mutlist=split(",",$mutStr);
	    foreach my $m ( @mutlist ) {
		if (! defined $self->{resDat}->{$inPDB}->{$m}->{$first_field}) {
		    push(@{$missing},($inPDB."::".$m));
		}
	    }
	}
    }

    return $missing;
}


## method: dumpIntOutput(outFname)
## write output to a table

sub dumpIntOutput {

    my $self=shift;
    my $outFname=shift;

    my $outF=&GenUtil::getOutputFile($outFname);

    my $field_list=();
    printf $outF "FILENAME";

    foreach my $field ( keys %{$self->{intDat}->{$self->{inPDBlist}->[0]}} ) {
	push(@{$field_list},$field);
	printf $outF "  %s",$field;
    }
    printf $outF "\n";

    foreach my $f ( @{$self->{inPDBlist}} ) {
	printf $outF "%s", $f;
	foreach my $field ( @{$field_list} ) {
	    printf $outF "  %s",$self->{intDat}->{$f}->{$field};
	}
	printf $outF "\n";
    }

    return;
}



## method: writeHtml()
## write descriptions to formatted html pages

sub writeHtml {

    my $self=shift;
    my $outFname=shift;

    foreach my $inPDB ( @{$self->{inPDBlist}} ) {

	die "Interface data not defined" if (! defined $self->{intDat}->{$inPDB});

	if (! defined $outFname) {
	    $outFname=$inPDB;
	    $outFname=~s/\.gz//g;
	    $outFname=~s/\.pdb//g;
	    $outFname.=".html";
	}
	my $outF=&GenUtil::getOutputFile($outFname);

	# Print the header
	printf $outF "\<html\>\n";
	printf $outF "\<head\>\n";
	printf $outF "\<title\>\n";
	printf $outF "Interface assessment for %s\n", $inPDB;
	printf $outF "\<\/title\>\n";
	printf $outF "\<\/head\>\n\n";
	printf $outF "\<body bgcolor=\"\#FFFFCC\" text=black link=blue vlink=purple alink=red\>\n";
	printf $outF "\<center\>\n";
	printf $outF "\<h2\>\n";
	printf $outF "Interface assessment for %s\n", $inPDB;
	printf $outF "\<\/h2\>\n";
	printf $outF "\<\/center\>\n";
	printf $outF "\<hr size=5\>\n\n";

	# List the mutations
	if (defined $self->{intDat}->{$inPDB}->{"Mutlist"}) {
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    my $mutstrA=$self->{intDat}->{$inPDB}->{"MutA"};
	    $mutstrA=~s/,/, /g;
	    printf $outF "\<p\>%d mutations in chain A:  %s\n",$self->{intDat}->{$inPDB}->{"NumMutA"},$mutstrA;
	    my $mutstrB=$self->{intDat}->{$inPDB}->{"MutB"};
	    $mutstrB=~s/,/, /g;
	    printf $outF "\<p\>%d mutations in chain B:  %s\n",$self->{intDat}->{$inPDB}->{"NumMutB"},$mutstrB;
	    printf $outF "\<\/h3\>\n";
	    printf $outF "\<\/left\>\n";
	    printf $outF "\<hr size=0\>\n\n";
	}

	# We'll use these in making tables
	my $table;
	my $row=();

	# List the misc. quantifiers
	if ((defined $self->{intDat}->{$inPDB}->{"ncont"}) ||
	    (defined $self->{intDat}->{$inPDB}->{"CharmmElec"})) {
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    printf $outF "Misc. quantifiers\n";
	    $table=();
	    if (defined $self->{intDat}->{$inPDB}->{"ncont"}) {
		$row=();
		push(@{$row},"NumContacts",$self->{intDat}->{$inPDB}->{"ncont"});
		push(@{$table},$row);
	    }
	    if (defined $self->{intDat}->{$inPDB}->{"CharmmElec"}) {
		$row=();
		push(@{$row},"CharmmElec",$self->{intDat}->{$inPDB}->{"CharmmElec"});
		push(@{$table},$row);
	    }
	    &writeHtmlTable($outF,$table);
	}

	# Predicted stabilities of monomers
	if (defined $self->{intDat}->{$inPDB}->{"AvsWT"}) {
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    printf $outF "Predicted stability differences of monomers\n";
	    $table=();
	    $row=();
	    push(@{$row},"AvsWT",$self->{intDat}->{$inPDB}->{"AvsWT"});
	    push(@{$table},$row);
	    $row=();
	    push(@{$row},"AvsIdeal",$self->{intDat}->{$inPDB}->{"AvsIdeal"});
	    push(@{$table},$row);
	    $row=();
	    push(@{$row},"BvsWT",$self->{intDat}->{$inPDB}->{"BvsWT"});
	    push(@{$table},$row);
	    $row=();
	    push(@{$row},"BvsIdeal",$self->{intDat}->{$inPDB}->{"BvsIdeal"});
	    push(@{$table},$row);
	    &writeHtmlTable($outF,$table);
	}

	# List d_*
	if (defined $self->{intDat}->{$inPDB}->{"d_score"}) {
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    printf $outF "Bound - free:\n";
	    $table=();
	    foreach my $k ( sort keys %{$self->{intDat}->{$inPDB}} ) {
		if (($k =~ /^d_/) && (lc($k) !~ /vspdb/)) {
		    $row=();
		    (my $tag=$k)=~s/^d_//;
		    push(@{$row},$tag,$self->{intDat}->{$inPDB}->{$k});
		    push(@{$table},$row);
		}
	    }
	    &writeHtmlTable($outF,$table);
	}

	# List i_*
	if (defined $self->{intDat}->{$inPDB}->{"i_score"}) {
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    printf $outF "Sums over interface residues:\n";
	    $table=();
	    foreach my $k ( sort keys %{$self->{intDat}->{$inPDB}} ) {
		if (($k =~ /^i_/) && ($k !~ /max|min/)) {
		    $row=();
		    (my $tag=$k)=~s/^i_//;
		    push(@{$row},$tag,$self->{intDat}->{$inPDB}->{$k});
		    push(@{$table},$row);
		}
	    }
	    &writeHtmlTable($outF,$table);
	}

	# List c_*
	if (defined $self->{intDat}->{$inPDB}->{"c_score"}) {
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    printf $outF "Sums over interface core residues:\n";
	    $table=();
	    foreach my $k ( sort keys %{$self->{intDat}->{$inPDB}} ) {
		if (($k =~ /^c_/) && ($k !~ /max|min/)) {
		    $row=();
		    (my $tag=$k)=~s/^c_//;
		    push(@{$row},$tag,$self->{intDat}->{$inPDB}->{$k});
		    push(@{$table},$row);
		}
	    }
	    &writeHtmlTable($outF,$table);
	}

	# List a_*
	if (defined $self->{intDat}->{$inPDB}->{"a_score"}) {
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    printf $outF "Averages over interface residues:\n";
	    $table=();
	    foreach my $k ( sort keys %{$self->{intDat}->{$inPDB}} ) {
		if ($k =~ /^a_/) {
		    $row=();
		    (my $tag=$k)=~s/^a_//;
		    push(@{$row},$tag,$self->{intDat}->{$inPDB}->{$k});
		    push(@{$table},$row);
		}
	    }
	    &writeHtmlTable($outF,$table);
	}

	# List b_*
	if (defined $self->{intDat}->{$inPDB}->{"b_score"}) {
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    printf $outF "Averages over interface core residues:\n";
	    $table=();
	    foreach my $k ( sort keys %{$self->{intDat}->{$inPDB}} ) {
		if (($k =~ /^b_/) && ($k !~ /max|min/)) {
		    $row=();
		    (my $tag=$k)=~s/^b_//;
		    push(@{$row},$tag,$self->{intDat}->{$inPDB}->{$k});
		    push(@{$table},$row);
		}
	    }
	    &writeHtmlTable($outF,$table);
	}

	# List max/min
	if (defined $self->{intDat}->{$inPDB}->{"i_maxfa_rep"}) {
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    printf $outF "Interface max/min residue values:\n";
	    $table=();
	    foreach my $k ( sort keys %{$self->{intDat}->{$inPDB}} ) {
		if (($k =~ /^i_/) && ($k =~ /max|min/)) {
		    $row=();
		    (my $tag=$k)=~s/^i_//;
		    push(@{$row},$tag,$self->{intDat}->{$inPDB}->{$k});
		    push(@{$table},$row);
		}
	    }
	    &writeHtmlTable($outF,$table);
	}

	# List max/min
	if (defined $self->{intDat}->{$inPDB}->{"c_maxfa_rep"}) {
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    printf $outF "Interface max/min residue values:\n";
	    $table=();
	    foreach my $k ( sort keys %{$self->{intDat}->{$inPDB}} ) {
		if (($k =~ /^c_/) && ($k =~ /max|min/)) {
		    $row=();
		    (my $tag=$k)=~s/^c_//;
		    push(@{$row},$tag,$self->{intDat}->{$inPDB}->{$k});
		    push(@{$table},$row);
		}
	    }
	    &writeHtmlTable($outF,$table);
	}

	# List the UNS
	my $unsName="outPDBs/".$inPDB;
	$unsName=~s/\.gz//g;
	$unsName=~s/\.pdb//g;
	$unsName.=".unsIntVsRef";
	if (&GenUtil::exists($unsName)) {
	    my $inUns=&GenUtil::getInputFile($unsName);
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    printf $outF "Unsatisfied Hbond donors/acceptors:\n";
	    $table=();
	    while (<$inUns>) {
		chomp;
		my @p=split(' ',$_);
		$row=();
		foreach my $e ( @p ) {
		    push(@{$row},$e);
		}
		push(@{$table},$row);
	    }
	    &writeHtmlTable($outF,$table);
	}

	# List effects of individual mutations
	if (defined $self->{resDat}->{$inPDB}) {
	    printf $outF "\<left\>\n";
	    printf $outF "\<h3\>\n";
	    printf $outF "Effects of individual mutations:\n";
	    my @resfields=();
	    push(@resfields,"ddG");
	    push(@resfields,"dG_vsIdeal_DEScontext");
	    push(@resfields,"dG_vsWT_DEScontext");
	    push(@resfields,"dG_vsIdeal_WTcontext");
	    push(@resfields,"dG_vsWT_WTcontext");
	    push(@resfields,"consensusStr");

	    $table=();
	    $row=();
	    push(@{$row},"Mutation");
	    foreach my $rf ( @resfields ) {
		push(@{$row},$rf);
	    }
	    push(@{$table},$row);

	    # Print the mutations sorted by resnum
	    foreach my $k ( sort { $self->{resDat}->{$inPDB}->{$a}->{"resnum"} <=> 
				       $self->{resDat}->{$inPDB}->{$b}->{"resnum"} }
			    keys %{ $self->{resDat}->{$inPDB}} ) {
		$row=();
		push(@{$row},$k);
		foreach my $rf ( @resfields ) {
		    if ($rf eq "consensusStr") {
			my $str=$self->{resDat}->{$inPDB}->{$k}->{"consensusStr"};
			my @aa_arr=split(/[0-9]+/,$str);
			my @ratio_arr=split(/[A-Z]+/,$str);
			shift(@ratio_arr);
			my $harr=();
			for (my $i=0; $i<=$#aa_arr; $i++) {
			    my $f={};
			    $f->{ratio}=$ratio_arr[$i];
			    $f->{aa}=$aa_arr[$i];
			    push(@{$harr},$f);
			}
			my $clean_arr=();
			foreach my $k ( sort { $b->{ratio} <=> $a->{ratio} } @{$harr} ) {
			    push(@{$clean_arr},($k->{aa}.":".$k->{ratio}));
			}
			push(@{$row},join(", ",@{$clean_arr}));
		    } else {
			push(@{$row},$self->{resDat}->{$inPDB}->{$k}->{$rf});
		    }
		}
		push(@{$table},$row);
	    }
	    &writeHtmlTable($outF,$table);
	}

	printf $outF "\<\/body\>\n\n";
	undef $outF;
    }


    return;
}


## function: writeHtmlTable(outF,table_data)
## write tables in assessment html pages

sub writeHtmlTable {

    my $outF=shift;
    my $table=shift;

    printf $outF "\<table border=2\>\n";
    foreach my $row ( @{$table} ) {
	printf $outF "\<tr\>";
	foreach my $ele ( @{$row} ) {
	    printf $outF "\<td\>";
	    printf $outF "\<center\>";
	    if (((length($ele) < 7)) || ($ele =~ /[A-Za-z]/)) {
		printf $outF " %s ",$ele;
	    } else {
		printf $outF " %.2f ",$ele;
	    }
	}
	printf $outF "\n";
    }
    printf $outF "\<\/table\>\n";

    return;
}


## method: _print_header(outdat,field_list)
## print header info

sub _print_header {

    my $self=shift;
    my $outdat=shift;
    my $field_list=shift;
    my $by_res=shift;

    $by_res=0 if (! defined $by_res);

    printf $outdat "fname ";
    printf $outdat "mut " if ($by_res);
    printf $outdat "%s\n",join(" ",@{$field_list});

    return;
}

## function: _setup_Ros_fields()
## expected fields when calling rosetta
sub _setup_Ros_fields {

    my $field_list=();

    # Measures of interface size
    push(@{$field_list},"delta_sasa");
    push(@{$field_list},"num_interface_res");
    push(@{$field_list},"num_interface_atoms");
    push(@{$field_list},"num_core_interface_res");
    push(@{$field_list},"ncont");

    # Hotspot-related
    push(@{$field_list},"num_hotspot_residues");
    push(@{$field_list},"best_hotspot_atr_ene");
    push(@{$field_list},"best_hotspot_SASApack");
    push(@{$field_list},"best_hotspot_SASAprob");

    # Stack-related
    push(@{$field_list},"num_stacks");
    push(@{$field_list},"best_stack_intrastack_ene");
    push(@{$field_list},"best_stack_plane_ene");
    push(@{$field_list},"best_stack_atr_ene");
    push(@{$field_list},"best_stack_SASApack");
    push(@{$field_list},"best_stack_SASAprob");

    # Measures of composition
    push(@{$field_list},"interface_num_aromatic_atoms");
    push(@{$field_list},"interface_aromatic_sasa");
    push(@{$field_list},"interface_fraction_polar");
    push(@{$field_list},"met_sc_atoms");

    # Measures of packing
    push(@{$field_list},"interface_num_tight_core_atoms");
    push(@{$field_list},"num_tight_core_atoms_per_sasa");
    push(@{$field_list},"interface_fraction_tight");
    push(@{$field_list},"delta_atr_per_sasa");

    # Bound-free energetic measures
    push(@{$field_list},"delta_atr");
    push(@{$field_list},"delta_rep");
    push(@{$field_list},"delta_sol");
    push(@{$field_list},"delta_hbond");
    push(@{$field_list},"delta_pair");
#    push(@{$field_list},"delta_gb");
#    push(@{$field_list},"delta_cst");
    push(@{$field_list},"delta_plane");
#    push(@{$field_list},"delta_h2o");
#    push(@{$field_list},"delta_h2o_hb");
    push(@{$field_list},"delta_score");

    # Residue energetic measures
    push(@{$field_list},"interface_atr");
    push(@{$field_list},"interface_rep");
    push(@{$field_list},"interface_sol");
    push(@{$field_list},"interface_hbond");
    push(@{$field_list},"interface_prob");
    push(@{$field_list},"interface_dun");
    push(@{$field_list},"interface_pair");
    push(@{$field_list},"interface_intrares");
    push(@{$field_list},"interface_total");

#    push(@{$field_list},"interface_gb");
#    push(@{$field_list},"interface_cst");
#    push(@{$field_list},"interface_h2o");
#    push(@{$field_list},"interface_h2o_hb");
#    push(@{$field_list},"interface_wsol");
#    push(@{$field_list},"interface_unf");

    push(@{$field_list},"interface_atr_VsPDB");
    push(@{$field_list},"interface_rep_VsPDB");
    push(@{$field_list},"interface_sol_VsPDB");
    push(@{$field_list},"interface_hbond_VsPDB");
    push(@{$field_list},"interface_prob_VsPDB");
    push(@{$field_list},"interface_dun_VsPDB");
    push(@{$field_list},"interface_pair_VsPDB");
    push(@{$field_list},"interface_intrares_VsPDB");
    push(@{$field_list},"interface_total_VsPDB");

    push(@{$field_list},"interface_SASApackVsPDB");
    push(@{$field_list},"interface_SASAprobVsPDB");

    push(@{$field_list},"interface_max_rep");
    push(@{$field_list},"interface_max_atr_VsPDB");
    push(@{$field_list},"interface_max_SASApackVsPDB");
    push(@{$field_list},"interface_min_SASAprobVsPDB");

    # Residue energetic measures averaged over interface residues
    push(@{$field_list},"interface_ave_atr");
    push(@{$field_list},"interface_ave_rep");
    push(@{$field_list},"interface_ave_sol");
    push(@{$field_list},"interface_ave_hbond");
    push(@{$field_list},"interface_ave_prob");
    push(@{$field_list},"interface_ave_dun");
    push(@{$field_list},"interface_ave_pair");
    push(@{$field_list},"interface_ave_intrares");
    push(@{$field_list},"interface_ave_total");

#    push(@{$field_list},"interface_ave_gb");
#    push(@{$field_list},"interface_ave_cst");
#    push(@{$field_list},"interface_ave_h2o");
#    push(@{$field_list},"interface_ave_h2o_hb");
#    push(@{$field_list},"interface_ave_wsol");
#    push(@{$field_list},"interface_ave_unf");

    push(@{$field_list},"interface_ave_atr_VsPDB");
    push(@{$field_list},"interface_ave_rep_VsPDB");
    push(@{$field_list},"interface_ave_sol_VsPDB");
    push(@{$field_list},"interface_ave_hbond_VsPDB");
    push(@{$field_list},"interface_ave_prob_VsPDB");
    push(@{$field_list},"interface_ave_dun_VsPDB");
    push(@{$field_list},"interface_ave_pair_VsPDB");
    push(@{$field_list},"interface_ave_intrares_VsPDB");
    push(@{$field_list},"interface_ave_total_VsPDB");

    push(@{$field_list},"interface_ave_SASApackVsPDB");
    push(@{$field_list},"interface_ave_SASAprobVsPDB");

    # Core residue energetic measures
    push(@{$field_list},"core_interface_atr");
    push(@{$field_list},"core_interface_rep");
    push(@{$field_list},"core_interface_sol");
    push(@{$field_list},"core_interface_hbond");
    push(@{$field_list},"core_interface_prob");
    push(@{$field_list},"core_interface_dun");
    push(@{$field_list},"core_interface_pair");
    push(@{$field_list},"core_interface_intrares");
    push(@{$field_list},"core_interface_total");

#    push(@{$field_list},"core_interface_gb");
#    push(@{$field_list},"core_interface_cst");
#    push(@{$field_list},"core_interface_h2o");
#    push(@{$field_list},"core_interface_h2o_hb");
#    push(@{$field_list},"core_interface_wsol");
#    push(@{$field_list},"core_interface_unf");

    push(@{$field_list},"core_interface_atr_VsPDB");
    push(@{$field_list},"core_interface_rep_VsPDB");
    push(@{$field_list},"core_interface_sol_VsPDB");
    push(@{$field_list},"core_interface_hbond_VsPDB");
    push(@{$field_list},"core_interface_prob_VsPDB");
    push(@{$field_list},"core_interface_dun_VsPDB");
    push(@{$field_list},"core_interface_pair_VsPDB");
    push(@{$field_list},"core_interface_intrares_VsPDB");
    push(@{$field_list},"core_interface_total_VsPDB");

    push(@{$field_list},"core_interface_SASApackVsPDB");
    push(@{$field_list},"core_interface_SASAprobVsPDB");

    push(@{$field_list},"core_interface_max_rep");
    push(@{$field_list},"core_interface_max_atr_VsPDB");
    push(@{$field_list},"core_interface_max_SASApackVsPDB");
    push(@{$field_list},"core_interface_min_SASAprobVsPDB");

    # Residue energetic measures averaged over interface residues
    push(@{$field_list},"core_interface_ave_atr");
    push(@{$field_list},"core_interface_ave_rep");
    push(@{$field_list},"core_interface_ave_sol");
    push(@{$field_list},"core_interface_ave_hbond");
    push(@{$field_list},"core_interface_ave_prob");
    push(@{$field_list},"core_interface_ave_dun");
    push(@{$field_list},"core_interface_ave_pair");
    push(@{$field_list},"core_interface_ave_intrares");
    push(@{$field_list},"core_interface_ave_total");

#    push(@{$field_list},"core_interface_ave_gb");
#    push(@{$field_list},"core_interface_ave_cst");
#    push(@{$field_list},"core_interface_ave_h2o");
#    push(@{$field_list},"core_interface_ave_h2o_hb");
#    push(@{$field_list},"core_interface_ave_wsol");
#    push(@{$field_list},"core_interface_ave_unf");

    push(@{$field_list},"core_interface_ave_atr_VsPDB");
    push(@{$field_list},"core_interface_ave_rep_VsPDB");
    push(@{$field_list},"core_interface_ave_sol_VsPDB");
    push(@{$field_list},"core_interface_ave_hbond_VsPDB");
    push(@{$field_list},"core_interface_ave_prob_VsPDB");
    push(@{$field_list},"core_interface_ave_dun_VsPDB");
    push(@{$field_list},"core_interface_ave_pair_VsPDB");
    push(@{$field_list},"core_interface_ave_intrares_VsPDB");
    push(@{$field_list},"core_interface_ave_total_VsPDB");

    push(@{$field_list},"core_interface_ave_SASApackVsPDB");
    push(@{$field_list},"core_interface_ave_SASAprobVsPDB");

    return $field_list;
}


## function: _setup_Wat_fields()
## expected fields when calling rosetta with explicit water
## note: these are in addition to the regular interface_ds fields
sub _setup_Wat_fields {

    my $field_list = &_setup_Ros_fields();

    # Put these fields at the front, so that "first_field" will
    # report missing data (ie. UNS data) if UNS data is required
    unshift(@{$field_list},"INT_bbhb_uns");
    unshift(@{$field_list},"INT_schb_uns");
    unshift(@{$field_list},"INT_uns_score");
    unshift(@{$field_list},"INT_group_uns");
    unshift(@{$field_list},"INT_group_uns_score");

    return $field_list;
}


## method: _call_mutlist()
## find the mutations relative to WT
sub _call_mutlist {
    my $self=shift;
    my $missing_inputs=shift;
    my $outdat=shift;
    my $print_header=shift;

    my $wtStruct=$self->{files}->{wtStruct};

    my $field_list=();
    push(@{$field_list},"Mutlist");
    push(@{$field_list},"NumMut");
    push(@{$field_list},"MutA");
    push(@{$field_list},"NumMutA");
    push(@{$field_list},"MutB");
    push(@{$field_list},"NumMutB");
    $self->_print_header($outdat,$field_list)
	if ((defined $print_header) && ($print_header));

    my $WTmol=Molecule->new($wtStruct);

    foreach my $inPDB (@{$missing_inputs}) {

	my $mol=Molecule->new($self->{cmap}->{$inPDB});
	my $mutlist=$mol->generateMutlist($WTmol);
	my $numMut=$#{$mutlist}+1;

	my $mutStr=&Molecule::generateMutStr($mutlist);
	$self->{intDat}->{$inPDB}->{"Mutlist"}=$mutStr;
	$self->{intDat}->{$inPDB}->{"NumMut"}=$numMut;

	my $Acid=$mol->activeChains()->[0]->{id};
	my $mutStrA=&Molecule::generateMutStr($mutlist,$Acid);
	my @Alist=split(/,/,$mutStrA);
	my $numMutA=$#Alist+1;
	$self->{intDat}->{$inPDB}->{"MutA"}=$mutStrA;
	$self->{intDat}->{$inPDB}->{"NumMutA"}=$numMutA;

	my $Bcid=$mol->activeChains()->[1]->{id};
	my $mutStrB=&Molecule::generateMutStr($mutlist,$Bcid);
	my @Blist=split(/,/,$mutStrB);
	my $numMutB=$#Blist+1;
	$self->{intDat}->{$inPDB}->{"MutB"}=$mutStrB;
	$self->{intDat}->{$inPDB}->{"NumMutB"}=$numMutB;

	printf $outdat "%s %s %d %s %s %s %s\n",$inPDB,$mutStr,
	    $numMut,$mutStrA,$numMutA,$mutStrB,$numMutB;
    }

    return;
}


## method: _call_Cha()
## call Charmm for assessment
sub _call_Cha {
    my $self=shift;
    my $missing_inputs=shift;
    my $outdat=shift;
    my $print_header=shift;

    my $field="CharmmElec";

    my $field_list=();
    push(@{$field_list},$field);
    $self->_print_header($outdat,$field_list)
	if ((defined $print_header) && ($print_header));

    my $charmm;
    my $logfile="cha.log";
    my %par;
#	&GenUtil::parsePar(\%par,"param=19");
    &GenUtil::parsePar(\%par,"param=22x");
    &GenUtil::parsePar(\%par,"cutnb=20.0");
    &GenUtil::parsePar(\%par,"cutoff=18.0");
    &GenUtil::parsePar(\%par,"cuton=16.0");
    &GenUtil::parsePar(\%par,"shake=1");
#	&GenUtil::parsePar(\%par,"gb=bgb");
    &GenUtil::parsePar(\%par,"gb=gbmva");
    &GenUtil::parsePar(\%par,"gbmvsa=0.015");
    &GenUtil::parsePar(\%par,"gbmvas=0.9");
    &GenUtil::parsePar(\%par,"gbmvap7=0.0001");
    &GenUtil::parsePar(\%par,"gbmvad=-0.1");
    &GenUtil::parsePar(\%par,"gbmvaig=32");
    &GenUtil::parsePar(\%par,"scalerad=mnina");
    &GenUtil::parsePar(\%par,"cmap");
    &GenUtil::parsePar(\%par,"xpar=par_cmap27.inp");
    &GenUtil::parsePar(\%par,"xtop=top_all22_prot_cmap.inp");
    &GenUtil::parsePar(\%par,"gbmvafrq=1");
    &GenUtil::parsePar(\%par,"gbmvaemp=999999.0");

    $charmm=&CHARMM::new($logfile);
    $charmm->loadParameters(%par);

    foreach my $inPDB (@{$missing_inputs}) {

	# Setup the bound structure
	my $mol=Molecule->new();
	$mol->readPDB($self->{cmap}->{$inPDB},ignorehet=>1);
	$mol->fixHistidine($self->{par}->{hsd},
		    $self->{par}->{hse},$self->{par}->{hsp});
	$mol->changeResName($self->{par}->{resmod});
	$mol->generateSegNames();
	$charmm->setupFromMolecule($mol,undef);

	# Get the energy of the complex
	$charmm->setupEnergy();
	my $cener=$charmm->getEnergy();

	my $chaCmd;

	# Find the center of mass of first protein
	my $Acid=$mol->{chain}->[0]->{id};
	$chaCmd = "coor stat mass sele iseg 1 end";
	$charmm->verbose($chaCmd);
	$chaCmd = "set xdir ?xave";
	$charmm->verbose($chaCmd);
	$chaCmd = "set ydir ?yave";
	$charmm->verbose($chaCmd);
	$chaCmd = "set zdir ?zave";
	$charmm->verbose($chaCmd);

	# Find the center of mass of second protein
	my $Bcid=$mol->{chain}->[1]->{id};
	$chaCmd = "coor stat mass sele iseg 2 end";
	$charmm->verbose($chaCmd);
	$chaCmd = "decr xdir by ?xave";
	$charmm->verbose($chaCmd);
	$chaCmd = "decr ydir by ?yave";
	$charmm->verbose($chaCmd);
	$chaCmd = "decr zdir by ?zave";
	$charmm->verbose($chaCmd);

	# Move the complex apart
	my $movedist=50.0;
	$chaCmd = "coor trans xdir \@xdir ydir \@ydir zdir \@zdir";
	$chaCmd .= " fact 4";
	$chaCmd .= " sele iseg 1 end";
	$charmm->verbose($chaCmd);

	# Get the energy of the unbound
	my $fener=$charmm->getEnergy();
	$charmm->clearEnergy();
	$charmm->verbose("dele atom sele all end");

	# Find the difference upon binding
	my $E=($cener->{elec}+$cener->{gb})-($fener->{elec}+$fener->{gb});
	$self->{intDat}->{$inPDB}->{$field}=$E;
	printf $outdat "%s %f\n",$inPDB,$E;

    }

    $charmm->finish();
    &GenUtil::remove($logfile) if (! $debug);

    return;
}


## method: _call_Ros()
## call Rosetta, using interface_ds
sub _call_Ros {
    my $self=shift;
    my $missing_inputs=shift;
    my $outdat=shift;
    my $print_header=shift;

    my $field_list = &_setup_Ros_fields();
    $self->_print_header($outdat,$field_list)
	if ((defined $print_header) && ($print_header));

    my $sArr=();
    my $outlookup={};
    foreach my $fname (@{$missing_inputs}) {
	my $c={};
	$c->{filename}=$self->{cmap}->{$fname};
	push(@{$sArr},$c);
	$outlookup->{$self->{cmap}->{$fname}}=$fname;
    }

    # Call Rosetta
    my $outdir="tmp_Ros_output";
# Don't use rosetta until interface_ds has been incorporated
#    &RosettaMolecule::useSvnRosettaExec();
    my $Rmol=RosettaMolecule->new($sArr->[0]->{filename});
    $Rmol->setParameter(pathfile=>"tmp_paths.txt");
    $Rmol->setParameter(outdir=>$outdir);
    $Rmol->setParameter(stdout=>"tmp_rosetta.out");
    $Rmol->setParameter(stderr=>"tmp_rosetta.err");
    $Rmol->setParameter(multi_chain=>1);
    $Rmol->setParameter(decoystats=>1);
    $Rmol->setParameter(interface_ds=>1);
    $Rmol->setParameter(tight_core=>1);
    $sArr=$Rmol->computeScore($sArr);

    # Parse output
    my $outfiles={};
    foreach my $f ( @{$sArr} ) {
	my $inpname=$f->{origfile};
	my $fname=$outlookup->{$inpname};
	my $type="C";
	$outfiles->{$fname}->{$type}=$f->{filename};
    }

    foreach my $fname ( keys %{$outfiles} ) {

	# Read in output PDB files, including details of scores
	my $Cout=$outfiles->{$fname}->{C};
	my $Cmol=RosettaMolecule->new($Cout);
	$Cmol->updateScore($Cout,1,1);

	printf $outdat "%s",$fname;
	foreach my $field ( @{$field_list} ) {
	    my $val=$Cmol->{"score"}->{$field};
	    if (! defined $val) {
		$val = 9999.9;
		$val = -9999.9 if ($field =~ "SASAprobVsPDB");
	    }
	    $self->{intDat}->{$fname}->{$field}=$val;
	    printf $outdat " %s",$val;
	}
	printf $outdat "\n";

    }

    if (! $debug) {
	# Clean up
	&GenUtil::remove($outdir);
    }

    return;
}


## method: _call_Wat()
## call Rosetta with explicit water
sub _call_Wat {
    my $self=shift;
    my $missing_inputs=shift;
    my $outdat=shift;
    my $print_header=shift;

    my $field_list = &_setup_Wat_fields();
    $self->_print_header($outdat,$field_list)
	if ((defined $print_header) && ($print_header));

    my $sArr=();
    my $outlookup={};
    foreach my $fname (@{$missing_inputs}) {
	my $r={};
	$r->{filename}=$self->{cmap}->{$fname};
	$outlookup->{$self->{cmap}->{$fname}}=$fname;
	push(@{$sArr},$r);
    }

    # Call Rosetta
    my $outdir="W_tmp_Ros_output";
    my $Rmol=RosettaMolecule->new($sArr->[0]->{filename});
    $Rmol->setParameter(pathfile=>"W_tmp_paths.txt");
    $Rmol->setParameter(outdir=>$outdir);
    $Rmol->setParameter(stdout=>"W_tmp_rosetta.out");
    $Rmol->setParameter(stderr=>"W_tmp_rosetta.err");
    $sArr=$Rmol->identifyInterfaceUNS($sArr);

    # For PDBs colored by LJatr versus avePDB and lists of UNS polars
    &GenUtil::makeDir($outPDBdir);

    # Parse output
    my $outfiles={};
    foreach my $f ( @{$sArr} ) {
	my $inpname=$f->{origfile};
	my $fname=$outlookup->{$inpname};
	my $type="C";
	$outfiles->{$fname}->{$type}=$f->{filename};
    }

    foreach my $fname ( keys %{$outfiles} ) {

	# Read in output PDB files, including details of scores
	my $Cout=$outfiles->{$fname}->{C};
	my $Cmol=RosettaMolecule->new($Cout);
	$Cmol->updateScore($Cout,1,1);
	&GenUtil::link($Cout,$outPDBdir.$fname);

	# Write a list of interface UNS to a separate file
	my $outUnsFname=$outPDBdir.$fname.".uns";
	$Cmol->writeUns($outUnsFname);

	printf $outdat "%s",$fname;
	foreach my $field ( @{$field_list} ) {
	    my $val=$Cmol->{"score"}->{$field};
	    if (! defined $val) {
		$val = 9999.9;
		$val = -9999.9 if ($field =~ "SASAprobVsPDB");
	    }
	    $self->{intDat}->{$fname}->{$field}=$val;
	    printf $outdat " %s",$val;
	}
	printf $outdat "\n";

    }

    if (! $debug) {
	# Clean up
	&GenUtil::remove($outdir);
    }

    return;
}


## method: _call_Mon()
## call Rosetta to assess the monomer stabilities
sub _call_Mon {
    my $self=shift;
    my $missing_inputs=shift;
    my $outdat=shift;
    my $print_header=shift;

    my $wtStruct=$self->{files}->{wtStruct};

    my $field_list=();
    push(@{$field_list},"AvsIdeal");
    push(@{$field_list},"AvsWT");
    push(@{$field_list},"BvsIdeal");
    push(@{$field_list},"BvsWT");

    $self->_print_header($outdat,$field_list)
	if ((defined $print_header) && ($print_header));

    die "Cannot find WT struct" if (! &GenUtil::exists($wtStruct));

    # Note: we'll do calculations using the WT structure, we'll
    # only use the designed structure to get list of mutations
    my $WTmol=Molecule->new($wtStruct);

    foreach my $inPDB (@{$missing_inputs}) {

	my $mutStr=$self->{intDat}->{$inPDB}->{"Mutlist"};
	my $mutlist = &Molecule::mutStrToMutlist($mutStr);

	my $Acid=$WTmol->activeChains()->[0]->{id};
	my $Bcid=$WTmol->activeChains()->[1]->{id};

	my $Aideal = $self->_computeEne($wtStruct,$Acid,"ideal",$mutlist);
	my $Ades = $self->_computeEne($wtStruct,$Acid,"des",$mutlist);
	my $Awt = $self->_computeEne($wtStruct,$Acid,"wt",$mutlist);

	my $Bideal = $self->_computeEne($wtStruct,$Bcid,"ideal",$mutlist);
	my $Bdes = $self->_computeEne($wtStruct,$Bcid,"des",$mutlist);
	my $Bwt = $self->_computeEne($wtStruct,$Bcid,"wt",$mutlist);

	my $Ai=$Ades-$Aideal;
	my $Aw=$Ades-$Awt;
	my $Bi=$Bdes-$Bideal;
	my $Bw=$Bdes-$Bwt;

	$self->{intDat}->{$inPDB}->{"AvsIdeal"}=$Ai;
	$self->{intDat}->{$inPDB}->{"AvsWT"}=$Aw;
	$self->{intDat}->{$inPDB}->{"BvsIdeal"}=$Bi;
	$self->{intDat}->{$inPDB}->{"BvsWT"}=$Bw;

	printf $outdat "%s %f %f %f %f\n",$inPDB,$Ai,$Aw,$Bi,$Bw;
    }

    return;
}


## method: _call_ResMinMax()
## call to find min/max of residue energetics
sub _call_ResMinMax {
    my $self=shift;
    my $missing_inputs=shift;
    my $outdat=shift;
    my $print_header=shift;

    my $inPDB=$missing_inputs->[0];
    my $mutStr=$self->{intDat}->{$inPDB}->{"Mutlist"};
    my @mutlist=split(",",$mutStr);
    my $mut=$mutlist[0];
    my $base_field_list=();
    my $field_list=();
    foreach my $k ( keys %{$self->{resDat}->{$inPDB}->{$mut}} ) {
	push(@{$base_field_list},$k);
	push(@{$field_list},"i_min_".$k);
	push(@{$field_list},"i_max_".$k);
    }

    $self->_print_header($outdat,$field_list)
	if ((defined $print_header) && ($print_header));

    foreach my $inPDB ( @{$missing_inputs} ) {
	my $mutStr=$self->{intDat}->{$inPDB}->{"Mutlist"};
	my @mutlist=split(",",$mutStr);
	foreach my $mut ( @mutlist ) {
	    foreach my $base_field ( @{$base_field_list} ) {
		my $minName="i_min_".$base_field;
		my $maxName="i_max_".$base_field;
		my $currval=$self->{resDat}->{$inPDB}->{$mut}->{$base_field};
		if ((! defined $self->{intDat}->{$inPDB}->{$minName}) ||
		    ($currval < $self->{intDat}->{$inPDB}->{$minName})) {
		    $self->{intDat}->{$inPDB}->{$minName} = $currval;
		}
		if ((! defined $self->{intDat}->{$inPDB}->{$maxName}) ||
		    ($currval > $self->{intDat}->{$inPDB}->{$maxName})) {
		    $self->{intDat}->{$inPDB}->{$maxName} = $currval;
		}
	    }
	}

	printf $outdat "%s",$inPDB;
	foreach my $field ( @{$field_list} ) {
	    printf $outdat " %s",$self->{intDat}->{$inPDB}->{$field};
	}
	printf $outdat "\n";

    }

    return;
}


## method: _call_ddG
## calculate ddG info
sub _call_ddG {
    my $self=shift;
    my $missing_inputs=shift;
    my $outdat=shift;
    my $print_header=shift;

    # We can look this up if it's already computed...
    printf STDERR "Computing free differences\n";
    $self->collectInterfaceDataType("DEScontext.out",1,\&{_call_free_dG});

    printf STDERR "Computing bound differences\n";

    my $field_list=();
    push(@{$field_list},"ddG");

    $self->_print_header($outdat,$field_list,1)
	if ((defined $print_header) && ($print_header));

    foreach my $inMut (@{$missing_inputs}) {

	my @p=split("::",$inMut);
	my $inPDB=$p[0];
	my $mut=$p[1];
	my $resnum=substr($mut,1,-1);
	my $mutlist=&Molecule::mutStrToMutlist($mut);

	my $des = $self->_computeEne($self->{cmap}->{$inPDB},undef,"des",$mutlist);
	my $wt = $self->_computeEne($self->{cmap}->{$inPDB},undef,"wt",$mutlist);

	my $CdG=$des-$wt;
	# We've already computed dG for the free elsewhere
	my $FdG=$self->{resDat}->{$inPDB}->{$mut}->{"dG_vsWT_DEScontext"};
	my $ddG=$CdG-$FdG;

	$self->{resDat}->{$inPDB}->{$mut}->{"ddG"}=$ddG;

	printf $outdat "%s %s %f\n",$inPDB,$mut,$ddG;

    }

    return;
}


## method: _call_free_dG
## calculate  info
sub _call_free_dG {
    my $self=shift;
    my $missing_inputs=shift;
    my $outdat=shift;
    my $print_header=shift;

    my $field_list=();
    push(@{$field_list},"dG_vsIdeal_DEScontext");
    push(@{$field_list},"dG_vsWT_DEScontext");

    $self->_print_header($outdat,$field_list,1)
	if ((defined $print_header) && ($print_header));

    foreach my $inMut (@{$missing_inputs}) {

	my @p=split("::",$inMut);
	my $inPDB=$p[0];
	my $mut=$p[1];
	my $resnum=substr($mut,1,-1);
	my $mutlist=&Molecule::mutStrToMutlist($mut);

	# Write the unbound PDB
	my $FreeFname="tmp_free.pdb";
	my $mol=RosettaMolecule->new($self->{cmap}->{$inPDB});
	my $Acid=$mol->{chain}->[0]->{id};
	my $Bcid=$mol->{chain}->[1]->{id};
	my ($Ax,$Ay,$Az)=$mol->centerOfMass($Acid);
	my ($Bx,$By,$Bz)=$mol->centerOfMass($Bcid);
	my $xdiff=$Ax-$Bx;
	my $ydiff=$Ay-$By;
	my $zdiff=$Az-$Bz;
	my $mag=$xdiff*$xdiff;
	$mag+=$ydiff*$ydiff;
	$mag+=$zdiff*$zdiff;
	$mag=sqrt($mag);
	my $movedist=50.0;
	$xdiff*=$movedist/$mag;
	$ydiff*=$movedist/$mag;
	$zdiff*=$movedist/$mag;
	$mol->move($xdiff,$ydiff,$zdiff,$Acid);
	$mol->writePDB($FreeFname,"coor");

	my $ideal=$self->_computeEne($FreeFname,undef,"ideal",$mutlist);
	my $des=$self->_computeEne($FreeFname,undef,"des",$mutlist);
	my $wt=$self->_computeEne($FreeFname,undef,"wt",$mutlist);

	my $i=$des-$ideal;
	my $w=$des-$wt;

	$self->{resDat}->{$inPDB}->{$mut}->{"dG_vsIdeal_DEScontext"}=$i;
	$self->{resDat}->{$inPDB}->{$mut}->{"dG_vsWT_DEScontext"}=$w;

	printf $outdat "%s %s %f %f\n",$inPDB,$mut,$i,$w;

    }

    return;
}


## method: _call_WT_dG
## calculate ddG info
sub _call_WT_dG {
    my $self=shift;
    my $missing_inputs=shift;
    my $outdat=shift;
    my $print_header=shift;

    my $wtStruct=$self->{files}->{wtStruct};

    my $field_list=();
    push(@{$field_list},"dG_vsIdeal_WTcontext");
    push(@{$field_list},"dG_vsWT_WTcontext");

    $self->_print_header($outdat,$field_list,1)
	if ((defined $print_header) && ($print_header));

    foreach my $inMut (@{$missing_inputs}) {

	my @p=split("::",$inMut);
	my $inPDB=$p[0];
	my $mut=$p[1];
	my $resnum=substr($mut,1,-1);
	my $mutlist=&Molecule::mutStrToMutlist($mut);

	my $ideal = $self->_computeEne($wtStruct,undef,"ideal",$mutlist);
	my $des = $self->_computeEne($wtStruct,undef,"des",$mutlist);
	my $wt = $self->_computeEne($wtStruct,undef,"wt",$mutlist);

	my $i=$des-$ideal;
	my $w=$des-$wt;

	$self->{resDat}->{$inPDB}->{$mut}->{"dG_vsIdeal_WTcontext"}=$i;
	$self->{resDat}->{$inPDB}->{$mut}->{"dG_vsWT_WTcontext"}=$w;

	printf $outdat "%s %s %f %f\n",$inPDB,$mut,$i,$w;
    }

    return;
}


## method: _call_consensus_seq()
## call to find SVN rosetta consensus residues
sub _call_consensus_seq {
    my $self=shift;
    my $missing_inputs=shift;
    my $outdat=shift;
    my $print_header=shift;

    my $field_list=();
    push(@{$field_list},"consensusStr");

    $self->_print_header($outdat,$field_list,1)
	if ((defined $print_header) && ($print_header));

    my $tot_counts=100;
    my $donePDBs={};

    # Call rosetta only once for all the mutations in a given PDB,
    # but store the results in the "per mutation" format
    foreach my $inMut (@{$missing_inputs}) {
	my @p=split("::",$inMut);
	my $inPDB=$p[0];
	if (! defined $donePDBs->{$inPDB}) {

	    my $pos_counts={};
	    my $mutStr=$self->{intDat}->{$inPDB}->{"Mutlist"};
	    my $mutlist = &Molecule::mutStrToMutlist($mutStr);

	    my $Nmol=RosettaMolecule->new($self->{cmap}->{$inPDB});

	    $Nmol->setAllDesType(undef,"DEFLT");
	    my $resfile=$self->{files}->{resfile};
	    $Nmol->readResfile($resfile) if (defined $resfile);
	    $Nmol->changeDesType("DEFLT","NATAA");

	    foreach my $mut ( @{$mutlist} ) {
	      foreach my $c ( @{$Nmol->activeChains()} ) {
		  my $prevtype = $Nmol->reportResDesType($c->{id},$mut->{resnum});
		  if (defined $prevtype) {
		    if ($prevtype eq "NATAA") {
		      $Nmol->setResDesType($c->{id},$mut->{resnum},"ALLAA");
		    } else {
		      my @p=split(' ',$prevtype);
		      if ($#p == 1) {
			die "Unexpected resfile entry" if ($p[0] ne "NOTAA");
			my $excl=$p[1];
			$Nmol->setResDesType($c->{id},$mut->{resnum},
				     "PIKAA  ".&Sequence::notRes($excl));
		      }
		    }
		  }
	      }
	    }

	    $Nmol->changeDesType("NOTAA","NATAA");

	    my $outdir="tmp_consensus_dir";
	    $Nmol->setParameter(outdir=>$outdir);
	    $Nmol->setParameter(nstruct=>$tot_counts);
	    $Nmol->setParameter(ex_rot=>2);
	    $Nmol->setParameter(use_input_sc=>0);

	    &RosettaMolecule::useSvnRosettaExec();
	    my $opt_out=$Nmol->design();

	    foreach my $mut ( @{$mutlist} ) {
		my $resnum=$mut->{resnum};
		$pos_counts->{$resnum}={};
		foreach my $aa ( keys %Sequence::_seqlong ) {
		    $pos_counts->{$resnum}->{$aa}=0;
		}
	    }

	    foreach my $outstruct ( @{$opt_out} ) {
		my $mol=Molecule->new($outstruct->{filename});
		foreach my $mut ( @{$mutlist} ) {
		    my $resnum=$mut->{resnum};
		    foreach my $c ( @{$mol->activeChains()} ) {
			my $Mres=$mol->getResidueInChain($resnum,$c);
			if (defined $Mres) {
			    my $aa=$Sequence::_seqabbrev{$Mres->{name}};
			    $pos_counts->{$resnum}->{$aa}++;
			}
		    }
		}
	    }

	    foreach my $m ( @{$mutlist} ) {
		my $cstr="";
		my $resnum=$m->{resnum};
		my $mut=$m->{natoneres}.$resnum.$m->{mutoneres};
		printf $outdat "%s %s ",$inPDB,$mut;
		foreach my $aa ( keys %Sequence::_seqlong ) {
		    my $counts=$pos_counts->{$resnum}->{$aa};
		    if ($counts > 0) {
			printf $outdat "%s%d",$aa,$counts;
			$cstr.=$aa.$counts;
		    }
		}
		printf $outdat "\n";
		$self->{resDat}->{$inPDB}->{$mut}->{"consensusStr"}=$cstr;
	    }

	    if ($debug) {
		my $newdirname="tmp_consensus_dir";
		$newdirname.="_".$inPDB;
		&GenUtil::rename($outdir,$newdirname);
	    } else {
		# Clean up
		&GenUtil::remove($outdir);
	    }
	    $donePDBs->{$inPDB}=1;
	}
    }

    return;
}


## function: _computeEne()
## call Rosetta to assess the monomer stabilities
sub _computeEne {
    my $self=shift;
    my $inputStruct=shift;
    my $cid=shift;
    my $mode=shift;
    my $mutlist=shift;

    my $Nmol=RosettaMolecule->new($inputStruct);

    $Nmol->setAllDesType(undef,"DEFLT");
    my $resfile=$self->{files}->{resfile};
    $Nmol->readResfile($resfile) if (defined $resfile);
    $Nmol->changeDesType("DEFLT","NATAA");

    if (lc($mode) eq "wt") {
	foreach my $mut ( @{$mutlist} ) {
	    foreach my $c ( @{$Nmol->activeChains()} ) {
		$Nmol->setResDesType($c->{id},$mut->{resnum},("PIKAA  ".$mut->{natoneres}));
	    }
	}

    } elsif (lc($mode) eq "des") {
	foreach my $mut ( @{$mutlist} ) {
	    foreach my $c ( @{$Nmol->activeChains()} ) {
		$Nmol->setResDesType($c->{id},$mut->{resnum},("PIKAA  ".$mut->{mutoneres}));
	    }
	}

    } else {

	foreach my $mut ( @{$mutlist} ) {
	    foreach my $c ( @{$Nmol->activeChains()} ) {
		my $prevtype = $Nmol->reportResDesType($c->{id},$mut->{resnum});
		if (defined $prevtype) {
		    if ($prevtype eq "NATAA") {
			$Nmol->setResDesType($c->{id},$mut->{resnum},"ALLAA");
		    } else {
			my @p=split(' ',$prevtype);
			if ($#p == 1) {
			    die "Unexpected resfile entry" if ($p[0] ne "NOTAA");
			    my $excl=$p[1];
			    $Nmol->setResDesType($c->{id},$mut->{resnum},
			      "PIKAA  ".&Sequence::notRes($excl));
			}
		    }
		}
	    }
	}
    }

    $Nmol->changeDesType("NOTAA","NATAA");

    # If a chain ID is specified, remove all but the target chain
    if ((defined $cid) && ($cid ne "")) {
	foreach my $c ( @{$Nmol->activeChains()} ) {
	    $Nmol->removeChain($c->{id}) if ($c->{id} ne $cid);
	}
    }

    my $outdir="tmp_computeEne_dir";
    $Nmol->setParameter(outdir=>$outdir);
    $Nmol->setParameter(ex_rot=>2);
    $Nmol->setParameter(use_input_sc=>0);

    my $opt_out=$Nmol->design();
    my $Rmol=RosettaMolecule->new($opt_out->[0]->{filename});

    if ($debug) {
	my $newdirname="tmp_computeEne_dir";
	$newdirname.="_".$mode;
	$newdirname.="_".$cid if (defined $cid);
	foreach my $mut ( @{$mutlist} ) {
	    $newdirname.="_".$mut->{resnum};
	}
	&GenUtil::rename($outdir,$newdirname);
    } else {
	# Clean up
	&GenUtil::remove($outdir);
    }

    return $Rmol->{score}->{score};

}


1;

