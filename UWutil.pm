# General utilities useful in the Bakerlab environment
# 
# 2004, John Karanicolas, Baker group, UW

package UWutil;

require 5.004;

use strict;

use GenUtil;
use PDButil;
use lib "/Users/johnk/laptop/";

use autouse ImportSimple => qw( SimpleGet );

use vars qw ( $jobFilename $CondorCpus $frontend $emailAddress $sdscAccount );

BEGIN {

  if ((exists $ENV{'FRONTEND'}) && ($ENV{'FRONTEND'} ne "")) {
    $frontend=$ENV{'FRONTEND'};
  }

  $CondorCpus=0;
  if (exists $ENV{'CONDOR_CPUS'}) {
    $CondorCpus=$ENV{'CONDOR_CPUS'}+0;
    if ($CondorCpus>0) {
      if (! defined $frontend ) {
	die "Please set environment variable FRONTEND to use Condor";
      }
    }
  }

  my $user="";
  if (exists $ENV{'USER'}) {
      $user=$ENV{'USER'};
  }
  if ($user eq "johnk") {
      $user="johnk4";
  } elsif ($user eq "luke") {
      $user="lajoachi";
  }

  $emailAddress=$user."\@u.washington.edu";

  $sdscAccount = "WAS227";

}


## function: reportEmailAddress()
## Return the email address of the user

sub reportEmailAddress {

  return $emailAddress;
}


## function: runCmd(cmd)
## Run the system command passed in

sub runCmd {
  my $cmd=shift;
  my $stdOut=shift;
  my $stdErr=shift;

  if ((defined $frontend) && (($frontend =~ "teragrid") || ($frontend =~ "bluegene") || ($frontend =~ "datastar"))) {
    printf STDERR "ERROR: inappropriate to run a single-node job here\n";
    exit(0);
  }

  $stdOut="std.out" if (! defined $stdOut);
  $stdErr="std.err" if (! defined $stdErr);

  $cmd.=" > ".$stdOut;
  $cmd.=" >& ".$stdErr;
  system($cmd);

  return;
}


## function: runCmdCluster(cmd [,cpus])
## Run the system command passed in, through Cluster

sub runCmdCluster {
  my $cmd=shift;
  my $cpus=shift;
  my $outdir=shift;
  my $stdOut=shift;
  my $stdErr=shift;
  my $allowTracking=shift;
  my $highMemRequirements = shift;

  die "Error: no job to run!" if ( ( ! defined $cmd) || ($cmd eq "") );

  $allowTracking=1 if ( ! defined $allowTracking);

  if ((defined $frontend) && ($frontend =~ "teragrid")) {
    &runCmdTeragrid($cmd, $cpus, $stdOut, $stdErr, $allowTracking, $highMemRequirements);
  } elsif ((defined $frontend) && ($frontend =~ "bluegene")) {
    &runCmdBluegene($cmd, $cpus, $stdOut, $stdErr, $allowTracking, $highMemRequirements);
  } elsif ((defined $frontend) && ($frontend =~ "datastar")) {
    &runCmdDatastar($cmd, $cpus, $stdOut, $stdErr, $allowTracking, $highMemRequirements);
  } else {
    $cpus=$CondorCpus if (! defined $cpus);
    if ((defined $frontend) && ($cpus > 0)) {
      &runCmdCondor($cmd, $cpus, $outdir, $stdOut, $stdErr, $allowTracking, $highMemRequirements);
    } else {
      &runCmd($cmd);
    }
  }

  return;
}


## function: nodesFromBF()
## Check for backfill opportunities

sub nodesFromBF {
  my $cpus=shift;

  my $bf_file="bf.out";
  my $bfcmd="show_bf > ".$bf_file;
  system($bfcmd);

  my ( $nodes, $minutes, $hours );

  # If the number of cpus is defined (generally the case), overwrite this with
  # a backflow opportunity. If the number of cpus is NOT defined, use the whole cluster
  if ( defined $cpus ) {

    my $bfin=&GenUtil::getInputFile($bf_file);
    my $got_opportunity=0;
    $nodes=0;
    while (<$bfin>) {
      my $inline=$_;
      chomp($inline);
      # Take the opportunity with the largest number of nodes
      if ( $inline =~ /^  Opportunity # [0-9]+:[ ]+([0-9]+) nodes for[ ]+([0-9]+):([0-9]+) hours/) {
	my $new_nodes=$1;
	my $new_hours=$2;
	my $new_minutes=$3;
	if ( $new_nodes > $nodes ) {
	  $nodes=$new_nodes;
	  $hours=$new_hours;
	  $minutes=$new_minutes;
	}
	$got_opportunity=1;
      }
    }
    undef $bfin;

    if ( $got_opportunity) {
      printf STDERR "Found a backflow opportunity for %d nodes for %d hours and %d minutes\n", $nodes, $hours, $minutes;
      $minutes -= 5;
      if ( $minutes < 0) {
	$hours--;
	$minutes += 60;
      }
    } else {
      printf STDERR "No backflow opportunities found\n";
      $nodes=20;
      $hours=17;
      $minutes=55;
      printf STDERR "Requesting %d nodes for %d hours and %d minutes\n", $nodes, $hours, $minutes;
    }

  } else {
    # If the number of cpus is NOT specified (ie. for huge jobs), use the whole cluster
    $nodes = &reportClusterCpus() / &reportClusterPpn();
    $hours=11;
    $minutes=59;
    printf STDERR "Requesting the whole cluster (%d nodes) for %d hours and %d minutes\n", $nodes, $hours, $minutes;
  }

  return ( $nodes, $minutes, $hours );
}



## function: runCmdTeragrid(cmd [,cpus])
## Run the system command passed in, on Teragrid

sub runCmdTeragrid {
  my $cmd=shift;
  my $cpus=shift;
  my $stdOut=shift;
  my $stdErr=shift;
  my $allowTracking=shift;

#  if ( ( defined $cpus ) && ( $cpus < 2 ) ) {
#    printf STDERR "Teragrid should not be used for jobs requiring less than 2 cpus!\n";
#    exit(0);
#  }

  $stdOut="teragrid.out" if (! defined $stdOut);
  $stdErr="teragrid.err" if (! defined $stdErr);

  my $ppn=&reportClusterPpn();
  my ( $nodes, $minutes, $hours ) = &nodesFromBF($cpus);

  my $PBSjobFname = "teragrid.job";
  my $cout=&GenUtil::getOutputFile($PBSjobFname);

  my $pwd=$ENV{PWD};
  my @n=split(/\//,$pwd);
  my $jobname=pop(@n);
  ( $jobname = $jobname ) =~ s/^[0-9]+//gi;

  printf $cout "\#SDSC Teragrid PBS Script\n";
  printf $cout "\#PBS -o %s\n", $stdOut;
  printf $cout "\#PBS -e %s\n", $stdErr;
  printf $cout "\#PBS -l nodes=%d\:ppn=%d\n", $nodes, $ppn;
  printf $cout "\#PBS -l walltime=%s\:%s\:00\n", &GenUtil::zPad($hours,2), &GenUtil::zPad($minutes,2);
  printf $cout "\#PBS -q dque\n";
  printf $cout "\#PBS -V\n";
  printf $cout "\#PBS -M %s\n", $emailAddress;
  printf $cout "\#PBS -A %s\n", $sdscAccount;
  printf $cout "\#PBS -N %s\n", $jobname;
  printf $cout "\n";

  my $sfile = &reportJobFilename();

  printf $cout "cd %s\n",$pwd;
# Note: SDSC PBS doesn't give the "$NP" variable, so we'll set it ourselves...
  printf $cout "set NP = \`wc -l \< \$PBS_NODEFILE\`\n\n",$sfile;
  printf $cout "chmod 744 ./%s\n",$sfile;
  printf $cout "./%s \$PBS_NODEFILE \$NP \$PWD\n\n",$sfile;
  undef $cout;

  my $sout=&GenUtil::getOutputFile($sfile);
  printf $sout "\#!/bin/tcsh\n\n";
  printf $sout "set nodefile = \$argv[1]\n";
  printf $sout "set np = \$argv[2]\n";
  printf $sout "set dir = \$argv[3]\n";
  printf $sout "\n";
  printf $sout "cd \$dir\n";
  printf $sout "mpirun -v -machinefile \$nodefile -np \$np %s\n\n",$cmd;
  undef $sout;

  my $submit_cmd="qsub ".$PBSjobFname;
  $submit_cmd.=" > qsub.stdout" if ( $allowTracking );

  my $syscmd="ssh ".$frontend." \"";
  $pwd=$ENV{PWD};
  chomp $pwd;
  $syscmd.="cd ".$pwd."; ";
  $syscmd.=$submit_cmd;
  $syscmd.="\"";
  system($syscmd);

  printf STDERR "Job submitted - perl exiting\n";
  exit(0);

  return;
}


## function: runCmdBluegene(cmd [,cpus])
## Run the system command passed in, on Bluegene

sub runCmdBluegene {
  my $cmd=shift;
  my $cpus=shift;
  my $stdOut=shift;
  my $stdErr=shift;
  my $allowTracking=shift;
  my $highMemRequirements = shift;

  $highMemRequirements=0 if (! defined $highMemRequirements);

  # Choose random node set (top64-1 to 8, bot64-1 to 8)
  my $rackname="bot";
  my $randID = int( rand(16) ) + 1;    # get an int between 1 and 16 (inclusive)
  if ( $randID > 8 ) {
    $rackname="top";
    $randID -= 8;
  }

  my $node_set = $rackname."64-".$randID;
  my $mode="CO";
  $cpus = 64;

#  if ( ! $highMemRequirements ) {
#    $mode="VN";
#    $cpus *= 2;
#  }

  my $hours = 11;
  my $minutes = 59;

  $stdOut="bluegene.out" if (! defined $stdOut);
  $stdErr="bluegene.err" if (! defined $stdErr);

  my $jobFname = &reportJobFilename();
  my $cout=&GenUtil::getOutputFile($jobFname);

  my @c=split(' ',$cmd);
  my $exec=shift(@c);
  my $sarg=join(" ",@c);

  my $pwd=$ENV{PWD};
  my @n=split(/\//,$pwd);
  my $jobname=pop(@n);
#  ( $jobname = $jobname ) =~ s/^[0-9]+//gi;

  printf $cout "\#\!/usr/bin/ksh\n";
  printf $cout "\#\@bg_partition = %s\n", $node_set;
  printf $cout "\#\@environment = COPY_ALL\;\n";
  printf $cout "\#\@job_type = BlueGene\n";
  printf $cout "\#\@class = parallel\n";
  printf $cout "\#\@output = %s\n", $stdOut;
  printf $cout "\#\@error = %s\n", $stdErr;
  printf $cout "\#\@account_no = %s\n", $sdscAccount;
  printf $cout "\#\@notification = complete\n";
  printf $cout "\#\@notify_user = %s\n", $emailAddress;
  printf $cout "\#\@wall_clock_limit = %s\:%s\:00\n", &GenUtil::zPad($hours,2), &GenUtil::zPad($minutes,2);
  printf $cout "\#\@job_name = %s\n", $jobname;
  printf $cout "\#\@queue\n";
  printf $cout "\n\n";

  printf $cout "mpirun -mode %s -np %s -exe %s -cwd %s\/ -args \" %s \"\n", $mode, $cpus, $exec, $pwd, $sarg;

  printf $cout "\n\n";
  undef $cout;

#  my $submit_cmd="llsubmit ".$jobFname;
#  $submit_cmd.=" > llsub.stdout" if ( $allowTracking );
#  system($submit_cmd);

  printf STDERR "Job NOT submitted - choose nodes to use, and submit manually!\n";
  exit(0);

  return;
}


## function: runCmdDatastar(cmd [,cpus])
## Run the system command passed in, on Datastar

sub runCmdDatastar {
  my $cmd=shift;
  my $cpus=shift;
  my $stdOut=shift;
  my $stdErr=shift;
  my $allowTracking=shift;

#  if ( ( defined $cpus ) && ( $cpus < 2 ) ) {
#    printf STDERR "Datastar should not be used for jobs requiring less than 2 cpus!\n";
#    exit(0);
#  }

  $stdOut="datastar.out" if (! defined $stdOut);
  $stdErr="datastar.err" if (! defined $stdErr);

  my $ppn=&reportClusterPpn();
  my ( $nodes, $minutes, $hours ) = &nodesFromBF($cpus);

  my $jobFname = "datastar.job";
  my $cout=&GenUtil::getOutputFile($jobFname);

  my $pwd=$ENV{PWD};
  my @n=split(/\//,$pwd);
  my $jobname=pop(@n);
#  ( $jobname = $jobname ) =~ s/^[0-9]+//gi;

  printf $cout "\#\!/usr/bin/ksh\n";
  printf $cout "\#\#\#\n";
  printf $cout "\#\@wall_clock_limit = %s\:%s\:00\n", &GenUtil::zPad($hours,2), &GenUtil::zPad($minutes,2);
  printf $cout "\#\@node = %d\n", $nodes;
  printf $cout "\#\@tasks_per_node = %d\n", $ppn;
  printf $cout "\#\@job_name = %s\n", $jobname;
  printf $cout "\#\@initialdir = %s\n", $pwd;
  printf $cout "\#\#\n";

  printf $cout "\#\@output = %s\n", $stdOut;
  printf $cout "\#\@error = %s\n", $stdErr;
  printf $cout "\#\@notification = never\n";
  printf $cout "\#\@notify_user = %s\n", $emailAddress;
  printf $cout "\#\@account_no = %s\n", $sdscAccount;
  printf $cout "\#\#\n";

  printf $cout "\#\@environment = COPY_ALL\; \\\n";
  printf $cout "\#AIXTHREAD_SCOPE=S\; \\\n";
  printf $cout "\#MP_ADAPTER_USE=dedicated\; \\\n";
  printf $cout "\#MP_CPU_USE=unique\; \\\n";
  printf $cout "\#MP_CSS_INTERRUPT=no\; \\\n";
  printf $cout "\#MP_EAGER_LIMIT=64K\; \\\n";
  printf $cout "\#MP_EUIDEVELOP=min\; \\\n";
  printf $cout "\#MP_LABELIO=yes\; \\\n";
  printf $cout "\#MP_POLLING_INTERVAL=100000\; \\\n";
  printf $cout "\#MP_PULSE=0\; \\\n";
  printf $cout "\#MP_SHARED_MEMORY=yes\; \\\n";
  printf $cout "\#MP_SINGLE_THREAD=yes\; \\\n";
  printf $cout "\#RT_GRQ=ON\; \\\n";
  printf $cout "\#SPINLOOPTIME=0\; \\\n";
  printf $cout "\#YIELDLOOPTIME=0\n";
  printf $cout "\#\; \\\n";
  printf $cout "\#\@network.MPI = sn_all, shared, US\n";
  printf $cout "\#\#\n";

  printf $cout "\#\@job_type = parallel\n";
  printf $cout "\#\@class = normal\n";
  printf $cout "\#\@node_usage = not_shared\n";
  printf $cout "\#\@queue = normal\n";
  printf $cout "\#\#\n";

  my $sfile = &reportJobFilename();

  printf $cout "\n\n";
  printf $cout "chmod 744 ./%s\n",$sfile;
  printf $cout "./%s\n\n",$sfile;
  undef $cout;

  my $sout=&GenUtil::getOutputFile($sfile);
  printf $sout "\#!/bin/tcsh\n\n";
  printf $sout "if \( ! -d lib \) then\n";
  printf $sout "  \\ln -s ~/working_rosetta++/lib\n";
  printf $sout "endif\n\n";
  printf $sout "poe %s\n\n",$cmd;
  undef $sout;

  my $submit_cmd="llsubmit ".$jobFname;
  $submit_cmd.=" > llsub.stdout" if ( $allowTracking );
  system($submit_cmd);

  printf STDERR "Job submitted - perl exiting\n";
  exit(0);

  return;
}


## function: runCmdCondor(cmd [,cpus])
## Run the system command passed in, via Condor

sub runCmdCondor {
  my $cmd=shift;
  my $cpus=shift;
  my $outdir=shift;
  my $stdOut=shift;
  my $stdErr=shift;
  my $allowTracking=shift;
  my $highMemRequirements = shift;

  $stdOut="condor.out" if (! defined $stdOut);
  $stdErr="condor.err" if (! defined $stdErr);

  $cpus=$CondorCpus if (! defined $cpus);

  die "Error: cannot run condor without a frontend"
    if ( ! defined $frontend);

  die "Error: cannot run condor with no cpus"
    if ($cpus < 1);

  $cpus=$CondorCpus if ($CondorCpus < $cpus);

  my @c=split(' ',$cmd);
  my $exec=shift(@c);
  my $sarg=join(" ",@c);


  if ( $exec =~ /rosetta/) {
#    $sarg .= " -delay_at_start -seed_offset \$(Process)";   # deprecated
    $sarg .= " -concatenate_pdbs -concatenate_pdbs_filename \$\(Process\)\/outpdbs" .
      " -concatenate_pdbs_listname ".$outdir."\$\(Process\)\/outlist.txt";

    # Append the process ID to the scorefile and to "infiles.list"
    my $old_scorefile_string="scorefile scorefile";
    my $new_scorefile_string="scorefile \$\(Process\)\/scorefile";
    ( $sarg = $sarg ) =~ s/$old_scorefile_string/$new_scorefile_string/;

    my $old_infiles_string="infiles.list";
    my $new_infiles_string="infiles.list\$\(Process\)";
    ( $sarg = $sarg ) =~ s/$old_infiles_string/$new_infiles_string/;

  }

  my $jobFname = &reportJobFilename();
  my $cout=&GenUtil::getOutputFile($jobFname);
  printf $cout "\#Simple condor script\n";
  printf $cout "Executable   =   %s\n",$exec;
  printf $cout "Universe     =   vanilla\n";

  if ( $cpus > 1 ) {
#    printf $cout "Output       =   %s\$\(Process\)\/%s\n",$outdir,$stdOut;
#    printf $cout "Error        =   %s\$\(Process\)\/%s\n",$outdir,$stdErr;
#    printf $cout "Log          =   %s\$\(Process\)\/condor.log\n",$outdir;

    # send output to /dev/null, for syd compliance
    printf $cout "Output       =   /dev/null\n";
    printf $cout "Error        =   /dev/null\n";
    printf $cout "Log          =   /dev/null\n";
  } else {
#    printf $cout "Output       =   %s%s\n",$outdir,$stdOut;
#    printf $cout "Error        =   %s%s\n",$outdir,$stdErr;
#    printf $cout "Log          =   %scondor.log\n",$outdir;

    # send output to /dev/null, for syd compliance
    printf $cout "Output       =   /dev/null\n";
    printf $cout "Error        =   /dev/null\n";
    printf $cout "Log          =   /dev/null\n";

  }

  printf $cout "Getenv       =   True\n";
#  printf $cout "Requirements =   Memory >= 1999\n" if ( $highMemRequirements );
  printf $cout "Notification =   Never\n";
  printf $cout "MaxHosts     =   100\n";
  printf $cout "Arguments    =   %s\n",$sarg;
  printf $cout "Initialdir   =   ./\n";
  printf $cout "Queue %d\n",$cpus;
  undef $cout;

  my $submit_cmd="condor_submit ".$jobFname;
  $submit_cmd.=" > condor.stdout" if ( $allowTracking );

  my $syscmd="";
  $syscmd.="ssh ".$frontend." \"";
  my $pwd=$ENV{PWD};
  chomp $pwd;
  $syscmd.="cd ".$pwd."; ";
  $syscmd.=$submit_cmd;
  $syscmd.="\"";
  system($syscmd);

  printf STDERR "Submitted condor job, now exiting\n";
  exit(0);

  return;
}


## function: setJobFilename()
## Set the filename jobs will have

sub setJobFilename {
  my $newJobFilename=shift;

  if ( defined $newJobFilename ) {
    $jobFilename = $newJobFilename;
  } else {
    if ((defined $frontend) && ($frontend =~ "teragrid")) {
      $jobFilename = "teragrid.csh";
    } elsif ((defined $frontend) && ($frontend =~ "datastar")) {
      $jobFilename = "datastar.csh";
    } elsif ((defined $frontend) && ($frontend =~ "bluegene")) {
      $jobFilename = "bluegene.job";
    } else {
      $jobFilename = "condor.job";
    }
  }

  return;
}


## function: reportJobFilename()
## Return the filename jobs will have

sub reportJobFilename {

  &setJobFilename() if ( ! defined $jobFilename );
  return $jobFilename;
}


## function: reportClusterCpus()
## Return the number of available cpus

sub reportClusterCpus {

  return 1 if (! defined $frontend);

  if (($frontend =~ "teragrid") || ($frontend =~ "datastar")) {
    return 504;
  } else {
    return &reportCondorCpus();
  }

  return $CondorCpus;
}


## function: reportClusterPpn()
## Return the number of cpus per node

sub reportClusterPpn {

  if ((defined $frontend) && ($frontend =~ "teragrid")) {
    return 2;
  } elsif ((defined $frontend) && ($frontend =~ "datastar")) {
    return 8;
  }

  return 1;
}


## function: reportCondorCpus()
## Return the number of CondorCpus

sub reportCondorCpus {

  return $CondorCpus;
}


## function: setCondorCpus(num)
## Set the number of CondorCpus (overrides the environment variable)

sub setCondorCpus {
  my $num=shift;
  $CondorCpus=$num if (defined $num);

  return;
}


## function: waitForClusterCompletion(num)

sub waitForClusterCompletion {
  my $Rmol=shift;

  if ((defined $frontend) && ($frontend =~ "teragrid")) {
    printf STDERR "BEWARE: teragrid version of waitForClusterCompletion not yet implemented\n";
    printf STDERR "Continuing without waiting....\n";
  } elsif ((defined $frontend) && ($frontend =~ "datastar")) {
    printf STDERR "BEWARE: datastar version of waitForClusterCompletion not yet implemented\n";
    printf STDERR "Continuing without waiting....\n";
  } elsif ((defined $frontend) && ($frontend =~ "bluegene")) {
    printf STDERR "BEWARE: BlueGene version of waitForClusterCompletion not yet implemented\n";
    printf STDERR "Continuing without waiting....\n";
  } else {
    &waitForCondorCompletion($Rmol);
  }

  return;
}

## function: checkForClusterJob(jobid)

sub checkForClusterJob {
  my $jobid=shift;

  return 0 if ( ! defined $jobid);

  if ((defined $frontend) && ($frontend =~ "teragrid")) {
    return &checkForTeragridJob($jobid);
  } elsif ((defined $frontend) && ($frontend =~ "datastar")) {
    printf STDERR "Checking for running datastar jobs is not yet implemented\n";
    printf STDOUT "Checking for running datastar jobs is not yet implemented\n";
    exit(1);
  } elsif ((defined $frontend) && ($frontend =~ "bluegene")) {
    printf STDERR "Checking for running BlueGene jobs is not yet implemented\n";
    printf STDOUT "Checking for running BlueGene jobs is not yet implemented\n";
    exit(1);
  } else {
    return &checkForCondorJob($jobid);
  }

  return;
}

## function: getClusterJobID()

sub getClusterJobID {

  if ((defined $frontend) && ($frontend =~ "teragrid")) {
    return &getTeragridJobID();
  } elsif ((defined $frontend) && ($frontend =~ "datastar")) {
    printf STDERR "Checking for running datastar jobs is not yet implemented\n";
    printf STDOUT "Checking for running datastar jobs is not yet implemented\n";
    exit(1);
  } elsif ((defined $frontend) && ($frontend =~ "bluegene")) {
    printf STDERR "Checking for running BlueGene jobs is not yet implemented\n";
    printf STDOUT "Checking for running BlueGene jobs is not yet implemented\n";
    exit(1);
  } else {
    return &getCondorJobID();
  }

  return;
}


## function: getTeragridJobID()

sub getTeragridJobID {
  my $fname=shift;

  $fname="qsub.stdout" if ( ! defined $fname);
  return undef if (! &GenUtil::zexists($fname));

  my $cin=&GenUtil::getInputFile($fname);
  my $inline;
  while (<$cin>) {
    $inline=$_ if (/sdsc\.teragrid\.org/);
  }
  undef $cin;
  chomp($inline);
  my @p=split(/\./,$inline);
  my $jobid=shift(@p);

  return $jobid;
}


## function: waitForCondorCompletion([jobid])
## Wait for a condor job to finish running, taking the jobid
## from condor.stdout if not passed in

sub waitForCondorCompletion {
  my $Rmol=shift;

  return if ($CondorCpus < 1);

  my $jobid = &getCondorJobID("condor.stdout");
  return if (! defined $jobid);

  my $done=0;
  while (! $done) {
    $done=1 if ( ! &checkForCondorJob($jobid) );
    sleep(900) if (! $done);
    $Rmol->_cleanOutput() if (defined $Rmol);
  }

  &GenUtil::remove("condor.stdout");
  &GenUtil::remove("condor.status");
  &GenUtil::remove("condor.log");

  return;
}


## function: getCondorJobID(fname)
## Read the condor job ID from the file condor.stdout

sub getCondorJobID {
  my $fname=shift;

  $fname="condor.stdout" if ( ! defined $fname);
  return undef if (! &GenUtil::zexists($fname));

  my $cin=&GenUtil::getInputFile($fname);
  my $inline;
  while (<$cin>) {
    $inline=$_ if (/submitted to cluster/);
  }
  undef $cin;
  chomp($inline);
  my @p=split(' ',$inline);
  my $jobid=pop(@p);
  ( $jobid = $jobid ) =~ s/\.$//gi;

  return $jobid;
}


## function: checkForTeragridJob(jobid)

sub checkForTeragridJob {
  my $jobid=shift;

  return 0 if ( ! defined $jobid);

  my $stat_file="qsub.status";
  my $syscmd="qstat > ".$stat_file;
  system($syscmd);

  my $found=0;
  if ( &GenUtil::zexists($stat_file) ) {
    my $stat=&GenUtil::getInputFile($stat_file);
    while (<$stat>) {
      chomp;
      my @s=split(' ',$_);
      if (defined $s[0]) {
	my @p=split(/\./,$s[0]);
	if ($p[0] =~ $jobid) {
	  $found=1;
	}
      }
    }
    undef $stat;
  }

  return $found;
}


## function: checkForCondorJob([jobid])
## Check whether a specific condor job is running

sub checkForCondorJob {
  my $jobid=shift;

  return 0 if ( ! defined $jobid);

  my $user="";
  if ((exists $ENV{'USER'}) && ($ENV{'USER'} ne "")) {
      $user=$ENV{'USER'};
  }

  my $stat_file="condor.status";

  my $syscmd="ssh ".$frontend." \"";
  my $pwd=$ENV{PWD};
  chomp $pwd;
  $syscmd.="cd ".$pwd."; ";
  $syscmd.="condor_q ".$user." > ".$stat_file;
  $syscmd.="\"";
  system($syscmd);

  my $found=0;
  if ( &GenUtil::zexists($stat_file) ) {
    my $stat=&GenUtil::getInputFile($stat_file);
    while (<$stat>) {
      chomp;
      my $inline=$_;
      if ($inline =~ $user) {
	my @s=split(' ',$inline);
	shift(@s) if ((defined $s[0]) && ($s[0] !~ /\./));
	if (((" ".$s[0]) =~ / $jobid\./) && (($s[5] eq "R") || ($s[5] eq "I"))) {
	  $found=1;
	}
      }
    }
    undef $stat;
  }

  return $found;
}


## function: PDBfname(PDBid)
## Given a PDB code, return a complete path for it (on /net/wwwpdb/)
## E.g. given "1PGX", returns "/net/wwwpdb/pg/1pgx.pdb"

sub PDBfname {
  my $PDBid=shift;

  return if (! defined $PDBid);

  # Check current directory first (include path if passed in)
  return $PDBid if (-e $PDBid);
  return $PDBid."pdb" if (-e ($PDBid."pdb"));
  return $PDBid.".pdb" if (-e ($PDBid.".pdb"));
  return $PDBid."gz" if (-e ($PDBid."gz"));
  return $PDBid.".gz" if (-e ($PDBid.".gz"));
  return $PDBid."pdb.gz" if (-e ($PDBid."pdb.gz"));
  return $PDBid.".pdb.gz" if (-e ($PDBid.".pdb.gz"));

  my $prefix=&PDButil::prefix($PDBid);

  # Next check local PDB dir
  my $home=$ENV{HOME};
  my $localPDBdir=$home."/PDBs/";
  if (-d $localPDBdir) {
      my $fname=$localPDBdir.lc($prefix).".pdb.gz";
      return $fname if (-e $fname);
  }

  # Next, check check /net/pdb/
  if (-d ("/net/wwwpdb/".substr($prefix,1,2)."/")) {
      my $fname="/net/wwwpdb/".substr($prefix,1,2)."/".$prefix.".pdb";
      return $fname if (-e $fname);
  }

  # Finally, try downloading from RCSB
  # Note: this prompts autouse loading LWP::Simple
  if (-d $localPDBdir) {
      my $url = "http://www.rcsb.org/pdb/files/".$prefix.".pdb";
      my $outname=$prefix.".pdb";
      my $outpdb=&GenUtil::getOutputFile($outname);
      printf $outpdb "%s\n",&SimpleGet($url);
      undef $outpdb;
      if ( (-s $outname) > 999) {
	  &GenUtil::compress($outname);
	  $outname.=".gz";
	  &GenUtil::copy($outname,$localPDBdir.lc($outname));
	  printf STDERR "Adding %s to local PDB dir\n",$outname;
      } else {
	  &GenUtil::remove($outname);
	  printf STDERR "PDB file %s not found in RCSB\n",$outname;
	  exit(1);
      }
      return $localPDBdir.$outname;
  }

  printf STDERR "PDB file %s not found\n",$PDBid;
  exit(1);

  return;
}


## function: ScopClaName()
## return name of local copy of parsable SCOP data file

sub ScopClaName {

  # Data dir
  my $scopDir;
  if (defined $ENV{MISCDATA}) {
      $scopDir=$ENV{MISCDATA};
  } else {
      die "Environment variable MISCDATA must be set";
  }

  my @s=glob($scopDir."/dir.cla.scop.txt*");
  if ($#s > 0) {
      die "Multiple SCOP data files are present!";
  }
  die "Could not find a valid SCOP file" if (! &GenUtil::zexists($s[0]));

  printf STDERR "Using SCOP data file: %s\n",$s[0];

  return $s[0];
}


1;
