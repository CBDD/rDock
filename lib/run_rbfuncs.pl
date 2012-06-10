#!/usr/bin/perl
# Perl functions for the run_rb* collection of automated docking scripts

use strict;

################################################################################
# sub get_prm_table
#
# Purpose: return a hash table of
# key=parameter file name, value = file title
# for all .prm files in a directory
#
# Usage: %prm_table = get_prm_table($dir)
#
# Arguments:
# $prmDir - directory to search for .prm files
#
# Return parameters:
# %prm_table - hash table
#
sub get_prm_table {
#Arguments
  my $prm_dir = shift;

#Local variables
  my @prm_list;
  my %prm_table;
  my $TITLEREC = "TITLE ";
  my $file;
  my $title;
  my @titleLines;

#Get the list of files ending in .prm
  @prm_list = glob "$prm_dir/*.prm";

#Read each file and extract the title record
  foreach $file (@prm_list) {
    open PRMHANDLE,$file;
    @titleLines = grep /^$TITLEREC/,<PRMHANDLE>;
    close PRMHANDLE;
#Extract just the file name from the full path
    my $prmName = substr($file, rindex($file,"/")+1);
#Check if a title record was found
    if (scalar @titleLines) {
      $title = substr $titleLines[0],length($TITLEREC);
      chomp $title;
      $prm_table{$prmName} = $title;
    }
    else {
      $prm_table{$prmName} = "No title";
    }
  }
  return %prm_table;
}
################################################################################
# sub get_dock_table
#
# Purpose: Modified version of get_prm_table
# specific for the RBT_LIGDB docking database
# Return a hash table of
# key=docking library, value = expression for all sd files in that library
#
# Assumes that all subdirs within $db_dir are separate vendor libraries,
# and that each subdir contains a set of compressed .sd.gz files
#
# Usage: %dock_table = get_dock_table($dir)
#
# Arguments:
# $db_dir - directory to search for docking libraries
#
# Return parameters:
# %dock_table - hash table#
#
sub get_dock_table {
#Arguments
  my $db_dir = shift;

#Local variables
  my ($lib,$sdfiles);
  my (@lib_list,@sd_list);
  my %dock_table;

  opendir DBDIR,$db_dir;
  @lib_list = readdir DBDIR;
  closedir DBDIR;

  foreach $lib (sort @lib_list) {
    if (-d "$db_dir/$lib") {
      $sdfiles = "$db_dir/$lib/*.sd*";
      @sd_list = glob $sdfiles;
      my $n = scalar(@sd_list);
      if ($n > 0) {
        $dock_table{"$lib ($n files)"} = $sdfiles;
      }
    }
  }

  return %dock_table;
}
#
################################################################################
# sub get_selection
#
# Purpose: allow a user to select from
# a list of parameter files.
#
# Usage: $receptor = get_selection(\%prm_table,"receptor")
#
# Arguments:
# $prm_table_ref - reference to hash table returned by get_prm_table
# $name          - descriptive name for items (e.g. "receptor" or "script"
#
# Return parameters:
# $item - selected items (key into %prm_table)
#
sub get_selection {
#Arguments
  my $prm_table_ref = shift;
  my $name = shift;

#Local variables
  my @items = sort keys %$prm_table_ref;
  my $nItems = scalar @items;
  my ($i,$itemNum,$item);
  my $inRange;

  print "\n\n$name selection:\n\n";
  $i=1;
  foreach $item (@items) {
    print "$i\t$item\t$$prm_table_ref{$item}\n";
    $i++;
  }
  do {
    print "\nEnter the $name number: ";
    $itemNum = <STDIN>;
    chomp $itemNum;
    print "\n";
    $inRange = (($itemNum >= 1) && ($itemNum <= $nItems));
    print "$itemNum is out of range\n" if (!$inRange);
  } until ($inRange);

  print "\nYou have selected the following $name\n";
  $item = $items[$itemNum-1];
  print "$itemNum\t$item\t$$prm_table_ref{$item}\n";
  return $item;
}

################################################################################
# sub get_multiple_selection
#
# Purpose: allow a user to multiply select from the docking library list
#
# Usage: $receptor = get_selection(\%prm_table,"receptor")
#
# Arguments:
# $prm_table_ref - reference to hash table returned by get_prm_table
# $name          - descriptive name for items (e.g. "receptor" or "script"
#
# Return parameters:
# $item - selected items (key into %prm_table)
#
sub get_multiple_selection {
#Arguments
  my $prm_table_ref = shift;
  my $name = shift;

#Local variables
  my @items = sort keys %$prm_table_ref;
  my $nItems = scalar @items;
  my ($i,$idstring,$itemNum,$item);
  my @itemNums;
  my @selectedItems;
  my $inRange;
  my $allInRange;

  print "\n\n$name selection:\n\n";
  $i=1;
  foreach $item (@items) {
    print "$i\t$item\t$$prm_table_ref{$item}\n";
    $i++;
  }
  do {
    print "\nEnter the $name number(s): ";
    my $idstring = <STDIN>;
    chomp $idstring;
    @itemNums = get_ids($idstring);
    print "\n";
    $allInRange = 1;
    foreach $itemNum (@itemNums) {
      $inRange = (($itemNum >= 1) && ($itemNum <= $nItems));
      print "$itemNum is out of range\n" if (!$inRange);
      $allInRange = $allInRange && $inRange;
    }
  } until ($allInRange);

  print "\nYou have selected the following $name(s)\n";
  foreach $itemNum (@itemNums) {
    $item = $items[$itemNum-1];
    push @selectedItems,$item;
    print "$itemNum\t$item\t$$prm_table_ref{$item}\n";
  }
  return @selectedItems;
}

################################################################################
# sub get_input
#
# Purpose: get user input, or returns default if no response given
#
# Usage: $nRuns = get_input("Enter no. of runs per ligand",10)
#
# Arguments:
# $question    - text of question to ask
# $defResponse - default answer
#
# Return parameters:
# $response - user response, or default

sub get_input {
#Arguments
  my $question = shift;
  my $defResponse = shift;
  print "$question [$defResponse]: ";
  my $response = <STDIN>;
  chomp $response;
  $response = $defResponse if ($response eq "");
  return $response;
}

# Based on //depot/intranet/1.0/lib/rbt_func.pl#1
# this function converts a list of ids in one string
# in the format 1,2,4-6,8,15-20
# and returns an @array with the unique ids
sub get_ids {
  my $idstring = shift;
  my %ids;

  foreach my $id (split (',',$idstring)){
    if ( grep (/-/,$id)){
      (my $low, my $up) = split ('-',$id);
      for ( my $co = $low; $co <= $up; $co++ ){
	$ids{$co}=$co;
      }
    }
    else {
      $ids{$id}=$id;
    }
  }
  
  return (sort {$a<=>$b} keys %ids);
}
################################################################################
# sub get_filter_table
#
# Purpose: return a hash table of
# key=filter file name, value = file title
# for all .filter files in a directory
#
# Usage: %filter_table = get_filter_table($dir)
#
# Arguments:
# $filterDir - directory to search for .filter files
# $tmp - temperature for first filter
# Note:
# The first two filters are not in a file, but are created here
#
# Return parameters:
# %filter_table - hash table
#
sub get_filter_table {
#Arguments
  my $filter_dir = shift;
  my $tmp = shift;

#Local variables
  my $pwd = $ENV{"PWD"};
  my @filter_list;
  my %filter_table;
  my $file;

#Get the list of files ending in .filter
  @filter_list = ((glob "$pwd/*.filter"), (glob "$filter_dir/*.filter"));
  $filter_table{"1no_other_filters"} = "";
  # $filter_table{"1threshold"} = "\tSCORE.INTER < $tmp";
  # $filter_table{"2cavity"} = "\t\tSCORE.RESTR.CAVITY < 1 ";

#Read each file and extract the title record
  foreach $file (@filter_list) {
#Extract just the file name from the full path
    my $filterName = substr($file, rindex($file,"/")+1);
#Check if a title record was found
    $filter_table{$filterName} = "No title";
  }
  return %filter_table;
}
1;
