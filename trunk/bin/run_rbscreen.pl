#!/usr/bin/perl
use File::stat;
use lib "$ENV{'RBT_ROOT'}/lib";
use strict;
require "run_rbfuncs.pl";

my $rbtroot = $ENV{'RBT_ROOT'};
my $rbthome = $ENV{'RBT_HOME'};
my $rbtligdb = $ENV{'RBT_LIGDB'};
my $recepDir = "$rbtroot/data/receptors";
my $scriptDir = "$rbtroot/data/scripts";
my $filterDir = "$rbtroot/data/filters";
my $libraryDir = "$rbtligdb/dock";
my $binDir = "$rbtroot/bin";

#override default directories
#Parse command line arguments
my $nArgs = scalar(@ARGV);

while (scalar(@ARGV)) {
  my $arg = shift @ARGV;
  if (index($arg,'-r')==0) {
      $recepDir = substr($arg,2);
  }
  elsif (index($arg,'-s')==0) {
      $scriptDir = substr($arg,2);
  }
  elsif (index($arg,'-l')==0) {
      $libraryDir = substr($arg,2);
  }
}


# Get the receptor, library and script tables
#my %receptor_table = get_prm_table($recepDir);
my %receptor_table = get_prm_table($rbthome);
my %script_table = get_prm_table($scriptDir);
my %library_table = get_dock_table($libraryDir);

#GET ALL THE INPUTS
#my $library = get_selection(\%library_table,"library");
my $isOK;
my @libraries;
my $receptor;
my $script;

do {
  @libraries = get_multiple_selection(\%library_table,"library");
  $isOK = get_input("\nIs this OK (Y/N)","Y");
} until ($isOK eq "Y");

do {
  $receptor = get_selection(\%receptor_table,"receptor");
  $isOK = get_input("\nIs this OK (Y/N)","Y");
} until ($isOK eq "Y");

# (SJ) check whether all the required files are in the dir
# DM 7 Feb 2002 - get PSF/CRD filenames from PRM file
print "\nLooking for all required target files in $rbthome...\n\n";
my $prm_file = "$rbthome/".$receptor;
my $prm_root = substr($prm_file, 0, rindex($prm_file,"."));
my $as_file = "${prm_root}.as";
my $vdw1_file = "${prm_root}_vdw1.grd";
my $vdw5_file = "${prm_root}_vdw5.grd";

#Set the default file names equal to their parameter names so that
#error messages make sense if the parameters are missing completely from
#the .prm file.
my $psf_file = "RECEPTOR_TOPOL_FILE";
my $crd_file = "RECEPTOR_COORD_FILE";
my $mol2_file;

# DM 12 May 2003
# Optional target files for experimental restraints
#
# These default restraint suffixes and filenames should be synchronised
# with the default parameter values in the relevant C++ class source files ($RBT_ROOT/fw/src/.cxx)
# They are used if the user does not explicitly define the parameter in the
# receptor .prm file
#
# Note: this code searches for specific SECTION names in the receptor .prm file
# thus forcing a naming convention for the restraint terms that is not mandatory in the C++ code itself
#
# Pharmacophore restraint file
# SECTION PHARMA
# CONSTRAINTS_FILE = explicit filename
# OPTIONAL_FILE = explicit filename
my $ph4_filename = ".const";#synchronise with RbtPharmaSF::_CONSTRAINTS_FILE
my $ph4_optfilename = "";
my $has_ph4 = 0;
# Tether restraint file
# SECTION TETHER
# REFERENCE_FILE = suffix to the root receptor.prm filename
my $tether_suffix = "_reference.sd";#synchronise with RbtTetherSF::_REFERENCE_FILE
my $has_tether = 0;
# Nmr restraint file
# SECTION NMR
# FILENAME = explicit filename (not suffix)
my $nmr_filename = "default.noe";#synchronise with RbtTetherSF::_FILENAME
my $has_nmr = 0;

open PRMHANDLE,$prm_file || die "Cannot read $prm_file";
while (<PRMHANDLE>) {
  my @fields = split(' ',$_);
  my $nfields = scalar(@fields);
  $psf_file = "$rbthome/".$fields[1] if (($fields[0] eq 'RECEPTOR_TOPOL_FILE') && ($nfields > 1));
  $crd_file = "$rbthome/".$fields[1] if (($fields[0] eq 'RECEPTOR_COORD_FILE') && ($nfields > 1));
  $mol2_file = "$rbthome/".$fields[1] if (($fields[0] eq 'RECEPTOR_FILE') && ($nfields > 1));
  chomp $psf_file;
  chomp $crd_file;
  $has_ph4 = 1 if (($nfields > 1) && ($fields[0] eq 'SECTION') && ($fields[1] eq 'PHARMA') );
  $ph4_filename = "$rbthome/".$fields[1] if (($nfields > 1) && ($fields[0] eq 'CONSTRAINTS_FILE') && $has_ph4);
  $ph4_optfilename = "$rbthome/".$fields[1] if (($nfields > 1) && ($fields[0] eq 'OPTIONAL_FILE') && $has_ph4);
  $has_tether = 1 if (($nfields > 1) && ($fields[0] eq 'SECTION') && ($fields[1] eq 'TETHER') );
  $tether_suffix = $fields[1] if (($nfields > 1) && ($fields[0] eq 'REFERENCE_FILE') && $has_tether);
  $has_nmr = 1 if (($nfields > 1) && ($fields[0] eq 'SECTION') && ($fields[1] eq 'NMR') );
  $nmr_filename = "$rbthome/".$fields[1] if (($nfields > 1) && ($fields[0] eq 'FILENAME') && $has_nmr);
}
close PRMHANDLE;

my @src_files = ($prm_file,$as_file, $vdw1_file, $vdw5_file);
# check for either MOL2 or PSF/CRD
if ($mol2_file ne "") {
  push (@src_files,"$mol2_file");
}
else {
  push (@src_files,"$psf_file");
  push (@src_files,"$crd_file");
}
push (@src_files,"${ph4_filename}") if ($has_ph4);
push (@src_files,"${ph4_optfilename}") if (($has_ph4) && ($ph4_optfilename ne ""));
push (@src_files,"${prm_root}${tether_suffix}") if ($has_tether);
push (@src_files,"${nmr_filename}") if ($has_nmr);

my $is_file_missing = 0;
foreach my $srcfile ( @src_files) {            # search through all the target files
    if(-e $srcfile) {
        print "Found: $srcfile";
        my $sb = stat($srcfile);
        #if ($sb->size < 1024) {print " WARNING size is suspiciously small ",$sb->size};
        print "\n";
    }
    else {
        print "ERROR: $srcfile is missing\n";
        $is_file_missing    = 1;
    }
}

if ($is_file_missing) {
  print "\nCannot continue - please generate all required target files and place in $rbthome\n";
  exit;
}

do {
  $script = get_selection(\%script_table,"docking script");
  $isOK = get_input("\nIs this OK (Y/N)","Y");
} until ($isOK eq "Y");

my $flags = get_input("Enter optional rbdock flags","-ap -an");
my $mode;
do
{
  $mode = get_input("High Throughput (HT) or Exhaustive Docking (ED)", "HT");
}until (($mode eq "HT") || ($mode eq "ED"));

my $HT;
my $filter;
my $FILTERFILE;
my $nRedock;
my $nStages;
my @targets;
my @runs;
my $redockTarget;
my $redockTargetDef = -20.0;#Default score threshold for redocking phase

if ($mode eq "ED")
{
  $HT = 0;
  my $run = get_input("Enter number of runs/ligand",50);
  $flags = $flags . " -n" . $run;
  $run--;
  $filter = "1\nif - SCORE.NRUNS $run 0.0 -1.0,\n0\n";
}
else
{
  $HT = 1;
  $nStages = get_input("Number of stages (excluding redocking, enter 0 for redocking only)",2);

  #GET AND STORE THE TARGET SCORES AND NUMBER OF RUNS FOR EACH STAGE

  if ($nStages > 0) {
    for (my $n = 0; $n < $nStages; $n++) {
      print "\nDOCKING STAGE $n\n";
      my $target = get_input("Enter target score",-20.0);
      my $run = get_input("Enter number of runs/ligand",20);
      if ($run > 0)
      {
        push @targets,$target;
        push @runs,$run;
      }
    }
    $nStages = @targets;
    print "\n";
    #Default score threshold for redocking is equal to target score for last docking stage
    $redockTargetDef = $targets[$nStages-1];
  }

  $nRedock = get_input("Enter number of redocking runs/ligand",0);
  if ($nRedock > 0) {
    $redockTarget = get_input("Enter target score for redocking",$redockTargetDef);
  }
  else
  {
    $redockTarget = $redockTargetDef;
  }
  my $nPhases = $nStages;
  if ($nRedock > 0)
  {
    $nPhases++;
  }
  $filter = "$nPhases\n";
  my $nr;
  for (my $n = 0; $n < $nStages ; $n++)
  {
    $nr = $runs[$n] - 1;
    $filter .= "if - $targets[$n] SCORE.INTER 1.0 if - SCORE.NRUNS $nr 0.0 -1.0,\n";
  }
  if ($nRedock > 0)
  {
    $nr = $nRedock - 1;
    $filter .= "if - SCORE.NRUNS $nr 0.0 -1.0,\n";
  }
  my %filter_table = get_filter_table($filterDir, $redockTarget);
  my @wfilters;
  do {
    print "\nSelect filters. Only the conformations that pass all the\n",
          "filters chosen AND have SCORE.INTER lower than ", $redockTarget,
          "\nAND SCORE.RESTR.CAVITY < 1.0 will be written down\n";
    @wfilters = get_multiple_selection(\%filter_table,"filters");
    $isOK = get_input("\nIs this OK (Y/N)","Y");
  } until ($isOK eq "Y");

  for (my $i = 0 ; $i < @wfilters ; $i++) # remove the no filters
  {
    if ($wfilters[$i] eq "1no_other_filters")
    {
      splice(@wfilters, $i, 1);
    }
  }
  
  my $nwfilters = @wfilters + 2;
  $filter .= "$nwfilters\n";
  $filter .= "- SCORE.INTER $redockTarget,\n";
  $filter .= "- SCORE.RESTR.CAVITY 1.0,\n";
  foreach my $wfilter (@wfilters) {
    open(FF, "$filterDir/$wfilter") || open(FF, "./$wfilter") ||
           die "cannot open $wfilter for reading: $!";
    my $line;
    while ($line = <FF>){ $filter .= $line;}
    $filter .= ",\n";
  }
}

      
    
  my $runRoot = get_input("Enter run name","rbdock");
    
  my $pwd = `pwd`;
  chomp $pwd;

  #Make a subdirectory for each docking stage
  #plus a tmp directory
  mkdir($runRoot,0755);
  mkdir("$runRoot/tmp",0755);
  mkdir("$runRoot/tmp/condor_logs",0755);
  mkdir("$runRoot/log",0755);
  mkdir("$runRoot/sd",0755);
  mkdir("$runRoot/errors",0755);
  chdir "$runRoot/tmp";

  my $filterFile = "${runRoot}_filter.txt";
  print "Writing $filterFile...\n";
open FILTERHANDLE,">$filterFile";
print FILTERHANDLE $filter;
#WRITE THE CUSTOM DOCKING SCRIPT FOR USE BY CONDOR
#This will be the Condor "Executable"
#gzipped SD file to dock is passed as $1 argument
#Run name is passed as $2 argument
my $nqsFile = "run_${runRoot}.csh";
print "Writing $nqsFile...\n";
open NQSHANDLE,">$nqsFile";
print NQSHANDLE "#!/bin/tcsh -f\n";
print NQSHANDLE "#CONDOR EXECUTABLE SCRIPT WRITTEN BY run_rbscreen.pl\n";

# DM 14 Mar 2003
#Create a tmp runtime directory on the local disk for rDock SD and log files
#These will be compressed before copying back to the submission directory
#on the file server
#There are now two tmp directories:
# /tmp/ribodock/<version> - stores the grid files etc that are common between runs
# /tmp/${PROG}.XXXXXX - stores the I/O SD and log files unique to each run

print NQSHANDLE "set PROG = \${0:t}\n";
print NQSHANDLE "set TMPDIR = `mktemp -d /tmp/\${PROG}.XXXXXX`\n";

# synchronizing target files
my $versionidx  = index $rbtroot,"ribodock";
my $targetdir   = "/tmp/".substr($rbtroot,$versionidx);
print NQSHANDLE "mkdir -pv -m0775 /tmp/ribodock\n";
print NQSHANDLE "mkdir -pv -m0775 $targetdir\n";
print NQSHANDLE "setenv RBT_HOME $targetdir\n";
print NQSHANDLE "rsync -avz ".join(" ",@src_files)." $targetdir\n";

# Copy the compressed sd file to the tmp directory and uncompress
# DM 28 Feb 2002 - use the /tmp directory on the local disk instead
# DM 14 Mar 2003 - use TMPDIR = /tmp/${PROG}.XXXXXX instead
print NQSHANDLE "set gzfile = \$1\n";
print NQSHANDLE "set runName = \$2\n";
print NQSHANDLE "set copied_gzfile = \${TMPDIR}/\${gzfile:t}\n";
print NQSHANDLE "set sdfile = \${copied_gzfile:r}\n";
print NQSHANDLE "echo Ligand_DB : \$gzfile\n";
print NQSHANDLE "echo Local copy: \$copied_gzfile\n";
print NQSHANDLE "echo Unzipped  : \$sdfile\n";
print NQSHANDLE "cp -f \$gzfile \$copied_gzfile\n";
print NQSHANDLE "gunzip -f \$copied_gzfile\n";

my $target = $redockTarget;
print NQSHANDLE "#\n";
print NQSHANDLE "set outRoot = \${TMPDIR}/\${runName}\n";
print NQSHANDLE "set logFile = \${outRoot}.log\n";
print NQSHANDLE "if (-e \$logFile) rm -f \$logFile\n";
print NQSHANDLE "\$RBT_ROOT/bin/rbdock -i\$sdfile -o\${outRoot} -r$receptor -p$script -t$filterFile $flags > \$logFile\n";

# Compress the SD and log files and copy back to submission directory
# The ..._errors.sd files are generated by the pharmacophore restraints term for ligands that
# do not fulfill the feature requirements
print NQSHANDLE "if (-e \${outRoot}.sd) then\n";
print NQSHANDLE "  gzip -9vf \${outRoot}.sd\n";
print NQSHANDLE "  cp \${outRoot}.sd.gz ../sd\n";
print NQSHANDLE "endif\n";
print NQSHANDLE "if (-e \${logFile}) then\n";
print NQSHANDLE "  gzip -9vf \${logFile}\n";
print NQSHANDLE "  cp \${logFile}.gz ../log\n";
print NQSHANDLE "endif\n";
print NQSHANDLE "if (-e \${outRoot}_errors.sd) then\n";
print NQSHANDLE "  gzip -9vf \${outRoot}_errors.sd\n";
print NQSHANDLE "  cp \${outRoot}_errors.sd.gz ../errors\n";
print NQSHANDLE "endif\n";

print NQSHANDLE "rm -rf \${TMPDIR}\n";
close NQSHANDLE;
chmod 0755,$nqsFile;

#PREPARE THE CONDOR JOB FILE
my $condorFile = "run_${runRoot}.cmd";
print "Writing $condorFile...\n";
open CMDHANDLE,">$condorFile";
print CMDHANDLE "#CONDOR COMMAND FILE WRITTEN BY run_rbscreen.pl\n";
print CMDHANDLE "universe        = vanilla\n";
print CMDHANDLE "initialdir      = $pwd/$runRoot/tmp\n";
print CMDHANDLE "environment     = RBT_ROOT=$rbtroot;RBT_LIGDB=$rbtligdb;RBT_HOME=$rbthome;LD_LIBRARY_PATH=/usr/local/lib:$rbtroot/lib;\n";
print CMDHANDLE "executable      = ${nqsFile}\n";
#Stop Condor filling your mailbox
print CMDHANDLE "Notification    = Error\n";
print CMDHANDLE "log             = run_${runRoot}.log\n";
# Added by DM, Enspiral Discovery, 28 Feb 2005
# Prevent rDock jobs running on dedicated Catalyst nodes
print CMDHANDLE "Requirements    = CATALYST_NODE =!= True\n"; 

#GET THE LIST OF COMPRESSED SD FILES TO DOCK
my $iFile=0;
foreach my $library (@libraries) {
  my @libFiles = glob "$library_table{$library}";
  foreach my $libFile (@libFiles) {
    $iFile++;
    my $runName = "${runRoot}_$iFile";
    print CMDHANDLE "output          = condor_logs/run_${runName}.out\n";
    print CMDHANDLE "error           = condor_logs/run_${runName}.err\n";
    print CMDHANDLE "arguments       = $libFile $runName\n";
    print CMDHANDLE "queue\n";
  }
}
close CMDHANDLE;

#WRITE THE CUSTOM CLEANUP SCRIPT
#This will be run by the user from the command line
#after all the jobs have finished
my $cleanFile = "clean_${runRoot}.csh";
print "Writing $cleanFile...\n";
open CLEANHANDLE,">$cleanFile";
print CLEANHANDLE "#!/bin/tcsh -f\n";
print CLEANHANDLE "#CLEANUP SCRIPT WRITTEN BY run_rbscreen.pl\n";
print CLEANHANDLE "set sdout = ../${runRoot}_raw_hits.sd\n";
print CLEANHANDLE "if (-e \${sdout}) then\n";
print CLEANHANDLE "  echo ARE YOU SURE YOU WANT TO DO THIS\n";
print CLEANHANDLE "  echo \${sdout} already exists\n";
print CLEANHANDLE "  echo To rerun the script remove \${sdout} manually first\n";
print CLEANHANDLE "  exit\n";
print CLEANHANDLE "endif\n";
print CLEANHANDLE "echo All compressed SD files in ../sd/ will be concatenated into a single compressed file \${sdout}.gz\n";
print CLEANHANDLE "touch \${sdout}\n";
print CLEANHANDLE "foreach file (../sd/*.sd.gz)\n";
print CLEANHANDLE "  echo Sorting and concatenating \${file}...\n";
print CLEANHANDLE "  zcat \${file} | \${RBT_ROOT}/bin/sdsort -n -fSCORE -s >> \${sdout}\n";
print CLEANHANDLE "end\n";
print CLEANHANDLE "echo Each block of consecutive records per ligand have been presorted by SCORE ";
print CLEANHANDLE "ready for filtering\n";
print CLEANHANDLE "echo Compressing \${sdout}...\n";
print CLEANHANDLE "gzip -9vf \${sdout}\n";
print CLEANHANDLE "echo DONE\n";
print CLEANHANDLE "echo You may now safely remove the ../sd/ directory if you wish to save disk space\n";
close CLEANHANDLE;
chmod 0755,$cleanFile;

print "To submit, cd $runRoot/tmp; condor_submit $condorFile\n";
print "When all jobs have completed, run $runRoot/tmp/$cleanFile to concatenate the output files\n";
