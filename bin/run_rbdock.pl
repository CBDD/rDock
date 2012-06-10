#!/usr/bin/perl
use lib "$ENV{'RBT_ROOT'}/lib";
require "run_rbfuncs.pl";

my $rbtroot = $ENV{'RBT_ROOT'};
my $rbthome = $ENV{'RBT_HOME'};
my $recepDir = "$rbtroot/data/receptors";
my $scriptDir = "$rbtroot/data/scripts";
my $libraryDir = "$rbtroot/data/libraries";
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
my %receptor_table = get_prm_table($recepDir);
my %script_table = get_prm_table($scriptDir);
my %library_table = get_prm_table($libraryDir);

#GET ALL THE INPUTS
my $defFiles = "*.{sd,sdf,mol,mdl}";
print join("\n",glob $defFiles);
my $libFiles = get_input("\n\nEnter expression for SD files to dock",$defFiles);
my $receptor = get_selection(\%receptor_table,"receptor");
my $script = get_selection(\%script_table,"docking script");
my $flags = get_input("Enter optional rbdock flags","-ap -an");
my $runs = get_input("Enter number of runs/ligand",20);

my $runRoot = get_input("Enter run name","rbdock");
my $queueSys = get_input("Enter queuing system (NQS,CODINE,CONDOR,NONE)","CONDOR");
my $queueName;
$queueName = get_input("Enter NQS queue","QRBD") if ($queueSys eq "NQS");

$pwd = `pwd`;
chomp $pwd;

#WRITE THE NQS SCRIPTS

#GET THE LIST OF SD FILES TO DOCK
my @libList = glob "${libFiles}";

my $subFile = "submit_$runRoot.csh";
open SUBHANDLE,">$subFile";
print SUBHANDLE "#!/bin/csh -ef\n";

my $iFile=0;
foreach $libFile (@libList) {
  $iFile++;
  #Write the NQS script file
  $runName = "${runRoot}_$iFile";
  $nqsFile = "run_${runName}.csh";
  my $logFile = "${runName}.log";

  print "Writing $nqsFile...\n";

  open NQSHANDLE,">$nqsFile";
  print NQSHANDLE "#!/bin/csh -ef\n";
  print NQSHANDLE "setenv RBT_ROOT $rbtroot\n";
  print NQSHANDLE "setenv RBT_HOME $rbthome\n";
  print NQSHANDLE "source \$RBT_ROOT/bin/setlibrbt.com\n";
  print NQSHANDLE "cd $pwd\n";
  print NQSHANDLE "if (-e $logFile) rm -f $logFile\n";
  print NQSHANDLE "\$RBT_ROOT/bin/rbdock -i$libFile -o${runName} -r$receptor -p$script -n$runs $flags > $logFile\n";
  close NQSHANDLE;

  chmod 0755,$nqsFile;
  print SUBHANDLE "qsub -q $queueName $nqsFile\n" if ($queueSys eq "NQS");
  print SUBHANDLE "qsub -e $pwd/$runRoot/tmp/ -o $pwd/$runRoot/tmp/ $nqsFile\n" if ($queueSys eq "CODINE");
  if ($queueSys eq "CONDOR") {
    my $condorFile = "run_${runName}.cmd";
    open CMDHANDLE,">$condorFile";
    print CMDHANDLE "#CONDOR FILE WRITTEN BY run_rbscreen.pl\n";
    print CMDHANDLE "universe        = vanilla\n";
    print CMDHANDLE "executable      = ${nqsFile}\n";
    print CMDHANDLE "Requirements    = CATALYST_NODE =!= True\n";
    print CMDHANDLE "output          = ${nqsFile}.out\n";
    print CMDHANDLE "error           = ${nqsFile}.err\n";
    print CMDHANDLE "log             = ${nqsFile}.log\n";
    print CMDHANDLE "queue\n";
    close NQSHANDLE;
    print SUBHANDLE "condor_submit $condorFile\n";
  }
}
close SUBHANDLE;
chmod 0755,$subFile;
















