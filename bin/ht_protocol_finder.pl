#!/usr/bin/perl -w

# Script that simulates the result of a high throughput protocol.

if (@ARGV < 4)
{
	print "Script that simulates the result of a high throughput protocol.\n";
	print "1st) exhaustive docking of a small representative part of the\n";
	print "     whole library.\n";
	print "2nd) Store the result of sdreport -t over that exhaustive dock.\n";
	print "     in file <sdreport_file> that will be the input of this\n";
	print "     script.\n";
	print "3rd) ht_protocol_finder.pl <sdreport_file> <output_file> <thr1max>";
	print " <thr1min> <ns1> <ns2>\n";
	print "     <ns1> and <ns2> are the number of steps in stage 1 and in\n";
	print "     stage 2. If not present, the default values are 5 and 15\n";
	print "     <thrmax> and <thrmin> setup the range of thresholds that will\n";
	print "     be simulated in stage 1. The threshold of stage 2 depends\n";
	print "     on the value of the threshold of stage 1.\n";
	print "     An input of -22 -24 will try protocols:\n";
	print "            5	-22	15	-27\n";
	print "            5	-22	15	-28\n";
	print "            5	-22	15	-29\n";
	print "            5	-23	15	-28\n";
	print "            5	-23	15	-29\n";
	print "            5	-23	15	-30\n";
	print "            5	-24	15	-29\n";
	print "            5	-24	15	-30\n";
	print "            5	-24	15	-31\n";
	print "     Output of the program is a 7 column values. First column\n";
	print "     represents the time. This is a percentage of the time it\n";
	print "     would take to do the docking in exhaustive mode, i.e. \n";
    print "     docking each ligand 100 times. Anything\n";
	print "     above 12 is too long.\n";
	print "     Second column is the first percentage. Percentage of\n";
	print "     ligands that pass the first stage.\n";
	print "     Third column is the second percentage. Percentage of\n";
	print "     ligands that pass the second stage.\n";
	print "     The four last columns represent the protocol.\n";
	print "     All the protocols tried are written at the end.\n";
	print "     The ones for which time is less than 12%, perc1 is\n";
	print "     less than 30% and perc2 is less than 5% but bigger than 1%\n";
	print "     will have a series of *** after, to indicate they are good choices\n";
	print "     WARNING! This is a simulation based in a small set.\n";
	print "     The numbers are an indication, not factual values.\n";
    exit (0);
}
open (IFILE, $ARGV[0]) || die "cannot open $ARGV[0] for reading: $!";
open (OFILE, ">$ARGV[1]") || die "cannot open $ARGV[1] for writing: $!";
<IFILE>;   # read first line
$line = <IFILE>;
$line =~ /\S+\s+(\S+)\s+\S+\s+(\S+)/;
$prevname = $1;
$inter = $2;
push @tmp, $inter; 
$i = 0;
while ($line = <IFILE>)
{
  $line =~ /\S+\s+(\S+)\s+\S+\s+(\S+)/;
  $name = $1;
  $inter = $2;
  if ($name ne $prevname)
  {
      $prevname = $name;
      push @inters, [ @tmp ];
      $runsperligand[$i++] = @tmp;
      $#tmp = -1;
  }
  push @tmp, $inter;
}
push @inters, [ @tmp ];
$runsperligand[$i++] = @tmp;
$totLigands = @inters;
$ta = $ARGV[2];
$tb = $ARGV[3];
if (@ARGV > 4)
{
	$n1 = $ARGV[4];
}
else
{
	$n1 = 5;
}
if (@ARGV > 5)
{
	$n2 = $ARGV[5];
}
else
{
	$n2 = 15;
}
printf OFILE "Command line args:\n\t";
foreach $arg (@ARGV)
{
	printf OFILE "%s\t", $arg;
}
printf OFILE "\n\n  TIME   PERC1   PERC2    N1  THR1    N2  THR2\n"; 
for ($t1 = $ta ; $t1 >= $tb ; $t1--)
{
    $tc = $t1 - 5;
    for ($t2 = $tc ; $t2 >= ($tc - 2) ; $t2--)
    {
        ($cref, $dref) = 
          getTables($t1, $t2);
        @ligBelowThr1 = @$cref;
        @ligBelowThr2 = @$dref;
        ($time,$p1,$p2) = simulation($n1,$t1,$n2,$t2);
        if (($time < 12.0) && ($p1 < 30.0) && ($p2 < 5.0) && ($p2 > 1.0))
        {
            printf OFILE "%6.3f, %6.3f, %6.3f, %4i, %4i, %4i, %4i ***\n", 
                         $time, $p1, $p2, $n1,$t1,$n2,$t2;
        }
        else
        {
            printf OFILE "%6.3f, %6.3f, %6.3f, %4i, %4i, %4i, %4i\n", 
                         $time, $p1, $p2, $n1,$t1,$n2,$t2;
        }
    }
}


sub simulation
{
  my(@params);
  my(@lt0, @lt1);
  my($thr1, $thr2,$ns1,$ns2,$tottime,$totnumHits,$total,$totnruns);
  @params = @_;
  $ns1 = $params[0];
  $thr1 = $params[1];
  $ns2 = $params[2];
  $thr2 = $params[3];
  $ntrials = 100;
  $tottime = 0;
  $totnumHitss1 = 0;
  $totnumHitss2 = 0;
  $total = $totLigands;
  $totnruns = 0;
  for ($i = 0 ; $i < $ntrials ; $i++)
  {
    $numHitss1 = 0;
    $numHitss2 = 0;
    $totnruns = 0;
    for ($j = 0 ; $j < $totLigands ; $j++)
    {
      ($passStage1, $nruns1) = stage($ns1, $ligBelowThr1[$j]);
      if ($passStage1)
      {
        $numHitss1++;
        ($passStage2, $nruns2) = stage($ns2, $ligBelowThr2[$j]);
        if ($passStage2)
        {
          $numHitss2++;
        }
        $totnruns += $nruns2;
      }
      $totnruns += $nruns1;
    }
    $time = $totnruns / $total;
    $tottime += $time;
    $totnumHitss1 += $numHitss1;
    $totnumHitss2 += $numHitss2;
  }
  $tottime /= $ntrials;
  $totnumHitss1 /= $ntrials;
  $totnumHitss2 /= $ntrials;
  $p1 = $totnumHitss1 * 100.0 / $totLigands;
  $p2 = $totnumHitss2 * 100.0 / $totLigands;
  return ($tottime, $p1, $p2);
}

sub getTables
{
  my(@params);
  my($totLigands, @lt1, @lt2);
  my($row,$inter, $belowThr1, $belowThr2);
  my($thr0, $thr1);
  @params = @_;
  $thr1 = $params[0];
  $thr2 = $params[1];
  $totLigands = 0;
  $belowThr1 = 0;
  $belowThr2 = 0;
  $totLigands = 0;
  for $i (0 .. $#inters)
  {
    $row = $inters[$i];
    for $j (0 .. $#{$row})
    {
      $inter = $row->[$j];
      if ($inter <= $thr1)
      {
        $belowThr1++;
      }
      if ($inter <= $thr2)
      {
        $belowThr2++;
      }
    }
    $lt1[$totLigands] = 
      $belowThr1 / $runsperligand[$i];
    $lt2[$totLigands] = 
      $belowThr2 / $runsperligand[$i];
    $totLigands++;
    $belowThr1 = 0;
    $belowThr2 = 0;
  }
  return (\@lt1, \@lt2);
}

  

sub stage
{
  my($nruns, $p, $i, $r);
  $nruns = $_[0];
  $p = $_[1];
  if ($nruns == 0)
  {
    return (1, 0);
  }
  for ($i = 0 ; $i < $nruns ; $i++)
  {
    $r = rand;
    if ($r < $p)
    {
      return (1, $i);
    }
  }
  return (0, $nruns);
}


