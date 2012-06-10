#!/usr/bin/perl
# Misc Perforce-related scripts to dynamically detect
# product, codeline and version
# Must be run in situ in the build directory
#
# e.g. if build path is ~/dev/rdock/enspiral_dev/build
#
# PRODUCT  = rdock
# CODELINE = enspiral_dev
# VERSION  = enspiral_dev.0
#
# Only difference between CODELINE and VERSION is that
# if CODELINE is non-numeric (e.g. enspiral_dev)
# then VERSION has a ".0" appended, to keep tmake happy
# with shared library names
#
# BASEVER is the lowest changelist for this CODELINE
# BUILDNUM is highest changelist for this CODELINE
#########################################################
my $arg = shift @ARGV;
if ($arg eq "PRODUCT") {
  print &product;
}
elsif ($arg eq "CODELINE") {
  print &codeline;
}
elsif ($arg eq "VERSION") {
  my $version = &codeline;
  chomp $version;
  $version .= ".0" if (index($version,".") == -1);
  print $version;
}
elsif ($arg eq "BUILDNUM") {
  print &buildnum;
}
elsif ($arg eq "BASEVER") {
  print &basever;
}
elsif ($arg eq "TEST") {
  if (&p4installed) {
    print "p4 is installed\n";
  }
  else {
    print "p4 is not installed\n";
  }
}

#No longer requires p4
sub product {
  my $pwd = `pwd`;
  chomp $pwd;
  my @pathElements = split("/", $pwd);
  my $n = scalar(@pathElements);
  return ($n > 2) ? $pathElements[$n-3] : "rdock";
}

#No longer requires p4
sub codeline {
  my $pwd = `pwd`;
  chomp $pwd;
  my @pathElements = split("/", $pwd);
  my $n = scalar(@pathElements);
  return ($n > 1) ? $pathElements[$n-2] : "2006.1";
}

#Requires p4, returns dummy value otherwise
sub buildnum {
  return &p4installed ? `p4 changes -i -m 1 ../...\#have | cut -f2,4 -d" "` : '902';
}
#Requires p4, returns dummy value otherwise
sub basever {
  return &p4installed ? `p4 changes ../... | tail -1|cut -d" " -f2` : 'UNKNOWN';
}

sub p4installed {
  my $p4exe = `which p4`;
  chomp $p4exe;
  return (-x $p4exe)
}
