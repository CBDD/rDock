package SDRecord;
# Perl module for reading/writing SD records from/to STDIN/STDOUT
# Methods:  new
#           readRec
#           writeRec
#           writeData
#           copy
#           addData
#
# Data:     LINES    - reference to array of record text lines
#           DATA     - reference to hash array of record data fields
#           DATAREF  - reference to hash array of line numbers for data fields

my $LINES = 'LINES';
my $DATA = 'DATA';
#DM 27 Sept 1999 - also store line number of beginning each data field
#so we can extract multi-line fields more easily from the LINES array
my $DATAREF = 'DATAREF';

############################################################
# Constructor
sub new {
  my $this = {};
  bless $this;
  return $this;
}

############################################################
# readRec() - read next SD record from STDIN
# Input params:
#   'LINES' => any value - save lines
#   'DATA'  => any value - save data
sub readRec {
  my $this = shift;
  my %params = @_;#input parameters
  my @lines;#array for storing all text lines
  my %data;#hash array for storing all data key,value pairs
  my %dataref;#hash array for storing line number references for each data field
  my ($fieldName,$ob,$cb);

  #clear current data
  delete $this->{$LINES};
  delete $this->{$DATA};
  delete $this->{$DATAREF};

  #read lines from STDIN to the next record delimiter
  #store in @lines array
  #DM 16 Nov 1999 - can now cope with DOS text files directly
  $/ = "\n";#Default Unix line separator
  while (<>) {
    chomp;#Get rid of Unix line separator
    $/ = "\r";
    chomp;#Get rid of DOS line separator
    $/ = "\n";#Reset to normal Unix
    last if ($_ eq '$$$$');#end of record
    push @lines,$_;
    #print "$_\n";
  }

  #check if we read anything
  if (scalar(@lines) > 0) {
    #store ref to @lines if required
    $this->{$LINES} = \@lines if (defined $params{$LINES});
    #search for data fields if required
    if (defined $params{$DATA}) {
      $fieldName = '';
      my $lineNum = 0;
      foreach $line (@lines) {
	$lineNum++;
	#DM 12 Jul 1999 - include the first three title lines as
	#pseudo data fields
	if ($lineNum <= 3) {
	  $data{"_TITLE$lineNum"} = $line;
	  #DM 20 Dec 1999 - include dimensionality (2D/3D) as pseudo data field
	  if ($lineNum == 2) {
	    $data{"_NDIM"} = substr($line,20,1);
	  }
	}
	#DM 05 Aug 1999 - include number of atoms as pseudo data field
	elsif ($lineNum == 4) {
	  my @fields = split(" ",$line);
	  $data{"_NATOMS"} = $fields[0] if (scalar(@fields)>0);
	}
	if (index($line,'>') == 0) {#found a data field
	  $ob = index($line,'<');#first open bracket
	  $cb = rindex($line,'>');#last close bracket
	  if (($ob != -1) && ($cb != -1)) { #save field name
	    $fieldName = substr($line,$ob+1,$cb-$ob-1);
	    $dataref{$fieldName} = $lineNum-1;#store ref to line number
	  }
	}
	#if field name defined, then store first line of value
	elsif ($fieldName ne '') {
	  $data{$fieldName} = $line;
	  $fieldName = '';#clear field name
	}
      }
      #store ref to %data in $this
      $this->{$DATA} = \%data;
      $this->{$DATAREF} = \%dataref;
    }
    return 1;
  }
  else {
    return 0;
  }
}

############################################################
# writeRec() - write current record to STDOUT (mol + data)
sub writeRec {
  my $this = shift;
  if (defined $this->{$LINES}) {
    foreach $line (@{$this->{$LINES}}) {
      print "$line\n";
    }
    print "\$\$\$\$\n";
  }
}

############################################################
# writeMol() - write current mol record to STDOUT
sub writeMol {
  my $this = shift;
  if (defined $this->{$LINES}) {
    foreach $line (@{$this->{$LINES}}) {
      print "$line\n";
      last if ($line eq "M  END");
    }
  }
}

############################################################
# writeData() - list data field/values to STDOUT
sub writeData {
  my $this = shift;
  my ($keys ,$value, $lineNum );
  if (defined $this->{$DATA}) {
    foreach $key (sort keys %{$this->{$DATA}}) {
      $value = $this->{$DATA}->{$key};
      $lineNum = $this->{$DATAREF}->{$key};
      #print "$key eq $value (line $lineNum)\n";
      print "\$$key eq \"$value\"\n";
    }
  }
}

############################################################
# copy() - create deep copy of SDRecord
# Input params:
#   'LINES' = any value - deep copy of lines
#   'DATA'  = any value - deep copy of data and dataref

sub copy {
  my $this = shift;
  my %params = @_;#input parameters
  my $clone = new SDRecord;
  if ((defined $params{$LINES}) && (defined $this->{$LINES})) {
    my @lines = @{$this->{$LINES}};
    $clone->{$LINES} = \@lines;
  }
  if (defined $params{$DATA}) {
    if (defined $this->{$DATA}) {
      my %data = %{$this->{$DATA}};
      $clone->{$DATA} = \%data;
    }
    # DM 24 Jul 2001 - DATAREF was missing from deep copy
    # We copy the dataref array using the same DATA input parameter
    # i.e. if the data is copied then so is dataref
    if (defined $this->{$DATAREF}) {
      my %dataref = %{$this->{$DATAREF}};
      $clone->{$DATAREF} = \%dataref;
    }
  }
  return $clone;
}

############################################################
# addData() - adds data to data hash array
sub addData {
  my $this = shift;
  my %params = @_;#input parameters
  #if the array is already defined then add to it
  if (defined $this->{$DATA}) {
    my ($keys ,$value );
    while (($key, $value) = each (%params)) {
      $this->{$DATA}->{$key} = $value;
    }
  }
  #if not defined, then create it
  else {
    $this->{$DATA} = \%params;
  }
}

############################################################
# addDataAndLines() - adds data to data hash array,
# and adds corresponding lines also so that record may be
# rewritten with the new fields
sub addDataAndLines {
  my $this = shift;
  my %params = @_;#input parameters

  while (($key, $value) = each (%params)) {
    #DM 23 Oct 2000
    #Check if data field already exists, if so replace value
    #Note: this only works for single line data fields
    if (defined $this->{$DATA}->{$key}) {
      my $keyRef = $this->{$DATAREF}->{$key};
      ${$this->{$LINES}}[$keyRef+1] = $value;
    }
    #else add the corresponding lines to the lines arrays
    else {
      my $keyRef = scalar(@{$this->{$LINES}});
      $this->{$DATAREF}->{$key} = $keyRef;
      push @{$this->{$LINES}},">  <$key>";
      push @{$this->{$LINES}},"$value";
      push @{$this->{$LINES}},"";
    }
  }

  #DM 23 Oct 2000 - add data after adding the lines
  $this->addData(%params);

}

1;
