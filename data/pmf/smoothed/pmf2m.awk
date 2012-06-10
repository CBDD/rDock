#!/bin/awk -f
{
	a[NR]=$1;
	b[NR]=$2;
}END{
	print "# name: x"
	print "# type: matrix"
	print "# rows: 1"
	print "# columns: "NR
	for(i=1;i<=NR;i++) 
		printf("%f ",a[i]);
	print ""
	print "# name: y"
	print "# type: matrix"
	print "# rows: 1"
	print "# columns: "NR
	for(i=1;i<=NR;i++) 
		printf("%f ",b[i]); 
	print " ";
}
