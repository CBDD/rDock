#!/bin/awk -f
BEGIN{
	true			= 1;
	false			= !true;
	hashmarkcount	= 0;
	iscomment		= true;	# true for the first lines
	i=0;
	j=0;
}{
	if("#" == $1) {
		hashmarkcount++;
		iscomment	= true;
	} else 
		iscomment	= false;
	if(hashmarkcount < 6 && !iscomment ){	#get the first vector (x)
		x[i]=$1;
		i++;
	}
	else if ( !iscomment ){					# so we are at y
		y[j]=$1;
		j++;
	}
		
}END{
	# first padding at the beginning
	for (j=0.2;j<x[0];j+=0.2) {
		printf("%6.3f  %6.4f\n", j,3.0);
	}
	# note i sould be equal to y 
	# we are considering 'i' number of data (
	# ys without x values will be discarded, xs without y will
	# get y value 0
	# we are treating 'i' as the length of the vectors
	for (j=0;j<i;j++) {
		printf("%6.3f  %6.4f\n", x[j],y[j]);
	}
	# zero padding at the end
	for(j=x[j-1]+0.2;j<=12.0;j+=0.2) {
		printf("%6.3f  %6.4f\n", j,0.0);
	}
}
