//Random number generator classes. Taken from "The C++ Programming Language",
//Third Edition, Bjarne Stroustrup, Addison-Wesley, 1997,p686.

#ifndef _RANDINT_H_
#define _RANDINT_H_

//uniform distribution in the interval [0,max]
class Randint {
	unsigned long randx;
public:
	Randint(long s=0) {randx=s;}
	void seed(long s) {randx=s;}
	//DM 22 Feb 1999 - return current seed
	long GetSeed() {return randx;}
	//magic numbers chosen to use 31 bits of a 32-bit long:
	int abs(int x) {return x&0x7fffffff;}
	static double max() {return 2147483648.0;}
	int draw() {return randx=randx*1103515245+12345;}
	
	double fdraw() {return abs(draw())/max();}
	
	int operator() () {return abs(draw());}
};

//uniform distribution in the interval [0:n[
class Urand : public Randint {
	int n;
public:
	Urand(int nn) {n=nn;}

	int operator() () {int r=n*((int) fdraw()); return (r==n) ? n-1 : r;}
};
#endif //_RANDINT_H_
