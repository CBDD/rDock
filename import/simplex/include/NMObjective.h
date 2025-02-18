//objective.h
//this example file declares the functions and variable n as required
//by the search implementations
//simple definitions for certain parameters may be given here as well.

#if !defined _userfile_
#define _userfile_

#include <stdlib.h>
#include <vector>

double getInitialStep();
//void fcn(double (*fun) (const vector<double>& , bool&), int , double *, double & , int & );

void initpt(int , double *&);

extern int n;
//extern double tolerance;
extern bool newInitialPoint;
extern char Fname[120];

/* Initial trial step length, reflecting the refinement of the search lattice
 */
#define initialStepLength 2.0 
/*getInitialStep()*/

/*  The search will halt if the trial step length is less than 
 *  or equal to the stoppingStepLength after a reduction of the 
 *  trial step length.
 */ 
//#define stoppingStepLength tolerance

/*  The maxCalls variable is not a hard limit on the number of function
 *  evaluations that may be performed by the search.  The search should
 *  always halt at a point from which none of the trial steps can find 
 *  improvement to be consistent with the formal theory regarding pattern
 *  searches.
 */
extern long maxCalls;

#endif









