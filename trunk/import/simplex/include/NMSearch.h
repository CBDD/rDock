/*NMSearch.h
 *declarations of Nelder Mead Simplex Search
 *Adam Gurson College of William & Mary 1999
 */

#ifndef _NMSearch_
#define _NMSearch_

#include <vector>
using std::vector;
#include <iostream>
using std::istream;

#include "NMObjective.h"
#include "subscrpt.h"

class RbtChromElement;
class RbtBaseSF;

//default definitions: should be defined by user in objective.h, but just in case...
/*#ifndef maxCalls
#define maxCalls -1 // means program will only terminate based on stopping criteria
#endif*/
//#define stoppingStepLength 10e-8 /*sqrt(fabs(3.0 * (4.0/3.0 - 1.0) - 1.0))*/
#ifndef DEBUG
#define DEBUG 0
#endif
class RbtMatrix
{
  public:
    RbtMatrix(int M, int N, const double& value = 0.0);
    ~RbtMatrix();
    RbtMatrix(const RbtMatrix &A);
    vector<double>& row(int idx) const;
    vector<double>&  operator[](int i) const;
  private:
    vector<vector<double> > *m;
    int dimX, dimY;
};

class NMSearch
{ 
 public:
 // constructors and destructors
   NMSearch(RbtChromElement *cc,
            RbtBaseSF *psf, 
            double Alpha=1.0,
            double Beta=0.5,
            double Gamma=2.0,
            double Sigma=0.5);
   // constructor which allows coefficient initialization
 
   NMSearch(const NMSearch& Original);
   // deep copy constructor

   ~NMSearch();
   // destructor

 // algorithmic routines

   void ExploratoryMoves();
   // use Nelder Mead to find function minimum

   void ReplaceSimplexPoint(int index, const vector<double>& newPoint);
   // replaces simplex point indexed at index with newPoint

   void CalculateFunctionValue(int index);
   // finds the f(x) value for the simplex point indexed at index and
   // replaces the proper value in simplexValues

   void SetAlpha(double newAlpha);
   void SetBeta(double newBeta);
   void SetGamma(double newGamma);
   void SetSigma(double newSigma);
   static void SetMaxCalls(int mc);
   static void SetStoppingLength(double sl);
   // these functions allow the user to set values of the reflection,
   // contraction, expansion, and shrinking coefficients

   bool Stop();
   // returns true if the stopping criteria have been satisfied

   void fcnCall(int n, vector<double> *x, double& f, int& flag);
   // indirection of function call for purposes of keeping an accurate
   // tally of the number of function calls

 // Simplex-altering functions

   void InitRegularTriangularSimplex(const vector<double> *basePoint,
                                     const double edgeLength);
   // deletes any existing simplex and replaces it with a regular
   // triangular simplex in the following manner:
   //
   // basePoint points to a point that will be the "origin" of the
   //    simplex points (it will be a part of the simplex)
   // edgeLength is the length of each edge of the "triangle"
   //
   // functionCalls is reset to 0 and ALL FUNCTION VALUES ARE CALCULATED.
   //
   // NOTE: basePoint is assumed to be of proper dimension

   void InitFixedLengthRightSimplex(const vector<double> *basePoint,
                                    const double edgeLength);
   // deletes any existing simplex and replaces it with a right-angle
   // simplex in the following manner:
   //
   // basePoint points to a point that will be the "origin" of the
   //    simplex points (it will be a part of the simplex)
   // edgeLength is to be the length of each simplex side extending
   //    from the basePoint along each positive coordinate direction.
   //
   // functionCalls is reset to 0 and ALL FUNCTION VALUES ARE CALCULATED.
   //
   // NOTE: basePoint is assumed to be of proper dimension

   void InitVariableLengthRightSimplex(const vector<double> *basePoint,
                                       const double* edgeLengths);
   // deletes any existing simplex and replaces it with a right-angle
   // simplex in the following manner:
   //
   // basePoint points to a point that will be the "origin" of the
   //    simplex points (it will be a part of the simplex)
   // edgeLengths points to an array of n doubles, where n is the
   //    dimension of the given search. x_1 will then be located
   //    a distance of edgeLengths[0] away from the basepoint along the
   //    the x_1 axis, x_2 is edgeLengths[1] away on the x_2 axis, etc.
   //
   // functionCalls is reset to 0 and ALL FUNCTION VALUES ARE CALCULATED.
   //
   // NOTE: basePoint and edgeLengths are assumed to be of proper dimension

   void InitGeneralSimplex(const RbtMatrix *plex);
   // deletes any existing simplex and replaces it with the one
   // pointed to by plex
   //
   // functionCalls is reset to 0 and ALL FUNCTION VALUES ARE CALCULATED.
   //
   // NOTE: THIS ASSUMES THAT plex IS OF PROPER DIMENSION

   void ReadSimplexFile(istream& fp);
   // may also pass cin as input stream if desired
   // input the values of each trial point
   // NOTE: THIS FUNCTION WILL ONLY ACCEPT n+1 POINTS
   //
   // functionCalls is reset to 0 and ALL FUNCTION VALUES ARE CALCULATED.

 // Query functions

   int GetFunctionCalls() const;
   // number of objective function evaluations

   //DM 28 Jun 2002 - update this one function to use a more STL-friendly syntax
   //i.e. return a const STL vector<double> by reference, rather than a pointer by reference
   //void GetMinPoint(vector<double>* &minimum) const;
   const vector<double>& GetMinPoint() const;
   // simplex point which generates the best objective function
   // value found thus far
   // USER SHOULD PASS JUST A NULL POINTER, WITHOUT PREALLOCATED MEMORY

   double GetMinVal() const;
   // best objective function value found thus far

   void GetCurrentSimplex(RbtMatrix* &plex) const;
   // performs a deep copy of the simplex to a Matrix pointer
   // points to a newly allocated chunk of memory upon return
   // USER SHOULD PASS JUST A NULL POINTER, WITHOUT PREALLOCATED MEMORY

   void GetCurrentSimplexValues(double* &simValues) const;
   // performs a deep copy of the simplexValues array to a double pointer
   // points to a newly allocated chunk of memory upon return
   // USER SHOULD PASS JUST A NULL POINTER, WITHOUT PREALLOCATED MEMORY

   int GetVarNo() const;
   // returns the dimension of the problem

   int GetTolHit() const;
   // returns toleranceHit

   void printSimplex() const;
   // prints out the simplex points by row, their corresponding f(x)
   // values, and the number of function calls thus far

 private:

   void FindMinMaxIndices();
   // sets minIndex to the simplex index of the point which generates
   // the lowest value of f(x)
   // sets maxIndex to the simplex index of the point which generates
   // the highest value of f(x)

   int SecondHighestPtIndex();
   // returns simplex index of the point which
   // generates the second highest value of f(x)

   void FindCentroid();
   // finds the centroid

   void FindReflectionPt();
   // finds the reflection point and sets its f(x) value

   void FindExpansionPt();
   // finds the expansion point and sets its f(x) value

   void FindContractionPt();
   // finds the contraction point and sets its f(x) value

   void ShrinkSimplex();
   // this function goes through the simplex and reduces the
   // lengths of the edges adjacent to the best vertex



   double (*fun)(const vector<double> &, bool&); 
   int dimensions;                // the number of dimensions
                                  //    (the dimension of the problem)
   RbtMatrix *simplex;            // the current simplex
   double *simplexValues;         // their corresponding f(x) values
   double alpha;                  // reflection coefficient
   double beta;                   // contraction coefficient
   double gamma;                  // expansion coefficient
   double sigma;                  // shrinking coefficient
   int minIndex;                  // index of point generating min f(x)
   int maxIndex;                  // index of point generating max f(x)
   vector<double> *centroid;      // the current centroid
   vector<double> *reflectionPt;  // the reflection point
   double reflectionPtValue;      // the value of f(reflectionPt)
   vector<double> *expansionPt;   // the expansion point
   double expansionPtValue;       // the value of f(expansionPt)
   vector<double> *contractionPt; // the contraction point
   double contractionPtValue;     // the value of f(contractionPt)
   double maxPrimePtValue;        // min(f(maxIndexPoint),reflectionPtValue)
   long functionCalls;            // tally of number of function calls
   int toleranceHit;              // 1 if stop due to tolerance, 0 if funcCalls
   int maxPrimePtId;              // set by FindContractionPt() and used in
                                  // ExploratoryMoves() to branch in pos. 3

   // the following vectors are simply extra storage space
   // that are used by functions that require a vector of
   // size = dimensions
   vector<double> *scratch, *scratch2;
   RbtChromElement * c;
   RbtBaseSF *pSF;
   static int maxCalls;
   static double stoppingStepLength;

};


#endif
