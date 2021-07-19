/*NMSearch.cc
 *class implementation of Nelder Mead Simplex Search
 *Adam Gurson College of William & Mary 1999
 *
 * modified slightly by Anne Shepherd (pls), 8/00
 * modified by Beatriz Garmendia-Doval, 12/01
 * RbtMatrix destructor added, Dave Morley, 06/02
 */

#include "NMSearch.h"
#include "RbtChromElement.h"
#include "RbtBaseSF.h"

int NMSearch::maxCalls = -1;
double NMSearch::stoppingStepLength = 1e-8;

RbtMatrix::RbtMatrix(int M, int N, const double& value/* = 0.0*/)
{
   vector<double> aux(N, 0.0);
   m = new vector<vector<double> >; 
   for (int i = 0 ; i < M ; i++)
      m->push_back(aux);
   dimX = M;
   dimY = N;
}

RbtMatrix::~RbtMatrix() {
  delete m;
}

RbtMatrix::RbtMatrix(const RbtMatrix &A)
{
   m = new vector<vector<double> >(*(A.m));
   dimX = A.dimX;
   dimY = A.dimY;
}

vector<double>& RbtMatrix::row(int idx) const
{
   return (*m)[idx];
}

vector<double>& RbtMatrix::operator[](int idx) const
{
   return (*m)[idx];
}

NMSearch::NMSearch(RbtChromElement *cc, RbtBaseSF *psf, 
                   double Alpha, double Beta,
                   double Gamma, double Sigma)
{
   c = cc;
   pSF = psf;
   dimensions = c->GetLength();
   functionCalls = 0;
   simplex = NULL;
   simplexValues = NULL;
   centroid = new vector<double>(dimensions,0.0);
   reflectionPt = new vector<double>(dimensions,0.0);
   expansionPt = new vector<double>(dimensions,0.0);
   contractionPt = new vector<double>(dimensions,0.0);
   alpha = Alpha;
   beta = Beta;
   gamma = Gamma;
   sigma = Sigma;
   scratch = new vector<double>(dimensions,0.0);
   scratch2 = new vector<double>(dimensions,0.0);
} // NMSearch() (special)

NMSearch::NMSearch(const NMSearch& Original)
{
   c = Original.c;
   pSF = Original.pSF;
   dimensions = Original.GetVarNo();
   if (Original.simplex != NULL)
     Original.GetCurrentSimplex(simplex);
   else
     simplex = NULL;
   if (Original.simplexValues != NULL)
     Original.GetCurrentSimplexValues(simplexValues);
   else
     simplexValues = NULL;
   alpha = Original.alpha;
   beta = Original.beta;
   gamma = Original.gamma;
   sigma = Original.sigma;
   minIndex = Original.minIndex;
   maxIndex = Original.maxIndex;
   centroid = new vector<double>(*(Original.centroid));
   reflectionPt = new vector<double>(*(Original.reflectionPt));
   reflectionPtValue = Original.reflectionPtValue;
   expansionPt = new vector<double>(*(Original.expansionPt));
   expansionPtValue = Original.expansionPtValue;
   contractionPt = new vector<double>(*(Original.contractionPt));
   contractionPtValue = Original.contractionPtValue;
   functionCalls = Original.functionCalls;
} // NMSearch() (copy constructor)

NMSearch::~NMSearch()
{
   if(simplex != NULL) delete simplex;
   if(simplexValues != NULL) delete [] simplexValues;
   delete centroid;
   delete reflectionPt;
   delete expansionPt;
   delete contractionPt;
   delete scratch;
   delete scratch2;
   //NOTE: Matrix and vector classes have their own destructors
} // ~NMSearch

// algorithmic routines

void NMSearch::ExploratoryMoves()
{
   double secondHighestPtValue; // used for contraction/reflection decision
   toleranceHit = 0;

   FindMinMaxIndices();
   do {
      if(DEBUG) printSimplex();
      FindCentroid();
      secondHighestPtValue = simplexValues[SecondHighestPtIndex()];
    // reflection step
      FindReflectionPt();

      // stop if at maximum function calls and update the simplex
      /*changed  8/8/00 to fix the problem of maxCalls == -1  --pls
        formerly read if(functionCalls <= maxCalls)
      */
      if (maxCalls > -1
	  && functionCalls >= maxCalls) {
	FindMinMaxIndices();
	ReplaceSimplexPoint(maxIndex, *reflectionPt);
        if(DEBUG) cout << "reflection\n";
	simplexValues[maxIndex] = reflectionPtValue;
	FindMinMaxIndices(); 
	return;
      } // if using call budget
      
      // possibility 1
      if(simplexValues[minIndex] > reflectionPtValue) {
	FindExpansionPt(); // expansion step
	
	if (reflectionPtValue > expansionPtValue) {
	  ReplaceSimplexPoint(maxIndex, *expansionPt);
          if(DEBUG) cout << "expansion\n";
	  simplexValues[maxIndex] = expansionPtValue;
	} // inner if
	else {
	  ReplaceSimplexPoint(maxIndex, *reflectionPt);
          if(DEBUG) cout << "reflection\n";
	  simplexValues[maxIndex] = reflectionPtValue;
	} // else         
      } // if for possibility 1
      
    // possibility 2

      else if( (secondHighestPtValue > reflectionPtValue        ) &&
               (   reflectionPtValue >= simplexValues[minIndex]) ) {
         ReplaceSimplexPoint(maxIndex, *reflectionPt);
         if(DEBUG) cout << "reflection\n";
         simplexValues[maxIndex] = reflectionPtValue;
      } // else if for possibility 2

    // possibility 3
      else if( reflectionPtValue >= secondHighestPtValue ) {
         FindContractionPt(); // contraction step
         if(maxPrimePtId == 0) {
           if( contractionPtValue > maxPrimePtValue ) {
             ShrinkSimplex();
           } // inner if
           else {
             ReplaceSimplexPoint(maxIndex, *contractionPt);
             if(DEBUG) cout << "contraction\n";
             simplexValues[maxIndex] = contractionPtValue;
           } // inner else
         } // maxPrimePtId == 0
         else if(maxPrimePtId == 1) {
           if( contractionPtValue >= maxPrimePtValue ) {
             ShrinkSimplex();
           } // inner if
           else {
             ReplaceSimplexPoint(maxIndex, *contractionPt);
             if(DEBUG) cout << "contraction\n";
             simplexValues[maxIndex] = contractionPtValue;
           } // inner else
         } // maxPrimePtId == 1
      } // else if for possibility 3

    // if we haven't taken care of the current simplex, something's wrong
      else {
         cerr << "Error in ExploratoryMoves() - "
              << "Unaccounted for case.\nTerminating.\n";
         return;
      }
   FindMinMaxIndices();
   } while (!Stop());   // while stopping criteria is not satisfied
} // ExploratoryMoves()

void NMSearch::ReplaceSimplexPoint(int index, const vector<double>& newPoint)
{
   for( int i = 0; i < dimensions; i++ ) {
      (*simplex)[index][i] = newPoint[i];
   } // for
} // ReplaceSimplexPoint()

void NMSearch::CalculateFunctionValue(int index)
{
   for (int i = 0 ; i < scratch->size() ; i++)
     (*scratch)[i] = (*simplex).row(index)[i];
   int success;
   fcnCall(dimensions, scratch, simplexValues[index], success);
   if(!success) cerr<<"Error calculating point in CalculateFunctionValue().\n";
} // CalculateFunctionValue()

void NMSearch::SetAlpha(double newAlpha)
{
   alpha = newAlpha;
} // SetAlpha()

void NMSearch::SetBeta(double newBeta)
{
   beta = newBeta; 
} // SetBeta()

void NMSearch::SetGamma(double newGamma)
{
   gamma = newGamma;
} // SetGamma()

void NMSearch::SetSigma(double newSigma)
{
   sigma = newSigma;
} // SetGamma()

void NMSearch::SetMaxCalls(int mc)
{
  maxCalls = mc;
}  // SetMaxCalls

void NMSearch::SetStoppingLength(double sl)
{
  stoppingStepLength = sl;
}  // SetStoppingLength

bool NMSearch::Stop()
{
   if(maxCalls > -1) {
      if(functionCalls >= maxCalls)
         return true;
   }

   double mean = 0.0;

   for( int i = 0; i <= dimensions; i++) {
      if( i != minIndex ) {
         mean += simplexValues[i];
      } // if
   } //for 

   mean /= (double)dimensions;

   // Test for the suggested Nelder-Mead stopping criteria
   double total = 0.0;
   for( int i = 0; i <= dimensions; i++ ) {
      total += pow( simplexValues[i] - mean ,2);
   } //for
   total /= ((double)dimensions + 1.0);
   total = sqrt(total);
   
   
   // printSimplex();
   if(total < stoppingStepLength) {
      toleranceHit = 1;
      return true;
   }
   else
      return false;
} // Stop()

void NMSearch::fcnCall(int n, vector<double> *x, double& f, int& flag)
{
    //TODO: Surround SetVector with try/catch
    c->SetVector(*x);
    c->SyncToModel();
    f = pSF->Score();
    flag = true;
    functionCalls++;
    //TODO: Check if GetVector is necessary. Some of the element
    //values may have been standardised inside SetVector
    //e.g. angles that have gone out of range.
    //By calling GetVector here we may change the value of *x
    //Does this impact the search?
    //c->GetVector(*x);
}

// Simplex-altering functions

void NMSearch::InitRegularTriangularSimplex(const vector<double> *basePoint,
                                            const double edgeLength)
{
  //  This routine constructs a regular simplex (i.e., one in which all of 
  //  the edges are of equal length) following an algorithm given by Jacoby,
  //  Kowalik, and Pizzo in "Iterative Methods for Nonlinear Optimization 
  //  Problems," Prentice-Hall (1972).  This algorithm also appears in 
  //  Spendley, Hext, and Himsworth, "Sequential Application of Simplex 
  //  Designs in Optimisation and Evolutionary Operation," Technometrics, 
  //  Vol. 4, No. 4, November 1962, pages 441--461.

   int i,j;
   double p, q, temp;
   RbtMatrix *plex = new RbtMatrix(dimensions+1,dimensions,0.0);
   for( int col = 0; col < dimensions; col++ ) {
      (*plex)[0][col] = (*basePoint)[col];
   }

   temp = dimensions + 1.0;
   q = ((sqrt(temp) - 1.0) / (dimensions * sqrt(2.0))) * edgeLength;
   p = q + ((1.0 / sqrt(2.0)) * edgeLength);

   for(i = 1; i <= dimensions; i++) { 
      for(j = 0; j <= i-2; j++) {
         (*plex)[i][j] = (*plex)[0][j] + q;
      } // inner for 1
      j = i - 1;
      (*plex)[i][j] = (*plex)[0][j] + p;
      for(j = i; j < dimensions; j++) {
            (*plex)[i][j] = (*plex)[0][j] + q;
      } // inner for 2
   } // outer for

   InitGeneralSimplex(plex);
   delete plex;
} // InitRegularTriangularSimplex()

void NMSearch::InitFixedLengthRightSimplex(const vector<double> *basePoint,
                                           const double edgeLength)
{
  // to take advantage of code reuse, this function simply turns
  // edgeLength into a vector of dimensions length, and then
  // calls InitVariableLengthRightSimplex()
   double* edgeLengths = new double[dimensions];
   for( int i = 0; i < dimensions; i++ ) {
      edgeLengths[i] = edgeLength;
   }
   InitVariableLengthRightSimplex(basePoint,edgeLengths);
   delete [] edgeLengths;
} // InitFixedLengthRightSimplex()

void NMSearch::InitVariableLengthRightSimplex(const vector<double> *basePoint,
                                              const double* edgeLengths)
{
   RbtMatrix *plex = new RbtMatrix(dimensions+1,dimensions,0.0);
   for( int i = 0; i < dimensions; i++ ) {
      // we're building the basePoint component-by-component into
      // the (n+1)st row
      (*plex)[dimensions][i] = (*basePoint)[i];

      // now fill in the ith row with the proper point
      for( int j = 0; j < dimensions; j++ ) {
         ((*plex)[i])[j] = (*basePoint)[j];
         if( i == j )
            (*plex)[i][j] += edgeLengths[i];
      }
   } // for
   InitGeneralSimplex(plex);
   delete plex;
} // InitVariableLengthRightSimplex()

void NMSearch::InitGeneralSimplex(const RbtMatrix *plex)
{
   functionCalls = 0;
   if( simplex != NULL ) { delete simplex; }
   if( simplexValues != NULL ) { delete [] simplexValues;}
   simplex = new RbtMatrix((*plex));
   simplexValues = new double[dimensions+1];

   int success;
   for( int i = 0; i <= dimensions; i++ ) {
      for (int j = 0 ; j < scratch->size() ; j++)
         (*scratch)[j] = (*plex).row(i)[j];
      fcnCall(dimensions, scratch, simplexValues[i], success);
      if(!success) cerr<<"Error with point #"<<i<<" in initial simplex.\n";
   } // for
   FindMinMaxIndices();
} // InitGeneralSimplex()

void NMSearch::ReadSimplexFile(istream& fp)
{
   if(!fp) {
      cerr<<"No Input Stream in ReadSimplexFile()!\n";
      return; // There's no file handle!!
   }

   RbtMatrix *plex = new RbtMatrix(dimensions+1,dimensions);
   for( int i = 0; i <= dimensions; i++ ) {
      for ( int j = 0; j < dimensions; j++ ) {
         fp >> (*plex)[i][j];
      } // inner for
   } // outer for
   InitGeneralSimplex(plex);
   delete plex;
} // ReadSimplexFile()

// Query functions

int NMSearch::GetFunctionCalls() const
{
   return functionCalls;
} // GetFunctionCalls()

//DM 28 Jun 2002 - update this one function to use a more STL-friendly syntax
//i.e. Return a const STL vector<double> by reference, rather than a pointer by reference
//void NMSearch::GetMinPoint(vector<double>* &minimum) const;
const vector<double>& NMSearch::GetMinPoint() const
{
  //minimum = new vector<double>((*simplex).row(minIndex).size());
  //for (int i = 0 ; i < minimum->size() ; i++)
  //  (*minimum)[i] = (*simplex).row(minIndex)[i];
  return (*simplex).row(minIndex);
} // GetMinPoint()

double NMSearch::GetMinVal() const
{
   return simplexValues[minIndex];
} // GetMinVal()

void NMSearch::GetCurrentSimplex(RbtMatrix* &plex) const
{
   plex = new RbtMatrix((*simplex));
} // GetCurrentSimplex()

void NMSearch::GetCurrentSimplexValues(double* &simValues) const
{
   simValues = new double[dimensions+1];
   for( int i = 0; i <= dimensions; i++ ) {
      simValues[i] = simplexValues[i];
   } // for
} // GetCurrentSimplexValues()

int NMSearch::GetVarNo() const
{
   return dimensions;
} // GetVarNo()

int NMSearch::GetTolHit() const
{
   return toleranceHit;
} // GetTolHit()

// private functions

void NMSearch::FindMinMaxIndices()
{
   if(simplexValues == NULL) {
      cerr << "Error in FindMinMaxIndices() - "
           << "The vector of simplexValues is NULL!!\n";
      return;
   }
   minIndex = 0;
   maxIndex = dimensions;
   double min = simplexValues[0];
   double max = simplexValues[dimensions];
   for( int i = 1; i <= dimensions; i++ ) {
      if( simplexValues[i] < min ) {
         min = simplexValues[i];
         minIndex = i;
      } // if
      if( simplexValues[dimensions-i] > max ) {
         max = simplexValues[dimensions-i];
         maxIndex = dimensions - i;
      } // if
   } // for
} // FindMinMaxIndices()

int NMSearch::SecondHighestPtIndex()
{
   if(simplexValues == NULL) {
      cerr << "Error in SecondHighestPtValue() - "
           << "The vector of simplexValues is NULL!!\n";
      return -1;
   }
   int secondMaxIndex = minIndex;
   double secondMax = simplexValues[minIndex];
   for( int i = 0; i <= dimensions; i++ ) {
      if(i != maxIndex) {
         if( simplexValues[i] > secondMax ) {
            secondMax = simplexValues[i];
            secondMaxIndex = i;
         } // inner if
      } // outer if
   } // for
   return secondMaxIndex;
} // SecondHighestPtValue()

void NMSearch::FindCentroid()
{
   for (int i = 0 ; i < centroid->size() ; i++)
     (*centroid)[i] = 0.0;
   for( int i = 0; i <= dimensions; i++ ) {
      if( i != maxIndex ) {
         for (int j = 0 ; j < centroid->size() ; j++)
         (*centroid)[j] = (*centroid)[j] + (*simplex).row(i)[j];
      } // if
   } // for
   for (int i = 0 ; i < centroid->size() ; i++)
     (*centroid)[i] = (*centroid)[i] * ( 1.0 / (double)dimensions );
} // FindCentroid()

void NMSearch::FindReflectionPt()
{ 
   for (int i = 0 ; i < reflectionPt->size() ; i++)
     (*reflectionPt)[i] = 0.0;
   for (int i = 0 ; i < centroid->size() ; i++)
     (*reflectionPt)[i] = ( (*centroid)[i] * (1.0 + alpha) ) -
                     ( alpha * (*simplex).row(maxIndex)[i] );
   int success;
   fcnCall(dimensions, reflectionPt, reflectionPtValue, success);
   if(!success) {
      cerr << "Error finding f(x) for reflection point at"
           << "function call #" << functionCalls << ".\n";
   } // if
} // FindReflectionPt()

void NMSearch::FindExpansionPt()
{
   for (int i = 0 ; i < expansionPt->size() ; i++)
     (*expansionPt)[i] = 0.0;
   for (int i = 0 ; i < centroid->size() ; i++)
     (*expansionPt)[i] = ( (*centroid)[i] * (1.0 - gamma) ) +
                    ( gamma * (*reflectionPt)[i] );
   int success;
   fcnCall(dimensions, expansionPt, expansionPtValue, success);
   if(!success) {
      cerr << "Error finding f(x) for expansion point at"
           << "function call #" << functionCalls << ".\n";
   } // if
} // FindExpansionPt()

void NMSearch::FindContractionPt()
{
  // need to first define maxPrimePt
  //DM 27 Jun 2002 - we can't simply delete maxPrimePt at
  //the end of function because if the 'else' clause is
  //executed we would actually delete reflectionPt instead
  //So we use tmpPt instead to actually create a temp vector
  //inside the 'if' clause, which we can safely delete at the end
  vector<double> *tmpPt(NULL);
  vector<double> *maxPrimePt(NULL);
  if(simplexValues[maxIndex] <= reflectionPtValue) {
    tmpPt = new vector<double>(*scratch); //->size());
    maxPrimePt = tmpPt;
    for (int i = 0 ; i < scratch->size() ; i++)
         (*maxPrimePt)[i] = (*simplex).row(maxIndex)[i];
      maxPrimePtValue = simplexValues[maxIndex];
      maxPrimePtId = 1;
   } // if
   else {
      maxPrimePt = reflectionPt;
      maxPrimePtValue = reflectionPtValue;
      maxPrimePtId = 0;
   } // else

   for (int i = 0 ; i < centroid->size() ; i++)
     (*contractionPt)[i] = ( (*centroid)[i] * (1.0 - beta) ) +
                      ( beta * (*maxPrimePt)[i] );

   int success;
   fcnCall(dimensions, contractionPt, contractionPtValue, success);
   if(!success) {
      cerr << "Error finding f(x) for contraction point at"
           << "function call #" << functionCalls << ".\n";
   } // if
   delete tmpPt;//DM 27 Jun 2002 garbage collection
} // FindContractionPt()

void NMSearch::ShrinkSimplex()
{
   // stop if at maximum function calls
  // changed 5/01 to reflect maxcalls = -1 possibility ---pls
  if ( (maxCalls != (-1)) 
       && (functionCalls >= maxCalls) ) {return;}

   vector<double> *lowestPt= scratch;
   for (int i = 0 ; i < scratch->size() ; i++)
      (*lowestPt)[i] = (*simplex).row(minIndex)[i];
   vector<double> *tempPt = scratch2;
   int success;
   for( int i = 0; i <= dimensions; i++ ) {
      if( i != minIndex ) {
         for( int j = 0; j < dimensions; j++ ) {
            (*tempPt)[j] = (*simplex).row(i)[j];
            (*tempPt)[j] = (*tempPt)[j] + 
                           ( sigma * ( (*lowestPt)[j]-(*tempPt)[j] ) );
            (*simplex)[i][j] = (*tempPt)[j];
         } // inner for
         fcnCall(dimensions,tempPt,simplexValues[i],success);
         if (!success) cerr << "Error shrinking the simplex.\n";
         
         // stop if at maximum function calls 
	 // changed 5/01 to reflect maxcalls = -1 possibility ---pls
	 if ( (maxCalls != (-1)) 
	      && (functionCalls >= maxCalls) ) {return;}
	 
      } // if
   } // outer for
} // ShrinkSimplex()

void NMSearch::printSimplex() const
{
  for( int i = 0; i <= dimensions; i++ ) {
     cout << "   Point:";
     for ( int j = 0; j < dimensions; j++ ) {
        cout << (*simplex)[i][j] << "\t";
     } // inner for
     cout << "Value:" << simplexValues[i] << "\n";
  } // outer for

  cout << "\nFCalls: " << functionCalls << endl << endl;
}
