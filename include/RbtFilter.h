/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtFilter.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Reads a file that contains a list of filters that determine
//when to terminate with a given ligand, and another filter to determine
//when to write the current conformation to the output.

#ifndef _RBT_FILTER_H_
#define _RBT_FILTER_H_

#include "RbtBaseObject.h"
#include "RbtContext.h"
#include "RbtFilterExpression.h"


class RbtFilterExpressionVisitor;

class RbtFilter : public RbtBaseObject
{
public:
  static RbtString _CT;
  RbtFilter(RbtString strfilter, RbtBool filter = false);
    ///////////////////
    // Destructor
    //////////////////
  virtual ~RbtFilter();
  //Override RbtObserver pure virtual
  //Notify observer that subject has changed
  virtual void Update(RbtSubject* theChangedSubject);
  
  void SetupReceptor();  //Called by Update when receptor is changed
  void SetupLigand();  //Called by Update when ligand is changed
  void SetupScore();  //Called by Update when either model has changed
  RbtBool Write();    //Output conformation?
  RbtBool Terminate();  //Finished with ligand?
  RbtModelPtr GetReceptor() const;
  RbtModelPtr GetLigand() const;
  void SetMaxNRuns(RbtInt n){maxnruns = n;}

  
    ////////////////////
    // Private methods
    ////////////////////

protected:
  RbtFilter();  // Default constructor disabled

private:
  void ReadFilters();

  RbtInt filteridx, nTermFilters, nWriteFilters, nruns, maxnruns;
  RbtFilterExpressionList terminationFilters;
  RbtFilterExpressionList writtingFilter;
  RbtModelPtr m_spReceptor;
  RbtModelPtr m_spLigand;
  RbtContextPtr contextp;
};

//Useful typedefs
typedef SmartPtr<RbtFilter> RbtFilterPtr;  //Smart pointer
typedef vector<RbtFilterPtr> RbtFilterList;
                    //Vector of smart pointers
typedef RbtFilterList::iterator RbtFilterListIter;
typedef RbtFilterList::const_iterator 
                                 RbtFilterListConstIter;

#endif //_RbtFilter_H_
