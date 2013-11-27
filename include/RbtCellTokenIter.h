/***********************************************************************
* The rDock program was developed from 1998 - 2006 by the software team 
* at RiboTargets (subsequently Vernalis (R&D) Ltd).
* In 2006, the software was licensed to the University of York for 
* maintenance and distribution.
* In 2012, Vernalis and the University of York agreed to release the 
* program as Open Source software.
* This version is licensed under GNU-LGPL version 3.0 with support from
* the University of Barcelona.
* http://rdock.sourceforge.net/
***********************************************************************/

//Defines the interface for a token iterator to be used by the parser

#ifndef _RBT_CELLTOKENITER_H_
#define _RBT_CELLTOKENITER_H_

#include "RbtConfig.h"
#include "RbtToken.h"
#include "RbtContext.h"
#include "RbtGPChromosome.h"
#include <stack>
using std::stack;

#include "RbtTokenIter.h"

class RbtCellTokenIter : public RbtTokenIter
{
public:

  static RbtString _CT;
    ///////////////////
    // Constructors
    ///////////////////
  RbtCellTokenIter(const RbtCellTokenIter&);
  RbtCellTokenIter(const RbtGPChromosomePtr, RbtContextPtr);
 
    ///////////////////
    // Destructor
    //////////////////
  virtual ~RbtCellTokenIter();
  
  virtual void copy(const RbtCellTokenIter &);
  //void First();
  void Next(RbtContextPtr);
  //RbtBool IsDone() const;
  RbtTokenPtr Current() ;
 
    ////////////////////
    // Private methods
    ////////////////////

private:
  RbtCellTokenIter();  // Default constructor disabled
  RbtGPChromosomePtr chrom;
  RbtTokenPtr current;
  stack<RbtInt> cells;
  RbtContextPtr contextp;
  
};

//Useful typedefs
typedef SmartPtr<RbtCellTokenIter> RbtCellTokenIterPtr;  //Smart pointer
typedef vector<RbtCellTokenIterPtr> RbtCellTokenIterList;//Vector of smart pointers
typedef RbtCellTokenIterList::iterator RbtCellTokenIterListIter;
typedef RbtCellTokenIterList::const_iterator RbtCellTokenIterListConstIter;

#endif //_RbtCellTokenIter_H_
