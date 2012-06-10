/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtCellTokenIter.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
