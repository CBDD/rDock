/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtStringTokenIter.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Abstract StringTokenIter class. Defines the interface for a token iterator
//to be used by the parser.

#ifndef _RBT_STRINGTOKENITER_H_
#define _RBT_STRINGTOKENITER_H_

#include "RbtConfig.h"
#include "RbtTokenIter.h"
#include "RbtGPChromosome.h"
#include "RbtContext.h"
#include <stack>
#include <fstream>


typedef SmartPtr<istream> istreamPtr;
class RbtStringTokenIter : public RbtTokenIter
{
public:

  static RbtString _CT;
    ///////////////////
    // Constructors
    ///////////////////
  RbtStringTokenIter(const RbtStringTokenIter&);
  RbtStringTokenIter(const istreamPtr, RbtContextPtr);
 
    ///////////////////
    // Destructor
    //////////////////
  virtual ~RbtStringTokenIter();
  
  void copy(const RbtStringTokenIter &);
  //void First();
  void Next(RbtContextPtr);
  //RbtBool IsDone() const;
  RbtTokenPtr Current(); 
 
//  RbtBool operator==(const RbtStringTokenIter&) const;
//  ostream& Print(ostream&) const;
//  friend ostream& operator<<(ostream& s, const RbtStringTokenIter &p);

    ////////////////////
    // Private methods
    ////////////////////

private:
  RbtStringTokenIter();  // Default constructor disabled
  RbtTokenPtr translate(RbtString);

  istreamPtr filep;
  RbtTokenPtr current;
  RbtString strtok;
  RbtContextPtr contextp;
  
};

//Useful typedefs
typedef SmartPtr<RbtStringTokenIter> RbtStringTokenIterPtr;  //Smart pointer
typedef vector<RbtStringTokenIterPtr> RbtStringTokenIterList;//Vector of smart pointers
typedef RbtStringTokenIterList::iterator RbtStringTokenIterListIter;
typedef RbtStringTokenIterList::const_iterator RbtStringTokenIterListConstIter;

#endif //_RbtStringTokenIter_H_
