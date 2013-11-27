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

//Abstract TokenIter class. Defines the interface for a token iterator
//to be used by the parser.

#ifndef _RBT_TOKENITER_H_
#define _RBT_TOKENITER_H_

#include "RbtTypes.h"
#include "RbtToken.h"
#include "RbtContext.h"


class RbtTokenIter 
{
public:
//  virtual void copy(const RbtTokenIter &) = 0;
  virtual void Next(RbtContextPtr) = 0;
  virtual RbtTokenPtr Current() = 0;
  virtual ~RbtTokenIter(){}
  
  //  virtual RbtTokenPtr IsDone() const = 0;
  //  virtual void First() = 0;
  //  RbtBool operator==(const RbtTokenIter&) const;
  //  ostream& Print(ostream&) const;
  //  friend ostream& operator<<(ostream& s, const RbtTokenIter &p);

};

//Useful typedefs
typedef SmartPtr<RbtTokenIter> RbtTokenIterPtr;  //Smart pointer
typedef vector<RbtTokenIterPtr> RbtTokenIterList;//Vector of smart pointers
typedef RbtTokenIterList::iterator RbtTokenIterListIter;
typedef RbtTokenIterList::const_iterator RbtTokenIterListConstIter;

#endif //_RbtTokenIter_H_
