/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtTokenIter.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
