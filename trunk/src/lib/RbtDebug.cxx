/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtDebug.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#ifdef _RBTOBJECTCOUNTER
#include "RbtDebug.h"

RbtObjectCounter::RbtObjectCounter() :
  nConstr(0),nCopyConstr(0),nDestr(0)
{
}

RbtObjectCounter::~RbtObjectCounter()
{
}

ostream& operator<<(ostream& s, const RbtObjectCounter& counter)
{
  return s << "Constr=" << counter.nConstr
	   << ", Copy constr=" << counter.nCopyConstr
	   << ", Destr=" << counter.nDestr
	   << ", Net total=" << counter.nConstr + counter.nCopyConstr - counter.nDestr;
}

//The one and only RbtObjectCounterMap definition
RbtObjectCounterMap Rbt::theObjectCounterMap;

#endif //_RBTOBJECTCOUNTER


