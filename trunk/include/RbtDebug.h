/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtDebug.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Tools for tracking of object construction and destruction for
//debugging purposes

#ifndef _RBTDEBUG_H_
#define _RBTDEBUG_H_

#ifdef _RBTOBJECTCOUNTER

#include <iostream.h>
#include <map>
#include "RbtTypes.h"

//Object counter class for debugging no. of object constructions and destructions
class RbtObjectCounter
{
 public:
  RbtObjectCounter();
  ~RbtObjectCounter();
  friend ostream& operator<<(ostream& s, const RbtObjectCounter& counter);
  RbtInt nConstr;
  RbtInt nCopyConstr;
  RbtInt nDestr;
};

//Useful typedefs
typedef std::map<RbtString,RbtObjectCounter> RbtObjectCounterMap;
typedef RbtObjectCounterMap::iterator RbtObjectCounterMapIter;

namespace Rbt
{
  //External declaration of the one and only object counter map
  extern RbtObjectCounterMap theObjectCounterMap;
}

#define _RBTOBJECTCOUNTER_CONSTR_(class) Rbt::theObjectCounterMap[(class)].nConstr++
#define _RBTOBJECTCOUNTER_COPYCONSTR_(class) Rbt::theObjectCounterMap[(class)].nCopyConstr++
#define _RBTOBJECTCOUNTER_DESTR_(class) Rbt::theObjectCounterMap[(class)].nDestr++

#define _RBTOBJECTCOUNTER_DUMP_(stream) \
(stream) << endl << "RBT OBJECT COUNTERS:" << endl << endl;\
for (RbtObjectCounterMapIter dumpIter = Rbt::theObjectCounterMap.begin();\
dumpIter != Rbt::theObjectCounterMap.end(); dumpIter++) {\
  (stream) << (*dumpIter).first << ": " << (*dumpIter).second << endl;}

#else //_RBTOBJECTCOUNTER

#define _RBTOBJECTCOUNTER_CONSTR_(class)
#define _RBTOBJECTCOUNTER_COPYCONSTR_(class)
#define _RBTOBJECTCOUNTER_DESTR_(class)
#define _RBTOBJECTCOUNTER_DUMP_(stream)

#endif //_RBTOBJECTCOUNTER

#endif //_RBTDEBUG_H_
