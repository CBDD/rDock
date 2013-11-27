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
