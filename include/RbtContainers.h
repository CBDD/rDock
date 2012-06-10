/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtContainers.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Miscellaneous container typedefs etc for Rbt framework

#ifndef _RBTCONTAINERS_H_
#define _RBTCONTAINERS_H_

#include "RbtTypes.h"

#include <vector>
using std::vector;

#include <map>
using std::map;
using std::multimap;
using std::pair;
using std::make_pair;

#include <list>
using std::list;

#include <set>
using std::set;

#include <algorithm>
//DM 9 Dec 1998 - add copy_if template to std namespace (not included in STL)
namespace std
{
  template <class InputIterator, class OutputIterator, class Predicate>
    OutputIterator copy_if(InputIterator first, InputIterator last,
			   OutputIterator result, Predicate pred) {
    for ( ; first != last; ++first)
      if (pred(*first)) {
				*result = *first;
				++result;
      }
    return result;
  }
}

// Container Typedefs
	
// double
typedef vector<RbtDouble> RbtDoubleList;
typedef RbtDoubleList::iterator RbtDoubleListIter;
typedef RbtDoubleList::const_iterator RbtDoubleListConstIter;

// int
typedef vector<RbtInt> RbtIntList;
typedef RbtIntList::iterator RbtIntListIter;
typedef RbtIntList::const_iterator RbtIntListConstIter;

// uint
typedef vector<RbtUInt> RbtUIntList;
typedef RbtUIntList::iterator RbtUIntListIter;
typedef RbtUIntList::const_iterator RbtUIntListConstIter;

// string
typedef vector<RbtString> RbtStringList;
typedef RbtStringList::iterator RbtStringListIter;
typedef RbtStringList::const_iterator RbtStringListConstIter;

// bool
typedef vector<RbtBool> RbtBoolVec;//vector of bools
typedef RbtBoolVec::iterator RbtBoolVecIter;
typedef RbtBoolVec::const_iterator RbtBoolVecConstIter;

//Segment is a named part of an RbtModel (usually an intact molecule)
//For now, a segment is defined as just an RbtString
//RbtSegmentMap holds a map of (key=unique segment name, value=number of atoms in segment)
typedef map<RbtSegment,RbtUInt> RbtSegmentMap;
typedef RbtSegmentMap::iterator RbtSegmentMapIter;
typedef RbtSegmentMap::const_iterator RbtSegmentMapConstIter;

//Map of (key=string, value=integer)
typedef map<RbtString,RbtInt> RbtStringIntMap;
typedef RbtStringIntMap::iterator RbtStringIntMapIter;
typedef RbtStringIntMap::const_iterator RbtStringIntMapConstIter;

//Map of (key=string, value=double)
typedef map<RbtString,RbtDouble> RbtStringDoubleMap;
typedef RbtStringDoubleMap::iterator RbtStringDoubleMapIter;
typedef RbtStringDoubleMap::const_iterator RbtStringDoubleMapConstIter;

//Added by DM, 3 Dec 1998
//Map of (key=string, value=bool)
typedef map<RbtString,RbtBool> RbtStringBoolMap;
typedef RbtStringBoolMap::iterator RbtStringBoolMapIter;
typedef RbtStringBoolMap::const_iterator RbtStringBoolMapConstIter;

//Map of (key=int, value=string)
typedef map<RbtInt,RbtString> RbtIntStringMap;
typedef RbtIntStringMap::iterator RbtIntStringMapIter;
typedef RbtIntStringMap::const_iterator RbtIntStringMapConstIter;

//Set of (key=uint)
typedef set<RbtUInt> RbtUIntSet;
typedef RbtUIntSet::iterator RbtUIntSetIter;
typedef RbtUIntSet::const_iterator RbtUIntSetConstIter;  

#endif //_RBTCONTAINERS_H_
