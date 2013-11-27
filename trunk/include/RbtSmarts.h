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

//SMARTS class. Defines a function that takes a model and a smarts string as
//input, and returns a list of lists of atoms that match the pattern

#ifndef _RBT_SMARTS_H_
#define _RBT_SMARTS_H_

#include "RbtModel.h"

namespace DT {
  RbtAtomListList QueryModel(RbtModelPtr spModel, const RbtString& strSmart, RbtString& strSmiles);
}

//Map of (key=integer, value=integer)
typedef map<RbtInt,RbtInt> RbtIntIntMap;
typedef RbtIntIntMap::iterator RbtIntIntMapIter;
typedef RbtIntIntMap::const_iterator RbtIntIntMapConstIter;

#endif
