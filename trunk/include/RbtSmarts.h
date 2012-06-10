/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtSmarts.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
