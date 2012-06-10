/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/GP/RbtGPTypes.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#ifndef _RBT_GPTYPES_H_
#define _RBT_GPTYPES_H_

#include "RbtConfig.h"

// As the chromosome is a list of integers, I don't think I need a
// new class to define it. If later I find out it is necessary, I
// just need to create a class called RbtGPChromosome with the
// operator[] and size() defined on it
typedef float RbtReturnType;
typedef SmartPtr<RbtReturnType> RbtReturnTypePtr;  //Smart pointer
typedef vector<RbtReturnTypePtr> RbtReturnTypeList;//Vector of smart pointers
typedef vector<RbtReturnTypeList> RbtReturnTypeArray;//Vector of smart pointers
typedef RbtReturnTypeList::iterator RbtReturnTypeListIter;
typedef RbtReturnTypeList::const_iterator RbtReturnTypeListConstIter;

#endif
