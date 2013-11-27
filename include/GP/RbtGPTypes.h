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
