/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtTypes.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Miscellaneous typedefs etc for Rbt framework

#ifndef _RBTTYPES_H_
#define _RBTTYPES_H_

#include <string>
using std::string;

//Typedefs (to aid portability to other platforms)

typedef double RbtDouble;
typedef int RbtInt;
typedef unsigned int RbtUInt;
typedef string RbtString;
typedef bool RbtBool;
//Segment is a named part of an RbtModel (usually an intact molecule)
//For now, a segment is defined as just an RbtString
typedef RbtString RbtSegment;

#endif //_RBTTYPES_H_
