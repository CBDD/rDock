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
