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

//Common includes for Rbt framework

#ifndef _RBTCONFIG_H_
#define _RBTCONFIG_H_

#include <iostream>

//For g++ 3.4.3 (libstdc++)
#include <strstream>
using std::cout;
using std::cerr;
using std::endl;
using std::ends;
using std::istream;
using std::ostream;
using std::ostrstream;
using std::istrstream;
using std::ios_base;

#ifndef _RBTTYPES_H_
#include "RbtTypes.h"
#endif //_RBTTYPES_H_

#ifndef _RBTSMARTPOINTER_H_
#include "RbtSmartPointer.h"
#endif //_RBTSMARTPOINTER_H_

#ifndef _RBTCONTAINERS_H_
#include "RbtContainers.h"
#endif //_RBTCONTAINERS_H_

#ifndef _RBTDEBUG_H_
#include "RbtDebug.h"
#endif //_RBTDEBUG_H_

#ifndef _RBT_H_
#include "Rbt.h"
#endif //_RBT_H_

#ifndef _RBTERROR_H_
#include "RbtError.h"
#endif //_RBTERROR_H_

#endif //_RBTCONFIG_H_
