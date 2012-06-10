/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtConfig.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
