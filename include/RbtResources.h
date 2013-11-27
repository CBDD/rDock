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

//Miscellaneous string constants

#ifndef _RBTRESOURCES_H_
#define _RBTRESOURCES_H_

#include "RbtTypes.h"
#include "VERSION"

const RbtString IDS_COPYRIGHT = "The rDock program is licensed under GNU-LGPLv3.0. http://rdock.sourceforge.net/";
#ifdef _DEBUG
const RbtString IDS_PRODUCT   = "libRbtD.so";
#else
const RbtString IDS_PRODUCT   = "libRbt.so";
#endif //_DEBUG

#endif //_RBTRESOURCES_H_
