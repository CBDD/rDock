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

#ifndef _RBTPMF_H_
#define _RBTPMF_H_

// #include <vector>		// STL vector
#include "RbtTypes.h"
/**
 * Basic data type to represent distance-value pairs of the PMF function
 */
typedef struct {
	float distance;
	float density;
} RbtPMFValue;

typedef enum {
		CF=1, CP, cF, cP, C3, CW, CO, CN, C0, NC,	// 10
		  NP, NA, ND, NR, N0, NS, OC, OA, OE, OR, 	// 20
		  OS, OD, OW, P,  SA, SD, HL, HH, Zn, CL, 	// 30
		  Mn, Mg,  F, Fe, Br, V,  PMF_UNDEFINED
} RbtPMFType;

typedef RbtPMFType PMFType;

namespace Rbt {
	RbtString	PMFType2Str(RbtPMFType aType);
	RbtPMFType	PMFStr2Type(RbtString anStr);
};

#endif // _RBTPMF_H_
