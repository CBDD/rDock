/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtPMF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
