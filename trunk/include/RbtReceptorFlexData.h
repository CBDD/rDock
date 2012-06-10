/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtReceptorFlexData.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Manages flexibility data for macromolecules
#ifndef _RBTRECEPTORFLEXDATA_
#define _RBTRECEPTORFLEXDATA_

#include "RbtFlexData.h"

class RbtReceptorFlexData : public RbtFlexData {
public:
    //Distance from docking site within which OH/NH3+ groups are flexible (A)
    static const RbtString& _FLEX_DISTANCE;
    //Dihedral mutation step length (deg)
    static const RbtString& _DIHEDRAL_STEP;
    RbtReceptorFlexData(RbtDockingSite* pDockSite);
	virtual void Accept(RbtFlexDataVisitor& v) {v.VisitReceptorFlexData(this);}
     
private:
    RbtReceptorFlexData();
};
#endif //_RBTRECEPTORFLEXDATA_
