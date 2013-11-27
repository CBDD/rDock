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
