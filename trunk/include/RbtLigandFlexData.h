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

//Manages flexibility data for ligands
#ifndef _RBTLIGANDFLEXDATA_
#define _RBTLIGANDFLEXDATA_

#include "RbtFlexData.h"

class RbtLigandFlexData : public RbtFlexData {
public:
    //Translation mutation step length (A)
    static const RbtString& _TRANS_STEP;
    //Whole-body orientation mutation step length (deg)
    static const RbtString& _ROT_STEP;
    //Dihedral angle mutation step length (deg)
    static const RbtString& _DIHEDRAL_STEP;
    //Enum for translation mode (0=FIXED, 1=TETHERED, 2=FREE)
    static const RbtString& _TRANS_MODE;
    //Enum for rotation mode (0=FIXED, 1=TETHERED, 2=FREE)
    static const RbtString& _ROT_MODE;
    //Enum for dihedral mode (0=FIXED, 1=TETHERED, 2=FREE)
    static const RbtString& _DIHEDRAL_MODE;
    //Max allowed distance from initial COM (for tethered translation only)
    static const RbtString& _MAX_TRANS;
    //Max allowed rotation from initial orientation
    //(for tethered rotation only)
    static const RbtString& _MAX_ROT;
    //Max allowed dihedral rotation from initial dihedrals
    //(for tethered dihedrals only)
    static const RbtString& _MAX_DIHEDRAL;
    RbtLigandFlexData(RbtDockingSite* pDockSite);
	virtual void Accept(RbtFlexDataVisitor& v) {v.VisitLigandFlexData(this);}
         
private:
    RbtLigandFlexData();
};
#endif //_RBTLIGANDFLEXDATA_
