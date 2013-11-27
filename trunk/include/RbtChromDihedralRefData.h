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

//Manages the fixed reference data for a single dihedral angle chromosome element
//Also provides methods to map the genotype (dihedral angle value) onto the
//phenotype (model coords)
//A single instance is designed to be shared between all clones of a given element
#ifndef RBTCHROMDIHEDRALREFDATA_H_
#define RBTCHROMDIHEDRALREFDATA_H_

#include "RbtAtom.h"
#include "RbtBond.h"
#include "RbtChromElement.h"

class RbtChromDihedralRefData {
	public:
        //Class type string
        static RbtString _CT;
       //Sole constructor
        //If the tetheredAtoms list is empty, then
        //  the end of the bond with the fewest pendant atoms is rotated (other half remains fixed)
        //else if the tetheredAtoms list is not empty, then
        //  the end of the bond with the fewest tethered atoms is rotated (other half remains fixed)
        RbtChromDihedralRefData(RbtBondPtr spBond,//Rotatable bond
                                    RbtAtomList tetheredAtoms,//Tethered atom list
                                    RbtDouble stepSize,//maximum mutation step size (degrees)
                                    RbtChromElement::eMode mode=RbtChromElement::FREE,//sampling mode
                                    RbtDouble maxDihedral=0.0);//max deviation from reference (tethered mode only)
         virtual ~RbtChromDihedralRefData();
        
        //Gets the maximum step size for this bond
		RbtDouble GetStepSize() const {return m_stepSize;}
        //Gets the sampling mode for this bond
        RbtChromElement::eMode GetMode() const {return m_mode;}
        //Gets the maximum deviation from the reference dihedral (tethered mode only)
        RbtDouble GetMaxDihedral() const {return m_maxDihedral;}
        //Gets the current dihedral angle (phenotype) for this bond
        //from the model coords
        RbtDouble GetModelValue() const;
        //Sets the phenotype (model coords) for this bond
        //to a given dihedral angle
		void SetModelValue(RbtDouble dihedralAngle);
        //Gets the initial dihedral angle for this bond
        //(initialised from model coords in RbtChromDihedralRefData constructor)
        RbtDouble GetInitialValue() const {return m_initialValue;}
								
	private:
        //Sets up the movable atom list for this bond
		void Setup(RbtBondPtr spBond, const RbtAtomList& tetheredAtoms);
		RbtAtom* m_atom1;
		RbtAtom* m_atom2;
		RbtAtom* m_atom3;
		RbtAtom* m_atom4;
		RbtAtomRList m_rotAtoms;
		RbtDouble m_stepSize;
        RbtDouble m_initialValue;
        RbtChromElement::eMode m_mode;
        RbtDouble m_maxDihedral;//max deviation from reference (tethered mode only)
};

typedef SmartPtr<RbtChromDihedralRefData> RbtChromDihedralRefDataPtr;//Smart pointer

#endif /*RBTCHROMDIHEDRALREFDATA_H_*/
