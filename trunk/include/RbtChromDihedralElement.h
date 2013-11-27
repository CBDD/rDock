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

//Chromosome element for a single dihedral angle
#ifndef RBTCHROMDIHEDRALELEMENT_H_
#define RBTCHROMDIHEDRALELEMENT_H_

#include "RbtChromElement.h"
#include "RbtChromDihedralRefData.h"
#include "RbtRand.h"

class RbtChromDihedralElement : public RbtChromElement {
	public:
    //Class type string
    static RbtString _CT;
    //Sole constructor
    //If the tetheredAtoms list is empty, then
    //  the end of the bond with the fewest pendant atoms is rotated (other half remains fixed)
    //else if the tetheredAtoms list is not empty, then
    //  the end of the bond with the fewest tethered atoms is rotated (other half remains fixed)
    RbtChromDihedralElement(RbtBondPtr spBond,//Rotatable bond
                                RbtAtomList tetheredAtoms,//Tethered atom list
                                RbtDouble stepSize,//maximum mutation step size (degrees)
                                RbtChromElement::eMode mode=RbtChromElement::FREE,//sampling mode
                                RbtDouble maxDihedral=0.0);//max deviation from reference (tethered mode only)
    
    virtual ~RbtChromDihedralElement();
    virtual void Reset();
	virtual void Randomise();
	virtual void Mutate(RbtDouble relStepSize);
	virtual void SyncFromModel();
	virtual void SyncToModel();
	virtual RbtChromElement* clone() const;
    virtual RbtInt GetLength() const {return 1;}
    virtual RbtInt GetXOverLength() const {return 1;}
    virtual void GetVector(RbtDoubleList& v) const;
    virtual void GetVector(RbtXOverList& v) const;
    virtual void SetVector(const RbtDoubleList& v, RbtInt& i) throw (RbtError);
    virtual void SetVector(const RbtXOverList& v, RbtInt& i) throw (RbtError);
    virtual void GetStepVector(RbtDoubleList& v) const;
    virtual RbtDouble CompareVector(const RbtDoubleList& v, RbtInt& i) const;
    virtual void Print(ostream& s) const;
    
    //Returns a standardised dihedral angle in the range [-180, +180}
    //This function operates in degrees
    static RbtDouble StandardisedValue(RbtDouble dihedralAngle);
	
    protected:
    //For use by clone()
    RbtChromDihedralElement(RbtChromDihedralRefDataPtr spRefData,
                                RbtDouble value);
    RbtChromDihedralElement();
    void CorrectTetheredDihedral();
    
	private:
	RbtChromDihedralRefDataPtr m_spRefData;//Fixed reference data
	RbtDouble m_value;//The genotype value
};

#endif /*RBTCHROMDIHEDRALELEMENT_H_*/
