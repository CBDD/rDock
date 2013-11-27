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

//Chromosome element for model centre of mass and orientation
#ifndef RBTCHROMPOSITIONLEMENT_H_
#define RBTCHROMPOSITIONELEMENT_H_

#include "RbtChromElement.h"
#include "RbtChromPositionRefData.h"
#include "RbtRand.h"

class RbtChromPositionElement : public RbtChromElement {
	public:
    //Class type string
    static RbtString _CT;

	RbtChromPositionElement(const RbtModel* pModel,
                                const RbtDockingSite* pDockSite,
                                RbtDouble transStepSize,//Angstroms
                                RbtDouble rotStepSize,//Radians
                                RbtChromElement::eMode transMode=
                                            RbtChromElement::FREE,
                                RbtChromElement::eMode rotMode=
                                            RbtChromElement::FREE,
                                RbtDouble maxTrans=0.0,//Angstroms
                                RbtDouble maxRot=0.0);//radians
    virtual ~RbtChromPositionElement();
    virtual void Reset();
	virtual void Randomise();
	virtual void Mutate(RbtDouble relStepSize);
	virtual void SyncFromModel();
	virtual void SyncToModel();
	virtual RbtChromElement* clone() const;
    virtual RbtInt GetLength() const {return m_spRefData->GetLength();}
    virtual RbtInt GetXOverLength() const {return m_spRefData->GetXOverLength();}
    virtual void GetVector(RbtDoubleList& v) const;
    virtual void GetVector(RbtXOverList& v) const;
    virtual void SetVector(const RbtDoubleList& v, RbtInt& i) throw (RbtError);
    virtual void SetVector(const RbtXOverList& v, RbtInt& i) throw (RbtError);
    virtual void GetStepVector(RbtDoubleList& v) const;
    virtual RbtDouble CompareVector(const RbtDoubleList& v, RbtInt& i) const;
    virtual void Print(ostream& s) const;

    //Returns a standardised rotation angle in the range [-M_PI, +M_PI}
    //This function operates in radians
    static RbtDouble StandardisedValue(RbtDouble rotationAngle);
	
    protected:
    //For use by clone()
    RbtChromPositionElement(RbtChromPositionRefDataPtr spRefData,
                                const RbtCoord& com,
                                const RbtEuler& orientation);
    RbtChromPositionElement();
    void RandomiseCOM();
    void RandomiseOrientation();
    void MutateCOM(RbtDouble relStepSize);
    void MutateOrientation(RbtDouble relStepSize);
    void CorrectTetheredCOM();
    void CorrectTetheredOrientation();
    
	private:
	RbtChromPositionRefDataPtr m_spRefData;//Fixed reference data
	RbtCoord m_com;//Centre of mass genotype value
    RbtEuler m_orientation;//Euler angle orientation genotype value
};

#endif /*RBTCHROMPOSITIONELEMENT_H_*/
