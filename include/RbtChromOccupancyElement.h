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

//Chromosome element for solvent model occupancy
#ifndef RBTCHROMOCCUPANCYELEMENT_H_
#define RBTCHROMOCCUPANCYELEMENT_H_

#include "RbtChromElement.h"
#include "RbtChromOccupancyRefData.h"
#include "RbtRand.h"

class RbtChromOccupancyElement : public RbtChromElement {
  public:
    //Class type string
    static RbtString _CT;
    //Sole constructor
    RbtChromOccupancyElement(RbtModel* pModel,
                             RbtDouble stepSize,
                             RbtDouble threshold);
    
    virtual ~RbtChromOccupancyElement();
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
   
    //Returns a standardised occupancy value in the range [0,1]
    static RbtDouble StandardisedValue(RbtDouble occupancy);
 
  protected:
    //For use by clone()
    RbtChromOccupancyElement(RbtChromOccupancyRefDataPtr spRefData,
                                RbtDouble value);
    RbtChromOccupancyElement();
  private:
    RbtChromOccupancyRefDataPtr m_spRefData;//Fixed reference data
    RbtDouble m_value;//The genotype value
};

#endif /*RBTCHROMOCCUPANCYELEMENT_H_*/
