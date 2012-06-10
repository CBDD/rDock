/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtChromOccupancyElement.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
