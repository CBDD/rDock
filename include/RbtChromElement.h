/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtChromElement.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Abstract interface for all chromosome elements (degrees of freedom)
#ifndef RBTCHROMELEMENT_H_
#define RBTCHROMELEMENT_H_

#include "RbtConfig.h"
#include "RbtRand.h"

//Typedefs for crossover data.
//To prevent splitting (for example) an orientation or position vector in two
//during crossover, we convert the chromosome to a vector of vector of doubles
//An orientation vector (3 Euler angles) is a single RbtXOverElement
//GetVector(RbtDoubleList&) returns a flat vector of doubles (for use by Simplex)
//GetVector(RbtXOverList&) returns a vector of RbtXOverElements (for use by crossover)
typedef vector<RbtDouble> RbtXOverElement;
typedef RbtXOverElement::iterator RbtXOverElementIter;
typedef RbtXOverElement::const_iterator RbtXOverElementConstIter;
typedef vector<RbtXOverElement> RbtXOverList;
typedef RbtXOverList::iterator RbtXOverListIter;
typedef RbtXOverList::const_iterator RbtXOverListConstIter;

class RbtChromElement {
	public:
    //Class type string
    static RbtString _CT;
    //Threshold used to assess equality of two chromosome elements
    static RbtDouble _THRESHOLD;
    
    //General enum for defining types of chromosome element flexibility.
    //For use by subclasses of RbtChromElement if desired.
    //There is no functionality in the base class that uses this enum.
    enum eMode {
            FIXED = 0,
            TETHERED = 1,
            FREE = 2
    };
    //Static methods to convert from mode enum to string and vice versa
    static eMode StrToMode(const RbtString& modeStr) throw (RbtError);//case insensitive
    static RbtString ModeToStr(eMode mode);//returns "FIXED", "TETHERED" or "FREE"ss

 	virtual ~RbtChromElement();
    //
    //INTERFACE METHODS - must be implemented in subclasses
    //
    //Resets the chromosome element to its initial value (extracted from
    //model coords at time of construction)
    virtual void Reset() = 0;
    //Randomises the chromosome element within the context of any restraints
	virtual void Randomise() = 0;
    //Mutates the chromosome element
    //relStepSize = fraction (0.0->1.0) of maximum defined step size
	virtual void Mutate(RbtDouble relStepSize) = 0;
    //Updates the chromosome element to match the current model coords
	virtual void SyncFromModel() = 0;
    //Updates the model coords to match the chromosome element
	virtual void SyncToModel() = 0;
    //Creates an independent clone
    //The current value of the cloned element should match the current value
    //of this element. 
	virtual RbtChromElement* clone() const = 0;
    //Gets the number of double values needed to represent this
    //chromosome element
    virtual RbtInt GetLength() const = 0;
    //Gets the number of vectors of double values needed to represent this
    //chromosome element for crossover purposes
    //e.g. an orientation vector consists of 3 doubles, but should be crossed over
    //as a single intact entity 
    virtual RbtInt GetXOverLength() const = 0;
    //Converts chromosome element to a vector of double values,
    //Number of double values appended should match GetLength().
    //v = vector of doubles to append to.
    virtual void GetVector(RbtDoubleList& v) const = 0;
    //Converts chromosome element to a vector of XOverElements
    //where each XOverElement is itself a vector of double values
    //v = vector of RbtXOverElements to append to.
    virtual void GetVector(RbtXOverList& v) const = 0;
    //Updates chromosome element from a vector of double values,
    //Number of double values read should match GetLength().
    //v = vector of doubles to extract from
    //i = index of next vector element to read (should be updated by method)
    virtual void SetVector(const RbtDoubleList& v, RbtInt& i) throw (RbtError) = 0;
    //Updates chromosome element from a vector of XOverElements,
    //Number of XOverElements read should match GetXOverLength().
    //v = vector of XOverElements to extract from
    //i = index of next vector element to read (should be updated by method)
    virtual void SetVector(const RbtXOverList& v, RbtInt& i) throw (RbtError) = 0;
    //Gets the vector of absolute step sizes that correspond to each double value.
    virtual void GetStepVector(RbtDoubleList& v) const = 0;
    //Gets the maximum relative difference between this element and another element
    //as represented by the vector. Differences should be normalised by the
    //step size.
    //This method is provided to allow better comparisons between values
    //that are difficult to compare by simple numerical differences. e.g.
    //Dihedral angles are cyclical, therefore -180 and + 179 only differ by 1 deg.
    virtual RbtDouble CompareVector(const RbtDoubleList& v, RbtInt& i) const = 0;
    //
    //IMPLEMENTED VIRTUAL METHODS
    //
    //Invalid operation in base class
    virtual void Add(RbtChromElement* pChromElement) throw (RbtError);
    //Prints details of element to stream (null implementation in base class)
    virtual void Print(ostream& s) const {};
    //
    //NON-VIRTUAL METHODS
    //
    RbtRand& GetRand() const {return m_rand;}
    void CauchyMutate(RbtDouble mean, RbtDouble variance);
    //Compares two chromosome elements. Returns -1 if the comparison is invalid
    //(unequal lengths), else returns the maximum relative pair-wise difference
    //as returned by CompareVector.
    RbtDouble Compare(const RbtChromElement& c) const;
    //Returns true if chromosome elements have near-equal values.
    //i.e. if Compare(c) < threshold
    //Returns false if comparison is invalid (unequal lengths)
    RbtBool Equals(const RbtChromElement& c, RbtDouble threshold) const;
    //Convenience method that calls SetVector(const RbtDoubleList& v, RbtInt& i)
    //with i initialised to zero
    void SetVector(const RbtDoubleList& v);
    //Convenience method that calls SetVector(const RbtXOverList& v, RbtInt& i)
    //with i initialised to zero
    void SetVector(const RbtXOverList& v);
    //operator== and operator!= are implemented by calling Equals with the
    //static _THRESHOLD value
    friend bool operator== (const RbtChromElement& c1, const RbtChromElement& c2);
    friend bool operator!= (const RbtChromElement& c1, const RbtChromElement& c2);
    //operator<< is implemented by calling Print()
    friend ostream& operator<<(ostream& s, const RbtChromElement& c);
    
    protected:
    RbtChromElement();
    RbtChromElement(const RbtChromElement& c);
    RbtChromElement& operator=(const RbtChromElement& c);
    //Helper method to check that index i is in range [0,v.size()}
    //and that v has sufficient elements remaining to satisfy
    //GetLength()
    RbtBool VectorOK(const RbtDoubleList& v, RbtInt i) const;
    //Helper method to check that index i is in range [0,v.size()}
    //and that v has sufficient elements remaining to satisfy
    //GetXOverLength()
    RbtBool VectorOK(const RbtXOverList& v, RbtInt i) const;
    
    private:
    RbtRand& m_rand;//Reference to singleton random number generator
};

typedef SmartPtr<RbtChromElement> RbtChromElementPtr;
typedef vector<RbtChromElement*> RbtChromElementList;//Vector of regular pointers
typedef RbtChromElementList::iterator RbtChromElementListIter;
typedef RbtChromElementList::const_iterator RbtChromElementListConstIter;

namespace Rbt {
    //2-point crossover
    void Crossover(RbtChromElement* pChr1, RbtChromElement* pChr2,
                RbtChromElement* pChr3, RbtChromElement* pChr4) throw (RbtError);
}
#endif /*RBTCHROMELEMENT_H_*/
