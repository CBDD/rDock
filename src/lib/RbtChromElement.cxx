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

#include "RbtChromElement.h"
#include "RbtRand.h"

RbtString RbtChromElement::_CT = "RbtChromElement";
RbtDouble RbtChromElement::_THRESHOLD = 1E-4;

RbtChromElement::eMode RbtChromElement::StrToMode(const RbtString& modeStr) throw (RbtError) {
    eMode retVal;
    RbtString modeStrUpper = modeStr;
    std::transform(modeStrUpper.begin(), modeStrUpper.end(), modeStrUpper.begin(), toupper);
     if (modeStrUpper == "FIXED") {
        retVal = FIXED;
    }
    else if (modeStrUpper == "TETHERED") {
        retVal = TETHERED;
    }
    else if (modeStrUpper == "FREE") {
        retVal = FREE;
    }
    else {
        throw RbtBadArgument(_WHERE_,
                              "Unknown chromosome sampling mode ("
                              + modeStr + ")");
    }
    return retVal;
}

RbtString RbtChromElement::ModeToStr(RbtChromElement::eMode mode) {
    RbtString retVal;
    switch (mode) {
    case FIXED:
        retVal = "FIXED";
        break;
    case TETHERED:
        retVal = "TETHERED";
        break;
    case FREE:
        retVal = "FREE";
        break;
    default:
        retVal = "UNDEFINED";
        break;
    }
    return retVal;
}

RbtChromElement::RbtChromElement()
        : m_rand(Rbt::GetRbtRand())
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChromElement::~RbtChromElement() {
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtChromElement::Add(RbtChromElement* pChromElement) throw (RbtError) {
    throw RbtInvalidRequest(_WHERE_,
        "Add(RbtChromElement*) invalid for non-aggregate chromosome element");
}

RbtBool RbtChromElement::VectorOK(const RbtDoubleList& v, RbtInt i) const {
    RbtInt length = GetLength();
    //if the element is empty then any vector is valid
    return (length == 0)
            || ( (i >= 0) && (i < v.size()) && (length <= (v.size() - i) ));
}

RbtBool RbtChromElement::VectorOK(const RbtXOverList& v, RbtInt i) const {
    RbtInt length = GetXOverLength();
    //if the element is empty then any vector is valid
    return (length == 0)
            || ( (i >= 0) && (i < v.size()) && (length <= (v.size() - i) ));
}

void RbtChromElement::CauchyMutate(RbtDouble mean, RbtDouble variance) {
    //Need to convert the Cauchy random variable to a positive number
    //and use this as the relative step size for mutation
    RbtDouble relStepSize = fabs(m_rand.GetCauchyRandom(mean, variance));
    Mutate(relStepSize);
}

RbtDouble RbtChromElement::Compare(const RbtChromElement& c) const {
    RbtDouble retVal(0.0);
    if (GetLength() != c.GetLength()) {
        retVal = -1.0;
    }
    else {
        RbtDoubleList v;
        RbtInt i(0);
        c.GetVector(v);
        retVal = CompareVector(v,i);
    }
    return retVal;
}

RbtBool RbtChromElement::Equals(const RbtChromElement& c, RbtDouble threshold) const {
    RbtDouble cmp = Compare(c);
    return ( (cmp >= 0.0) && (cmp < threshold) );
}

void RbtChromElement::SetVector(const RbtDoubleList& v) {
    RbtInt i(0);
    SetVector(v,i);
}

void RbtChromElement::SetVector(const RbtXOverList& v) {
    RbtInt i(0);
    SetVector(v,i);
}

bool operator== (const RbtChromElement& c1, const RbtChromElement& c2) {
    return c1.Equals(c2, RbtChromElement::_THRESHOLD);
}

bool operator!= (const RbtChromElement& c1, const RbtChromElement& c2) {
    return !(c1 == c2);
}

ostream& operator<<(ostream& s, const RbtChromElement& c) {
    c.Print(s);
    return s;   
}


void Rbt::Crossover(RbtChromElement* pChr1, RbtChromElement* pChr2,
                        RbtChromElement* pChr3, RbtChromElement* pChr4)
                        throw (RbtError) {
  //Check all chromosomes have the same crossover length
  RbtInt length1 = pChr1->GetXOverLength();
  if ( (length1 != pChr2->GetXOverLength())
        || (length1 != pChr3->GetXOverLength())
        || (length1 != pChr4->GetXOverLength()) ) {
    throw RbtBadArgument(_WHERE_,"Crossover: mismatch in chromosome lengths");
  }
  //Extract the gene vectors from each parent
  RbtXOverList v1, v2;
  pChr1->GetVector(v1);
  pChr2->GetVector(v2);
  //2-point crossover
  //In the spirit of STL, ixbegin is the first gene to crossover, ixend is one after the last gene to crossover
  RbtRand& rand = pChr1->GetRand();
  RbtInt ixbegin = rand.GetRandomInt(length1);
  //if ixbegin is 0, we need to avoid selecting the whole chromosome
  RbtInt ixend = (ixbegin == 0) ? rand.GetRandomInt(length1-1)+1 
                                : rand.GetRandomInt(length1-ixbegin)+ixbegin+1;

  //cout << "XOVER: ixbegin = " << ixbegin << ", ixend = " << ixend << endl;
  std::swap_ranges(v1.begin()+ixbegin,v1.begin()+ixend,v2.begin()+ixbegin);
  //Now we can update the two children
  pChr3->SetVector(v1);
  pChr4->SetVector(v2);                      
}

