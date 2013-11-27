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

#include "RbtChrom.h"

RbtString RbtChrom::_CT = "RbtChrom";

RbtChrom::RbtChrom() : RbtChromElement() {
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChrom::RbtChrom(const RbtModelList& modelList)
    : RbtChromElement(), m_modelList(modelList) {
    for (RbtModelListConstIter iter = m_modelList.begin();
                               iter != m_modelList.end();
                               ++iter) {
        if ((*iter).Ptr()) {
            Add((*iter)->GetChrom());
        }
    }     
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtChrom::~RbtChrom() {
    for (RbtChromElementListIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        delete *iter;
    }
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

//pure virtual from RbtChromElement
void RbtChrom::Reset() {
    for (RbtChromElementListIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        (*iter)->Reset();
    }
}

void RbtChrom::Randomise() {
    for (RbtChromElementListIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        (*iter)->Randomise();
    }
}

void RbtChrom::Mutate(RbtDouble relStepSize) {
    for (RbtChromElementListIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        (*iter)->Mutate(relStepSize);
    }   
}

void RbtChrom::SyncFromModel() {
    for (RbtChromElementListIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        (*iter)->SyncFromModel();
    }    
}

void RbtChrom::SyncToModel() {
    for (RbtChromElementListIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        (*iter)->SyncToModel();
    }
    //Force an update of all the pseudo atom coords in each model
    for (RbtModelListIter iter = m_modelList.begin();
                          iter != m_modelList.end();
                          ++iter) {
        (*iter)->UpdatePseudoAtoms();
    }
}

RbtChromElement* RbtChrom::clone() const {
    RbtChromElement* clone = new RbtChrom();
    for (RbtChromElementListConstIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        clone->Add((*iter)->clone());
    }
    return clone;        
}

void RbtChrom::Add(RbtChromElement* pChromElement) throw (RbtError) {
    if (pChromElement) {
        m_elementList.push_back(pChromElement);
    }
}

RbtInt RbtChrom::GetLength() const {
    RbtInt retVal(0);
    for (RbtChromElementListConstIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        retVal += (*iter)->GetLength();
    }
    return retVal;
}

RbtInt RbtChrom::GetXOverLength() const {
    RbtInt retVal(0);
    for (RbtChromElementListConstIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        retVal += (*iter)->GetXOverLength();
    }
    return retVal;
}

void RbtChrom::GetVector(RbtDoubleList& v) const {
    for (RbtChromElementListConstIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        (*iter)->GetVector(v);
    }    
}

void RbtChrom::GetVector(RbtXOverList& v) const {
    for (RbtChromElementListConstIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        (*iter)->GetVector(v);
    }    
}

void RbtChrom::SetVector(const RbtDoubleList& v, RbtInt& i) throw (RbtError) {
    if (VectorOK(v,i)) {
        for (RbtChromElementListIter iter = m_elementList.begin();
                iter != m_elementList.end(); ++iter) {
            (*iter)->SetVector(v, i);
        }
    }
    else {
        throw RbtBadArgument(_WHERE_, "Index i out of range or insufficient elements remaining");
    }
}

void RbtChrom::SetVector(const RbtXOverList& v, RbtInt& i) throw (RbtError) {
    if (VectorOK(v,i)) {
        for (RbtChromElementListIter iter = m_elementList.begin();
                iter != m_elementList.end(); ++iter) {
            (*iter)->SetVector(v, i);
        }
    }
    else {
        throw RbtBadArgument(_WHERE_, "Index i out of range or insufficient elements remaining");
    }
}

void RbtChrom::GetStepVector(RbtDoubleList& v) const {
    for (RbtChromElementListConstIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        (*iter)->GetStepVector(v);
    }       
}

//Returns the maximum difference of any of the underlying chromosome elements
//or -1 if there is a mismatch in vector sizes
RbtDouble RbtChrom::CompareVector(const RbtDoubleList& v, RbtInt& i) const {
    RbtDouble retVal(0.0);
    if (VectorOK(v,i)) {
        RbtChromElementListConstIter iter = m_elementList.begin();
        //Exit the while loop as soon as a negative value is returned
        while ( (iter != m_elementList.end()) && (retVal >= 0.0) ) {
            RbtChromElement* element = *(iter++);
            RbtDouble cmp = element->CompareVector(v,i);
            retVal = (cmp < 0) ? cmp : std::max(retVal, cmp);
        }
    }
    else {
        retVal = -1.0;
    }
    return retVal;
}

void RbtChrom::Print(ostream& s) const {
    s << "CHROM" << endl;
    RbtInt i(0);
    for (RbtChromElementListConstIter iter = m_elementList.begin();
            iter != m_elementList.end(); ++iter) {
        s << "ELEMENT " << (i++) << ": ";
        (*iter)->Print(s);
    }
}
