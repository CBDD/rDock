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

#include "RbtGenome.h"
#include "RbtBaseSF.h"
#include "RbtChromElement.h"

RbtString RbtGenome::_CT("RbtGenome");

RbtGenome::RbtGenome(RbtChromElement* pChr)
        : m_chrom(pChr->clone()),
        m_score(0.0),
        m_RWFitness(0.0)
{
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtGenome::RbtGenome(const RbtGenome& g)
        : m_chrom( (g.m_chrom)->clone() ),
        m_score(g.m_score),
        m_RWFitness(g.m_RWFitness)
{
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtGenome& RbtGenome::operator=(const RbtGenome& g) {
  if (&g != this) {
    m_chrom = (g.m_chrom)->clone();
    m_score = g.m_score;
    m_RWFitness = g.m_RWFitness;
  }
  return *this;
}

RbtGenome::~RbtGenome() {
    delete m_chrom;
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

RbtGenome* RbtGenome::clone() const {
    return new RbtGenome(*this);
}

void RbtGenome::SetScore(RbtBaseSF* pSF) {
  if (pSF != NULL) {
    m_chrom->SyncToModel();
    m_score = -pSF->Score();
  }
  else {
    m_score = 0.0;
  }
  SetRWFitness(0.0, 0.0);
}

RbtDouble RbtGenome::SetRWFitness(RbtDouble sigmaOffset, RbtDouble partialSum) {
    //Apply sigma truncation to the raw score
    m_RWFitness = std::max(0.0, GetScore()-sigmaOffset);
    //The fitness value we store is the partial sum from previous RbtGenome elements
    //Not normalised at this point.
    m_RWFitness += partialSum;
    return m_RWFitness;
}

void RbtGenome::NormaliseRWFitness(RbtDouble total) {
    if (total > 0.0) {
        m_RWFitness /= total;
    }
}

void RbtGenome::Print(ostream& s) const {
  s << *m_chrom << endl;
  s << "Score: " << GetScore() << "; RWFitness: " << GetRWFitness() << endl;
}

ostream& operator<<(ostream& s, const RbtGenome& g) {
  g.Print(s);
  return s;
}
