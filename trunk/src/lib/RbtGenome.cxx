/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

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
