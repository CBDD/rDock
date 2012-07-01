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

#include "RbtRandPopTransform.h"
#include "RbtPopulation.h"
#include "RbtWorkSpace.h"
#include "RbtChrom.h"

RbtString RbtRandPopTransform::_CT("RbtRandPopTransform");
RbtString RbtRandPopTransform::_POP_SIZE("POP_SIZE");
RbtString RbtRandPopTransform::_SCALE_CHROM_LENGTH("SCALE_CHROM_LENGTH");

RbtRandPopTransform::RbtRandPopTransform(const RbtString& strName) :
  RbtBaseBiMolTransform(_CT,strName) {
  AddParameter(_POP_SIZE, 50);
  AddParameter(_SCALE_CHROM_LENGTH, true);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtRandPopTransform::~RbtRandPopTransform() {
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}


////////////////////////////////////////
//Protected methods
///////////////////
void RbtRandPopTransform::SetupReceptor() {}

void RbtRandPopTransform::SetupLigand() {}

void RbtRandPopTransform::SetupSolvent() {}

void RbtRandPopTransform::SetupTransform() {
    //Construct the overall chromosome for the system
    m_chrom = new RbtChrom(GetWorkSpace()->GetModels());
}

////////////////////////////////////////
//Private methods
///////////////////
//Pure virtual in RbtBaseTransform
//Actually apply the transform
void RbtRandPopTransform::Execute() {
    if (m_chrom.Ptr() == NULL) {
        return;
    }
    RbtBaseSF* pSF = GetWorkSpace()->GetSF();
    if (pSF == NULL) {
        return;
    }
    RbtInt popSize = GetParameter(_POP_SIZE);
    RbtBool bScale = GetParameter(_SCALE_CHROM_LENGTH);
    if (bScale) {
        RbtInt chromLength = m_chrom->GetLength();
        popSize *= chromLength;
    } 
    if (GetTrace() > 3) {
        cout << _CT << ": popSize=" << popSize << endl;
    }
    RbtPopulationPtr pop = new RbtPopulation(m_chrom, popSize, pSF);
    pop->Best()->GetChrom()->SyncToModel();
    GetWorkSpace()->SetPopulation(pop);
}
