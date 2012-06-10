/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtRandPopTransform.cxx#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

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
