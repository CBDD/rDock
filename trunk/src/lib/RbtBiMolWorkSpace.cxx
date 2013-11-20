/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtBiMolWorkSpace.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtBiMolWorkSpace.h"
#include "RbtBaseSF.h"
#include "RbtAnnotationHandler.h"
#include "RbtBaseMolecularFileSource.h"

using namespace std;

////////////////////////////////////////
//Constructors/destructors

RbtBiMolWorkSpace::RbtBiMolWorkSpace() : RbtWorkSpace(2) {
  _RBTOBJECTCOUNTER_CONSTR_("RbtBiMolWorkSpace");
}

RbtBiMolWorkSpace::~RbtBiMolWorkSpace() {
  _RBTOBJECTCOUNTER_DESTR_("RbtBiMolWorkSpace");
}
	
////////////////////////////////////////
//Public methods
////////////////

//Model handling
RbtModelPtr RbtBiMolWorkSpace::GetReceptor() const {return GetModel(RECEPTOR);}
RbtModelPtr RbtBiMolWorkSpace::GetLigand() const {return GetModel(LIGAND);}
RbtModelList RbtBiMolWorkSpace::GetSolvent() const {return GetModels(SOLVENT);}
RbtBool RbtBiMolWorkSpace::hasSolvent() const {
	return GetNumModels() > SOLVENT;
}

void RbtBiMolWorkSpace::SetReceptor(RbtModelPtr spReceptor) {
  SetModel(RECEPTOR,spReceptor);
  //DM 30 Oct 2001 - no need to clear the ligand population, there may be scenarios where we want
  //to continue a ligand docking with a revised receptor model
}

void RbtBiMolWorkSpace::SetLigand(RbtModelPtr spLigand) {
 //DM 30 Oct 2001 - invalidate the current population
  ClearPopulation();
  //DM 26 Aug 2002 - better to call SetModel AFTER defining the ligand as flexible
  SetModel(LIGAND,spLigand);
}


void RbtBiMolWorkSpace::SetSolvent(RbtModelList solventList) {
	if (hasSolvent()) {
		SetModels(SOLVENT,solventList);
	}
	else {
		AddModels(solventList);
	}
}

void RbtBiMolWorkSpace::RemoveSolvent() {
	if (hasSolvent()) {
		RemoveModels(SOLVENT);
	}
}

void RbtBiMolWorkSpace::UpdateModelCoordsFromChromRecords(RbtBaseMolecularFileSource* pSource,
                                            RbtInt iTrace) {
    RbtInt nModels = GetNumModels();
    for (RbtInt iModel = 0; iModel < nModels; iModel++) {
        //if (iModel == LIGAND) continue;
        RbtModelPtr spModel = GetModel(iModel);
        if (spModel.Ptr()) {
            RbtChromElementPtr spChrom = spModel->GetChrom();
            if (spChrom.Ptr()) {
                RbtInt chromLength = spChrom->GetLength();
                if (chromLength > 0) {
                    ostringstream ostr;
                    ostr << "CHROM." << iModel;
                    RbtString chromField = ostr.str();
                    if (pSource->isDataFieldPresent(chromField)) {
                        //TODO: Move this code to RbtVariant class
                        //Concatenate the multi-record value into a single string
                        RbtString chromRecord = Rbt::ConvertListToDelimitedString(
                                                        pSource->GetDataValue(chromField), ",");
                        //Now split into string values and convert to doubles
                         RbtStringList chromValues = Rbt::ConvertDelimitedStringToList(
                                                        chromRecord, ",");
                         if (chromLength == chromValues.size()) {
                            RbtDoubleList chromVec;
                            for (RbtStringListConstIter iter = chromValues.begin();
                                                        iter != chromValues.end();
                                                        ++iter) {
                                RbtString chromValue(*iter);
                                chromVec.push_back(atof(chromValue.c_str()));
                            }
                            spChrom->SetVector(chromVec);
                            spChrom->SyncToModel();
                            spModel->UpdatePseudoAtoms();
                         }
                         else if (iTrace > 0) {
                            cout << "Mismatched chromosome sizes for model #" << iModel << endl;
                            cout << chromField << " record in " << pSource->GetFileName() \
                                 << " has " << chromValues << " elements" << endl;
                            cout << "Expected number of elements is " << chromLength << endl;
                            cout << "Model chromosome not updated" << endl;
                         }
                    }
                }
            }
        }
    }
}

//Model I/O
//Saves ligand to file sink
void RbtBiMolWorkSpace::Save(RbtBool bSaveScores) {
  SaveLigand(GetSink(),bSaveScores);
}

//Model I/O
//Saves ligand to history file sink
void RbtBiMolWorkSpace::SaveHistory(RbtBool bSaveScores) {
  SaveLigand(GetHistorySink(),bSaveScores);
}

//Private method to save ligand to file sink, optionally with score attached
void RbtBiMolWorkSpace::SaveLigand(RbtMolecularFileSinkPtr spSink, RbtBool bSaveScores) {
  if (spSink.Null()) //Check we have a valid sink
    return;
  RbtModelPtr spLigand(GetLigand());
  if (spLigand.Null()) //Check we have a ligand
    return;
  
  RbtBaseSF* pSF(GetSF());
  //If we have a scoring function, clear any current score data
  if (pSF) {
    RbtString strSFName(pSF->GetFullName());
    spLigand->ClearAllDataFields(strSFName);
    spLigand->ClearAllDataFields(RbtAnnotationHandler::_ANNOTATION_FIELD);
  }
  //Save current score components as model data fields
  if (bSaveScores && pSF) {
    RbtStringVariantMap scoreMap;
    pSF->ScoreMap(scoreMap);
    for (RbtStringVariantMapConstIter vIter = scoreMap.begin(); vIter != scoreMap.end(); vIter++) {
      spLigand->SetDataValue((*vIter).first,(*vIter).second);
    }
    //Save the chromosome values for all models for later retrieval
    spLigand->ClearAllDataFields("CHROM.");
    RbtInt nModels = GetNumModels();
    for (RbtInt iModel = 0; iModel < nModels; iModel++) {
        RbtModelPtr spModel = GetModel(iModel);
        if (spModel.Ptr()) {
            RbtChromElementPtr spChrom = spModel->GetChrom();
            if (spChrom.Ptr() && (spChrom->GetLength() > 0) ) {
                spChrom->SyncFromModel();
                RbtDoubleList chromVec;
                spChrom->GetVector(chromVec);
                ostringstream ostr;
                ostr << "CHROM." << iModel;
                spLigand->SetDataValue(ostr.str(),RbtVariant(chromVec,72,8));
            }
        }
    }
  }
  spSink->SetModel(spLigand);
  if (hasSolvent()) {
  	spSink->SetSolvent(GetSolvent());
  }
  else {
  	spSink->SetSolvent(RbtModelList());
  }
  spSink->Render();
}

