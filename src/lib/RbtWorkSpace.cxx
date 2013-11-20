/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtWorkSpace.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtWorkSpace.h"
#include "RbtBaseSF.h"
#include "RbtBaseTransform.h"

//Static data members
RbtString RbtWorkSpace::_CT("RbtWorkSpace");
RbtString RbtWorkSpace::_NAME("NAME");

////////////////////////////////////////
//Constructors/destructors

//Create an empty model container of the right size
RbtWorkSpace::RbtWorkSpace(RbtUInt nModels) : m_models(nModels),m_SF(NULL),m_transform(NULL) {
	AddParameter(_NAME,_CT);
#ifdef _DEBUG
	cout << "RbtWorkSpace::RbtWorkSpace(): Created model list of size " << m_models.size() << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_("RbtWorkSpace");
}

RbtWorkSpace::~RbtWorkSpace() {
  _RBTOBJECTCOUNTER_DESTR_("RbtWorkSpace");
}

////////////////////////////////////////
//Public methods
////////////////

//Get/set workspace name
RbtString RbtWorkSpace::GetName() const {return GetParameter(_NAME);}
void RbtWorkSpace::SetName(const RbtString& strName) {SetParameter(_NAME,strName);}
		
//Model handling
		
//Returns number of models in workspace
RbtUInt RbtWorkSpace::GetNumModels() const {return m_models.size();}

//Returns vector of all models in workspace
RbtModelList RbtWorkSpace::GetModels() const {return m_models;}

//Returns a specific (numbered) model
//Throws RbtBadArgument if iModel out of range
RbtModelPtr RbtWorkSpace::GetModel(RbtUInt iModel) const throw (RbtError)
{
	if (iModel >= m_models.size()) {
		throw RbtBadArgument(_WHERE_,"iModel out of range");
	}
	return m_models[iModel];
}

//Replace an existing model, notify all observers
//Throws RbtBadArgument if iModel out of range
void RbtWorkSpace::SetModel(RbtUInt iModel,RbtModelPtr spModel) throw (RbtError)
{
	if (iModel >= m_models.size()) {
		throw RbtBadArgument(_WHERE_,"iModel out of range");
	}
	m_models[iModel] = spModel;
	Notify();//Notify observers of change in state
}

//Returns vector of models, starting from index iModel
RbtModelList RbtWorkSpace::GetModels(RbtUInt iModel) const throw (RbtError)
{
	RbtModelList retVal;
	if (iModel >= m_models.size()) {
		throw RbtBadArgument(_WHERE_,"iModel out of range");
	}
	std::copy(m_models.begin()+iModel, m_models.end(), std::back_inserter(retVal));
	return retVal;
}

//Append a number of models to the workspace, increasing the total number of models
void RbtWorkSpace::AddModels(RbtModelList modelList)
{
	std::copy(modelList.begin(),modelList.end(),std::back_inserter(m_models));
	Notify();//Notify observers of change in state
}

//Replace a number of existing models
//iModel is the index of the first model to replace
//Throws RbtBadArgument if iModel out of range or modelList too large
void RbtWorkSpace::SetModels(RbtUInt iModel, RbtModelList modelList) throw (RbtError)
{
	if (iModel > (m_models.size() - modelList.size())) {
		throw RbtBadArgument(_WHERE_,"iModel out of range");
	}
	else if (modelList.size() > m_models.size()) {
		throw RbtBadArgument(_WHERE_,"modelList too large");
	}
	std::copy(modelList.begin(),modelList.end(),m_models.begin()+iModel);	
	Notify();//Notify observers of change in state
}

//Removes a number of models from the workspace
//Removes from index iModel to end of model list
void RbtWorkSpace::RemoveModels(RbtUInt iModel) throw (RbtError) {
	if (iModel >= m_models.size()) {
		throw RbtBadArgument(_WHERE_,"iModel out of range");
	}
	m_models.erase(m_models.begin()+iModel,m_models.end());
	Notify();
}

//Model I/O
//Get/set the molecular file sink (for outputting ligands)
RbtMolecularFileSinkPtr RbtWorkSpace::GetSink() const {return m_spSink;}
void RbtWorkSpace::SetSink(RbtMolecularFileSinkPtr spSink) {m_spSink = spSink;}

//Saves models to file sink
//Base workspace does nothing
//It is up to subclasses to decide what action to take
void RbtWorkSpace::Save(RbtBool bSaveScores) {}

//Get/set the history file sink
RbtMolecularFileSinkPtr RbtWorkSpace::GetHistorySink() const {return m_spHisSink;}
void RbtWorkSpace::SetHistorySink(RbtMolecularFileSinkPtr spHisSink) {m_spHisSink = spHisSink;}
//Saves models to file sink
//Base workspace does nothing
//It is up to subclasses to decide what action to take
void RbtWorkSpace::SaveHistory(RbtBool bSaveScores) {}

//Scoring function handling
//SetSF automatically registers the scoring function with the workspace
RbtBaseSF* RbtWorkSpace::GetSF() const {return m_SF;}
void RbtWorkSpace::SetSF(RbtBaseSF* pSF) {
	if (m_SF) {
		m_SF->Unregister();
	}
	m_SF = pSF;
	if (m_SF) {
		m_SF->Register(this);
        Notify();
	}
}

//Get/Set the transform
//SetTransform automatically registers the transform with the workspace
RbtBaseTransform* RbtWorkSpace::GetTransform() const {return m_transform;}
void RbtWorkSpace::SetTransform(RbtBaseTransform* pTransform) {
	if (m_transform) {
		m_transform->Unregister();
	}
	m_transform = pTransform;
	if (m_transform) {
		m_transform->Register(this);
        Notify();
	}
}

//Run the simulation!
void RbtWorkSpace::Run() {
	if (m_transform) {
		m_transform->Go();
	}
}

//Population handling
//DM 30 Oct 2001 - provide a home for persistent populations
void RbtWorkSpace::SetPopulation(RbtPopulationPtr population) {
  m_population = population;
}

RbtPopulationPtr RbtWorkSpace::GetPopulation() const {
  return m_population;
}

void RbtWorkSpace::ClearPopulation() {
    m_population.SetNull();
}


//Docking site handling
//DM 09 Apr 2002 - workspace now manages the docking site
RbtDockingSitePtr RbtWorkSpace::GetDockingSite() const {
  return m_spDockSite;
}

void RbtWorkSpace::SetDockingSite(RbtDockingSitePtr spDockSite) {
  m_spDockSite = spDockSite;
}

//Filter handling
//BGD 27 Jan 2003 - SetFilter automotically registers the Filter 
//with the workspace
RbtFilterPtr RbtWorkSpace::GetFilter() const {
  return m_spFilter;
}

void RbtWorkSpace::SetFilter(RbtFilterPtr spFilter) 
{
	if (m_spFilter)
	{
		m_spFilter->Unregister();
	}
  m_spFilter = spFilter;
	if (m_spFilter)
	{
		m_spFilter->Register(this);
	}
}
