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

#include "RbtBaseMolecularFileSink.h"

////////////////////////////////////////
//Constructors/destructors
//RbtBaseMolecularFileSink::RbtBaseMolecularFileSink(const char* fileName, RbtModelPtr spModel) :
//  RbtBaseFileSink(fileName), m_spModel(spModel), m_bMultiConf(false)
//{
//  Reset();
//  _RBTOBJECTCOUNTER_CONSTR_("RbtBaseMolecularFileSink");  
//}

//DM 21 Apr 1999 - Default is to initialise m_bUseModelSegmentNames to true in the constructor
//so that the first model is rendered using the model segment names
//Default in calls to SetModel is to set m_bUseModelSegmentNames to false
//so that subsequent Render()s use numeric segment IDs
RbtBaseMolecularFileSink::RbtBaseMolecularFileSink(const RbtString& fileName, RbtModelPtr spModel,
						   RbtBool bUseModelSegmentNames) :
  RbtBaseFileSink(fileName), m_spModel(spModel), m_bMultiConf(false), m_bUseModelSegmentNames(bUseModelSegmentNames)
{
  Reset();
  _RBTOBJECTCOUNTER_CONSTR_("RbtBaseMolecularFileSink");  
}

//Default destructor
RbtBaseMolecularFileSink::~RbtBaseMolecularFileSink()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtBaseMolecularFileSink");
}


////////////////////////////////////////
//Public methods
////////////////
  
//Get/set the RbtModel that is linked to the sink
RbtModelPtr RbtBaseMolecularFileSink::GetModel() const
{
  return m_spModel;
}

void RbtBaseMolecularFileSink::SetModel(RbtModelPtr spModel, RbtBool bUseModelSegmentNames)
{
  m_spModel = spModel;
  m_bUseModelSegmentNames = bUseModelSegmentNames;
}

//Get/set the multiconformer status
RbtBool RbtBaseMolecularFileSink::GetMultiConf() const
{
  return m_bMultiConf;
}

void RbtBaseMolecularFileSink::SetMultiConf(RbtBool bMultiConf)
{
  //If the MultiConf status of the sink is changing, we need to ensure
  //that the cache gets written first
  if (m_bMultiConf != bMultiConf) {
    Write();
    Reset();
  }
  m_bMultiConf = bMultiConf;
}

//Force writing the file following several multi-conf Renders
void RbtBaseMolecularFileSink::WriteMultiConf()
{
  Write();
  Reset();
}

RbtModelList RbtBaseMolecularFileSink::GetSolvent() const {
	return m_solventList;
}

void RbtBaseMolecularFileSink::SetSolvent(RbtModelList solventList) {
	m_solventList = solventList;
}

////////////////////////////////////////
//Protected methods
///////////////////
//Reset the atom, residue and segment numbering
void RbtBaseMolecularFileSink::Reset()
{
  if ( (!m_bMultiConf) || isCacheEmpty()) {
    m_nAtomId = 0;
    m_nSubunitId = 0;
    m_nSegmentId = 0;
  }
}

////////////////////////////////////////
//Private methods
/////////////////
