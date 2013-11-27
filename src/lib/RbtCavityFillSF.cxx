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

#include "RbtCavityFillSF.h"
#include "RbtDockingSite.h"
#include "RbtWorkSpace.h"

//Static data members
RbtString RbtCavityFillSF::_CT("RbtCavityFillSF");
	
//NB - Virtual base class constructor (RbtBaseSF) gets called first,
//implicit constructor for RbtBaseInterSF is called second
RbtCavityFillSF::RbtCavityFillSF(const RbtString& strName)
  : RbtBaseSF(_CT,strName) {
  //Add parameters
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtCavityFillSF::~RbtCavityFillSF() {
#ifdef _DEBUG
	cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}
	

void RbtCavityFillSF::SetupReceptor() {
  m_spGrid = RbtFFTGridPtr();
  RbtDockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
  if (spDS.Null())
    return;

  RbtInt iTrace = GetTrace();

  //Recreate the cavity grid
  RbtCavityList cavList = spDS->GetCavityList();
  if (cavList.empty())
    return;
  RbtVector gridStep = cavList.front()->GetGridStep();
  RbtDouble border = 10.0;
  RbtCoord minCoord = spDS->GetMinCoord()-border;
  RbtCoord maxCoord = spDS->GetMaxCoord()+border;
  RbtVector recepExtent = maxCoord-minCoord;
  RbtUInt nX = int(recepExtent.x/gridStep.x)+1;
  RbtUInt nY = int(recepExtent.y/gridStep.y)+1;
  RbtUInt nZ = int(recepExtent.z/gridStep.z)+1;
  m_spGrid = new RbtFFTGrid(minCoord,gridStep,nX,nY,nZ);

  RbtAtomList atomList = GetReceptor()->GetAtomList();
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    RbtDouble r = (**iter).GetVdwRadius();
    m_spGrid->SetSphere((**iter).GetCoords(),r+0.3,-1.0,true);
  }

  if (iTrace > 1) {
    cout << endl << "EXCLUDE RECEPTOR VOLUME" << endl;
    cout << "N(excluded)=" << m_spGrid->Count(-1.0) << endl;
    cout << "N(unallocated)=" << m_spGrid->Count(0.0) << endl;
  }
  
  for (RbtCavityListConstIter iter = cavList.begin(); iter != cavList.end(); ++iter) {
    const RbtCoordList coordList = (*iter)->GetCoordList();
    for (RbtCoordListConstIter iter2 = coordList.begin(); iter2 != coordList.end(); ++iter2) {
      m_spGrid->SetValue(*iter2,1.0);
    }
  }

  if (iTrace > 1) {
    cout << endl << "DEFINE CAVITY VOLUME" << endl;
    cout << "N(excluded)=" << m_spGrid->Count(-1.0) << endl;
    cout << "N(cavity)=" << m_spGrid->Count(1.0) << endl;
    cout << "N(unallocated)=" << m_spGrid->Count(0.0) << endl;
  }

  //RbtDouble largeR = 4.0;
  //Map with a large solvent sphere
  //m_spGrid->SetAccessible(largeR,1.0,-1.0,0.0,false);
  //if (iTrace > 1) {
  //  cout << endl << "MAP WITH LARGE SPHERE" << endl;
  //  cout << "N(excluded)=" << m_spGrid->Count(-1.0) << endl;
  //  cout << "N(cavity)=" << m_spGrid->Count(1.0) << endl;
  //  cout << "N(unallocated)=" << m_spGrid->Count(0.0) << endl;
  //}
}

void RbtCavityFillSF::SetupLigand() {
  m_ligAtomList.clear();
  if (GetLigand().Null())
    return;

  m_ligAtomList = GetLigand()->GetAtomList();
}

void RbtCavityFillSF::SetupScore() {

}

RbtDouble RbtCavityFillSF::RawScore() const {
  //Check grid is defined
  if (m_spGrid.Null())
    return 0.0;

  RbtInt iTrace = GetTrace();
  RbtFFTGridPtr gridCopy = new RbtFFTGrid(*m_spGrid);
  
  for (RbtAtomListConstIter iter = m_ligAtomList.begin(); iter != m_ligAtomList.end(); ++iter) {
    RbtDouble r = (*iter)->GetVdwRadius();
    gridCopy->SetSphere((*iter)->GetCoords(),r+0.3,-1.0,true);
  }

  if (iTrace > 1) {
    cout << endl << "EXCLUDE LIGAND VOLUME" << endl;
    cout << "N(excluded)=" << gridCopy->Count(-1.0) << endl;
    cout << "N(cavity)=" << gridCopy->Count(1.0) << endl;
    cout << "N(unallocated)=" << gridCopy->Count(0.0) << endl;
  }

  //Map with a small solvent sphere
  gridCopy->SetAccessible(2.0,1.0,-1.0,0.0,false);
  gridCopy->SetAccessible(1.5,1.0,-1.0,2.0,false);

  if (iTrace > 1) {
    cout << endl << "VOID DETECTION" << endl;
    cout << "N(excluded)=" << gridCopy->Count(-1.0) << endl;
    cout << "N(cavity)=" << gridCopy->Count(1.0) << endl;
    cout << "N(voids)=" << gridCopy->Count(2.0) << endl;
    cout << "N(unallocated)=" << gridCopy->Count(0.0) << endl;
  }

  gridCopy->ReplaceValue(-1.0,0.0);
  gridCopy->ReplaceValue(1.0,0.0);
  ofstream dumpFile("void.grd");
  gridCopy->PrintInsightGrid(dumpFile);
  dumpFile.close();

  //Find the contiguous regions of cavity grid points
  RbtDouble minSize = 50;
  RbtFFTPeakMap peakMap = gridCopy->FindPeaks(2.0,minSize);

  //Convert peaks to cavity format
  RbtCavityList cavityList;
  RbtVector gridStep = gridCopy->GetGridStep();
  for (RbtFFTPeakMapConstIter pIter = peakMap.begin(); pIter != peakMap.end(); pIter++) {
    RbtFFTPeakPtr spPeak((*pIter).second);
    RbtCoordList coordList = gridCopy->GetCoordList(spPeak->points);
    RbtCavityPtr spCavity = RbtCavityPtr(new RbtCavity(coordList,gridStep));
    cavityList.push_back(spCavity);
  }

  //Sort cavities by volume
  std::sort(cavityList.begin(),cavityList.end(),Rbt::RbtCavityPtrCmp_Volume());

  if (iTrace > 0) {
    for (RbtCavityListConstIter cIter = cavityList.begin(); cIter != cavityList.end(); cIter++) {
      cout << (**cIter) << endl;
    }
  }

  return 0.0;
}

//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtCavityFillSF::ParameterUpdated(const RbtString& strName) {
  RbtBaseSF::ParameterUpdated(strName);
}

