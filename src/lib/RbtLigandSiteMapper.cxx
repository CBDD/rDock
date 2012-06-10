/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtLigandSiteMapper.cxx#4 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtLigandSiteMapper.h"
#include "RbtFFTGrid.h"
#include "RbtMdlFileSource.h"
#include "RbtWorkSpace.h"

//Static data member for class type
RbtString RbtLigandSiteMapper::_CT("RbtLigandSiteMapper");
RbtString RbtLigandSiteMapper::_REF_MOL("REF_MOL");
RbtString RbtLigandSiteMapper::_VOL_INCR("VOL_INCR");
RbtString RbtLigandSiteMapper::_SMALL_SPHERE("SMALL_SPHERE");
RbtString RbtLigandSiteMapper::_GRIDSTEP("GRIDSTEP");
RbtString RbtLigandSiteMapper::_RADIUS("RADIUS");
RbtString RbtLigandSiteMapper::_MIN_VOLUME("MIN_VOLUME");
RbtString RbtLigandSiteMapper::_MAX_CAVITIES("MAX_CAVITIES");

RbtLigandSiteMapper::RbtLigandSiteMapper(const RbtString& strName) : RbtSiteMapper(_CT,strName) {
  //Add parameters
  AddParameter(_REF_MOL,"ref.sd");
  AddParameter(_VOL_INCR,0.0);
  AddParameter(_SMALL_SPHERE,1.5);
  AddParameter(_GRIDSTEP,0.5);
  AddParameter(_RADIUS,10.0);
  AddParameter(_MIN_VOLUME,100.0);//Min cavity volume in A^3
  AddParameter(_MAX_CAVITIES,99);//Max number of cavities to return
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtLigandSiteMapper::~RbtLigandSiteMapper() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

RbtCavityList RbtLigandSiteMapper::operator() () {
  RbtCavityList cavityList;
  RbtModelPtr spReceptor(GetReceptor());
  if (spReceptor.Null()) return cavityList;

  RbtString strRef = GetParameter(_REF_MOL);
  RbtDouble rVolIncr = GetParameter(_VOL_INCR);
  RbtDouble smallR = GetParameter(_SMALL_SPHERE);
  RbtDouble step = GetParameter(_GRIDSTEP);
  RbtDouble radius = GetParameter(_RADIUS);
  RbtDouble minVol = GetParameter(_MIN_VOLUME);
  RbtInt maxCavities = GetParameter(_MAX_CAVITIES);
  RbtInt iTrace = GetTrace();

  //Grid values
  const RbtDouble recVal = -1.0;//Receptor volume
  const RbtDouble cavVal = 1.0;//Ligand expanded volume

  //Convert from min volume (in A^3) to min size (number of grid points)
  RbtDouble minSize = minVol/(step*step*step);

  //Get the reference atom list from the mol/SD file
  //Simple extension would be to read all the records from the SD file
  //to construct a docking site from the superimposition of all the reference ligands
  RbtString strMolFile = Rbt::GetRbtFileName("data/ligands",strRef);
  RbtMolecularFileSourcePtr spMdlFileSource(new RbtMdlFileSource(strMolFile,false,false,true));
  //Only include non-H reference atoms in the mapping
  RbtAtomList refAtomList = Rbt::GetAtomList(spMdlFileSource->GetAtomList(),std::not1(Rbt::isAtomicNo_eq(1)));

  RbtFFTGridPtr spGrid;

  //Only include non-H receptor atoms in the mapping
  RbtAtomList atomList = Rbt::GetAtomList(spReceptor->GetAtomList(),std::not1(Rbt::isAtomicNo_eq(1)));
  RbtVector gridStep(step,step,step);

  //Construct the grid to cover the ligand coords + radius
  RbtCoordList refCoordList;
  Rbt::GetCoordList(refAtomList,refCoordList);
  RbtDouble border = radius+smallR+step;
  RbtCoord minCoord = Rbt::Min(refCoordList) - border;
  RbtCoord maxCoord = Rbt::Max(refCoordList) + border;
  RbtVector extent = maxCoord-minCoord;
  RbtUInt nX = int(extent.x/gridStep.x)+1;
  RbtUInt nY = int(extent.y/gridStep.y)+1;
  RbtUInt nZ = int(extent.z/gridStep.z)+1;
  spGrid = RbtFFTGridPtr(new RbtFFTGrid(minCoord,gridStep,nX,nY,nZ));
  spGrid->SetAllValues(recVal);

  //Clear the spheres around each ligand atom
  for (RbtCoordListConstIter iter = refCoordList.begin(); iter != refCoordList.end(); iter++) {
    spGrid->SetSphere((*iter),radius,0.0,true);
  }

  if (iTrace > 1) {
    cout << endl << "INITIALISATION" << endl;
    cout << "Radius=" << radius << endl;
    cout << "Border=" << border << endl;
    cout << "N(excluded)=" << spGrid->Count(recVal) << endl;
    cout << "N(unallocated)=" << spGrid->Count(0.0) << endl;
  }

  //Now exclude the receptor volume
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    RbtDouble r = (**iter).GetVdwRadius();
    spGrid->SetSphere((**iter).GetCoords(),r+rVolIncr,recVal,true);
  }

  if (iTrace > 1) {
    cout << endl << "EXCLUDE RECEPTOR VOLUME" << endl;
    cout << "N(excluded)=" << spGrid->Count(recVal) << endl;
    cout << "N(unallocated)=" << spGrid->Count(0.0) << endl;
  }

  //Map with a small solvent sphere
  spGrid->SetAccessible(smallR,0.0,recVal,cavVal,false);

  if (iTrace > 1) {
    cout << endl << "FINAL CAVITIES" << endl;
    cout << "N(excluded)=" << spGrid->Count(recVal) << endl;
    cout << "N(unallocated)=" << spGrid->Count(0.0) << endl;
    cout << "N(cavities)=" << spGrid->Count(cavVal) << endl;
    cout << endl << "Min cavity size=" << minSize << endl;
  }

  //Find the contiguous regions of cavity grid points
  RbtFFTPeakMap peakMap = spGrid->FindPeaks(cavVal,minSize);

  //Convert peaks to cavity format
  for (RbtFFTPeakMapConstIter pIter = peakMap.begin(); pIter != peakMap.end(); pIter++) {
    RbtFFTPeakPtr spPeak((*pIter).second);
    RbtCoordList coordList = spGrid->GetCoordList(spPeak->points);
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

  //Limit the number of cavities if necessary
  if (cavityList.size() > maxCavities) {
    if (iTrace > 0) {
      cout << endl << cavityList.size() << " cavities identified - limit to " << maxCavities << " largest cavities" << endl;
    }
    cavityList.erase(cavityList.begin()+maxCavities,cavityList.end());
  }

  return cavityList;
}
