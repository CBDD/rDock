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

#include "RbtSphereSiteMapper.h"
#include "RbtFFTGrid.h"

//Static data member for class type
RbtString RbtSphereSiteMapper::_CT("RbtSphereSiteMapper");
RbtString RbtSphereSiteMapper::_VOL_INCR("VOL_INCR");
RbtString RbtSphereSiteMapper::_SMALL_SPHERE("SMALL_SPHERE");
RbtString RbtSphereSiteMapper::_LARGE_SPHERE("LARGE_SPHERE");
RbtString RbtSphereSiteMapper::_GRIDSTEP("GRIDSTEP");
RbtString RbtSphereSiteMapper::_CENTER("CENTER");
RbtString RbtSphereSiteMapper::_RADIUS("RADIUS");
RbtString RbtSphereSiteMapper::_MIN_VOLUME("MIN_VOLUME");
RbtString RbtSphereSiteMapper::_MAX_CAVITIES("MAX_CAVITIES");

RbtSphereSiteMapper::RbtSphereSiteMapper(const RbtString& strName) : RbtSiteMapper(_CT,strName) {
  //Add parameters
  AddParameter(_VOL_INCR,0.0);
  AddParameter(_SMALL_SPHERE,1.5);
  AddParameter(_LARGE_SPHERE,4.0);
  AddParameter(_GRIDSTEP,0.5);
  AddParameter(_CENTER,RbtCoord());
  AddParameter(_RADIUS,10.0);
  AddParameter(_MIN_VOLUME,100.0);//Min cavity volume in A^3
  AddParameter(_MAX_CAVITIES,99);//Max number of cavities to return
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSphereSiteMapper::~RbtSphereSiteMapper() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

RbtCavityList RbtSphereSiteMapper::operator() () {
  RbtCavityList cavityList;
  RbtModelPtr spReceptor(GetReceptor());
  if (spReceptor.Null()) return cavityList;

  RbtDouble rVolIncr = GetParameter(_VOL_INCR);
  RbtDouble smallR = GetParameter(_SMALL_SPHERE);
  RbtDouble largeR = GetParameter(_LARGE_SPHERE);
  RbtDouble step = GetParameter(_GRIDSTEP);
  RbtCoord center = GetParameter(_CENTER);
  RbtDouble radius = GetParameter(_RADIUS);
  RbtDouble minVol = GetParameter(_MIN_VOLUME);
  RbtInt maxCavities = GetParameter(_MAX_CAVITIES);
  RbtInt iTrace = GetTrace();

  //Grid values
  const RbtDouble recVal = -1.0;//Receptor volume
  const RbtDouble larVal = -0.75;//Accessible to large sphere
  const RbtDouble excVal = -0.5;//Excluded from calculation
  const RbtDouble borVal = -0.25;//Border region (used for mapping large sphere only)
  const RbtDouble cavVal = 1.0;//Cavities

  //Convert from min volume (in A^3) to min size (number of grid points)
  RbtDouble minSize = minVol/(step*step*step);
  RbtFFTGridPtr spReceptorGrid;
  //Only include non-H receptor atoms in the mapping
  RbtAtomList atomList = Rbt::GetAtomList(spReceptor->GetAtomList(),std::not1(Rbt::isAtomicNo_eq(1)));
  RbtVector gridStep(step,step,step);

  //We extend the grid by 2*largeR on each side to eliminate edge effects in the cavity mapping
  RbtCoord minCoord=center-radius;
  RbtCoord maxCoord=center+radius;
  RbtDouble border = 2.0*(largeR + step);
  minCoord -= border;
  maxCoord += border;
  RbtVector recepExtent = maxCoord-minCoord;
  RbtUInt nX = int(recepExtent.x/gridStep.x)+1;
  RbtUInt nY = int(recepExtent.y/gridStep.y)+1;
  RbtUInt nZ = int(recepExtent.z/gridStep.z)+1;
  spReceptorGrid = RbtFFTGridPtr(new RbtFFTGrid(minCoord,gridStep,nX,nY,nZ));
  center = spReceptorGrid->GetGridCenter();

  //Initialise the grid with a zero-value region in the user-specified sphere
  //surrounded by a border region of thickness = large sphere radius
  spReceptorGrid->SetAllValues(excVal);
  spReceptorGrid->SetSphere(center,radius+largeR,borVal,true);
  spReceptorGrid->SetSphere(center,radius,0.0,true);
  if (iTrace > 1) {
    cout << endl << "INITIALISATION" << endl;
    cout << "Center=" << center << endl;
    cout << "Radius=" << radius << endl;
    cout << "Border=" << border << endl;
    cout << "N(excluded)=" << spReceptorGrid->Count(excVal) << endl;
    cout << "N(border)=" << spReceptorGrid->Count(borVal) << endl;
    cout << "N(unallocated)=" << spReceptorGrid->Count(0.0) << endl;
  }

  //Set all vdW volume grid points to the value -1
  //Add increment to all vdw radii
  //Iterate over all receptor atoms as we want to include those whose centers are outside
  //the active site, but whose vdw volumes overlap the active site region.
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    RbtDouble r = (**iter).GetVdwRadius();
    spReceptorGrid->SetSphere((**iter).GetCoords(),r+rVolIncr,recVal,true);
  }

  if (iTrace > 1) {
    cout << endl << "EXCLUDE RECEPTOR VOLUME" << endl;
    cout << "N(receptor)=" << spReceptorGrid->Count(recVal) << endl;
    cout << "N(excluded)=" << spReceptorGrid->Count(excVal) << endl;
    cout << "N(border)=" << spReceptorGrid->Count(borVal) << endl;
    cout << "N(unallocated)=" << spReceptorGrid->Count(0.0) << endl;
  }
  
  //Now map the solvent accessible regions with a large sphere
  //We first map the border region, which will also sweep out and exclude regions of the user-specified inner region
  //This is the first key step for preventing edge effects. 
  spReceptorGrid->SetAccessible(largeR,borVal,recVal,larVal,false);
  if (iTrace > 1) {
    cout << endl << "EXCLUDE LARGE SPHERE (Border region)" << endl;
    cout << "N(receptor)=" << spReceptorGrid->Count(recVal) << endl;
    cout << "N(large sphere)=" << spReceptorGrid->Count(larVal) << endl;
    cout << "N(excluded)=" << spReceptorGrid->Count(excVal) << endl;
    cout << "N(border)=" << spReceptorGrid->Count(borVal) << endl;
    cout << "N(unallocated)=" << spReceptorGrid->Count(0.0) << endl;
  }
  spReceptorGrid->SetAccessible(largeR,0.0,recVal,larVal,false);
  if (iTrace > 1) {
    cout << endl << "EXCLUDE LARGE SPHERE (Unallocated inner region)" << endl;
    cout << "N(receptor)=" << spReceptorGrid->Count(recVal) << endl;
    cout << "N(large sphere)=" << spReceptorGrid->Count(larVal) << endl;
    cout << "N(excluded)=" << spReceptorGrid->Count(excVal) << endl;
    cout << "N(border)=" << spReceptorGrid->Count(borVal) << endl;
    cout << "N(unallocated)=" << spReceptorGrid->Count(0.0) << endl;
  }

  //Finally with a smaller radius. This is the region we want to search for peaks in,
  //so set to a positive value
  //But first we need to replace all non-zero values with the receptor volume value
  //otherwise SetAccessible will not work properly. This is the other key step
  //for preventing edge effects
  spReceptorGrid->ReplaceValue(borVal,recVal);
  spReceptorGrid->ReplaceValue(excVal,recVal);
  spReceptorGrid->ReplaceValue(larVal,recVal);
  spReceptorGrid->SetAccessible(smallR,0.0,recVal,cavVal,false);

  if (iTrace > 1) {
    cout << endl << "FINAL CAVITIES" << endl;
    cout << "N(receptor)=" << spReceptorGrid->Count(recVal) << endl;
    cout << "N(large sphere)=" << spReceptorGrid->Count(larVal) << endl;
    cout << "N(excluded)=" << spReceptorGrid->Count(excVal) << endl;
    cout << "N(border)=" << spReceptorGrid->Count(borVal) << endl;
    cout << "N(unallocated)=" << spReceptorGrid->Count(0.0) << endl;
    cout << "N(cavities)=" << spReceptorGrid->Count(cavVal) << endl;
    cout << endl << "Min cavity size=" << minSize << endl;
  }

  //Find the contiguous regions of cavity grid points
  RbtFFTPeakMap peakMap = spReceptorGrid->FindPeaks(cavVal,minSize);

  //Convert peaks to cavity format
  for (RbtFFTPeakMapConstIter pIter = peakMap.begin(); pIter != peakMap.end(); pIter++) {
    RbtFFTPeakPtr spPeak((*pIter).second);
    RbtCoordList coordList = spReceptorGrid->GetCoordList(spPeak->points);
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
