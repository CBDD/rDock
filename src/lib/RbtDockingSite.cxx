/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtDockingSite.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtDockingSite.h"
#include "RbtFileError.h"
#include <cstring> 

//Less than operator for sorting coords
class RbtCoordCmp {
public:
  RbtBool operator()(const RbtCoord& c1, const RbtCoord& c2) const {
    if (c1.x < c2.x)
      return true;
    else if (c1.x == c2.x) {
      if (c1.y < c2.y)
	return true;
      else if (c1.y == c2.y) {
	if (c1.z < c2.z)
	  return true;
      }
    }
    return false;
  }
};

//Static data members
RbtString RbtDockingSite::_CT("RbtDockingSite");

//STL predicate for selecting atoms within a defined distance range from nearest cavity coords
//Uses precalculated distance grid
RbtBool RbtDockingSite::isAtomInRange::operator() (RbtAtom* pAtom) const {
  const RbtCoord& c = pAtom->GetCoords();
  if (!m_pGrid->isValid(c)) {
    return false;
  }
  RbtDouble dist = m_pGrid->GetSmoothedValue(c);
  return (dist >= m_minDist && dist <= m_maxDist);
}

////////////////////////////////////////
//Constructors/destructors
RbtDockingSite::RbtDockingSite(const RbtCavityList& cavList, RbtDouble border) : m_cavityList(cavList),m_border(border) {
  if (!m_cavityList.empty()) {
    RbtCoordList minCoords;
    RbtCoordList maxCoords;
    std::transform(m_cavityList.begin(),m_cavityList.end(),std::back_inserter(minCoords),Rbt::ExtractCavityMinCoord);
    std::transform(m_cavityList.begin(),m_cavityList.end(),std::back_inserter(maxCoords),Rbt::ExtractCavityMaxCoord);
    m_minCoord = Rbt::Min(minCoords);
    m_maxCoord = Rbt::Max(maxCoords);
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
  }
}

RbtDockingSite::RbtDockingSite(istream& istr) {
  Read(istr);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtDockingSite::~RbtDockingSite() {
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

//Insertion operator
ostream& operator<<(ostream& s, const RbtDockingSite& site) {
  site.Print(s);
  return s;
}

//Virtual function for dumping docking site details to an output stream
//Derived classes can override if required
void RbtDockingSite::Print(ostream& s) const {
  s << "Total volume " << GetVolume() << " A^3" << endl;
  for(RbtInt i =0; i < m_cavityList.size(); i++) {
    s << "Cavity #" << i+1 << "\t" << *(m_cavityList[i]) << endl;
  }
}

//Public methods
void RbtDockingSite::Write(ostream& ostr) {
  //Make sure grid has been calculated
  if (m_spGrid.Null()) {
    CreateGrid();
  }
  //Write the class name as a title so we can check the authenticity of streams
  //on read
  const char* const header = _CT.c_str();
  RbtInt length = strlen(header);
  Rbt::WriteWithThrow(ostr, (const char*) &length, sizeof(length));
  Rbt::WriteWithThrow(ostr, header, length);
  
  //DM 4 Apr 2002 - write overall min, max coords of all cavities, plus border
  m_minCoord.Write(ostr);
  m_maxCoord.Write(ostr);
  Rbt::WriteWithThrow(ostr, (const char*) &m_border, sizeof(m_border));

  //Write the number of cavities
  RbtInt nCav = m_cavityList.size();
  Rbt::WriteWithThrow(ostr, (const char*) &nCav, sizeof(nCav));

  //Write each cavity
  for (RbtCavityListConstIter cIter = m_cavityList.begin(); cIter != m_cavityList.end(); cIter++) {
    (*cIter)->Write(ostr);
  }

  //DM 4 Apr 2002 - write the distance grid
  if ( (nCav > 0) && !m_spGrid.Null()) {
    m_spGrid->Write(ostr);
  }
}

void RbtDockingSite::Read(istream& istr) {
  m_cavityList.clear();
  m_minCoord=RbtCoord();
  m_maxCoord=RbtCoord();
  m_spGrid = RbtRealGridPtr();
  m_border = 0.0;

  //Read title
  RbtInt length;
  Rbt::ReadWithThrow(istr, (char*) &length, sizeof(length));
  char* header = new char [length+1];
  Rbt::ReadWithThrow(istr, header, length);
  //Add null character to end of string
  header[length] = '\0';
  //Compare title with class name
  RbtBool match = (_CT == header);
  delete [] header;
  if (!match) {
    throw RbtFileParseError(_WHERE_,"Invalid title string in " + _CT + "::Read()");
  }

  //DM 4 Apr 2002 - read overall min, max coords of all cavities, plus border
  m_minCoord.Read(istr);
  m_maxCoord.Read(istr);
  Rbt::ReadWithThrow(istr, (char*) &m_border, sizeof(m_border));

  //Read the number of cavities
  RbtInt nCav;
  Rbt::ReadWithThrow(istr, (char*) &nCav, sizeof(nCav));
  m_cavityList.reserve(nCav);
  //Read each cavity
  for (RbtInt i=0; i < nCav; i++) {
    RbtCavityPtr spCavity = RbtCavityPtr(new RbtCavity(istr));
    m_cavityList.push_back(spCavity);
  }

  if (nCav > 0) {
    //DM 4 Apr 2002 - read the distance grid
    m_spGrid = RbtRealGridPtr(new RbtRealGrid(istr));
  }
  
}

//Return distance grid, calculated on demand if necessary
RbtRealGridPtr RbtDockingSite::GetGrid() {
  if (m_spGrid.Null()) {
    CreateGrid();
  }
  return m_spGrid;
}

//returns total volume of all cavities in A^3
RbtDouble RbtDockingSite::GetVolume() const {
  RbtDouble vol(0.0);
  for (RbtCavityListConstIter iter = m_cavityList.begin(); iter != m_cavityList.end(); iter++) {
    vol += (*iter)->GetVolume();
  }
  return vol;
}

//Returns the combined coord lists of all the cavities
void RbtDockingSite::GetCoordList(RbtCoordList& retVal) const {
  retVal.clear();
  for(RbtCavityListConstIter iter = m_cavityList.begin(); iter != m_cavityList.end(); iter++) {
    const RbtCoordList& cavCoords = (*iter)->GetCoordList();
    retVal.reserve(retVal.size()+cavCoords.size());
    std::copy(cavCoords.begin(),cavCoords.end(),std::back_inserter(retVal));
    //Sort the coords so we can remove any dups
    std::sort(retVal.begin(),retVal.end(),RbtCoordCmp());
    RbtCoordListIter uniqIter = std::unique(retVal.begin(),retVal.end());
    retVal.erase(uniqIter,retVal.end());
    //cout << "Cav = " << cavCoords.size() << "; total = " << retVal.size() << endl;
  }
}


//Filters an atom list according to distance from the cavity coords
//Only returns atoms within minDist and maxDist from cavity
RbtAtomList RbtDockingSite::GetAtomList(const RbtAtomList& atomList, RbtDouble minDist, RbtDouble maxDist) throw (RbtError) {
  if (maxDist > m_border) {
    throw RbtBadArgument(_WHERE_,"maxDist is greater than grid border; recalculate grid");
  }
  RbtAtomList newAtomList;
  if (m_spGrid.Null()) {
    CreateGrid();
  }
  if (m_spGrid.Null()) return newAtomList;
  isAtomInRange bInRange(m_spGrid,minDist,maxDist);
  std::copy_if(atomList.begin(),atomList.end(),std::back_inserter(newAtomList),bInRange);
  return newAtomList;
}

//Filters an atom list according to distance from the cavity coords
//Only returns atoms within maxDist from cavity
//This version does not require the cavity grid
RbtAtomList RbtDockingSite::GetAtomList(const RbtAtomList& atomList, RbtDouble maxDist) throw (RbtError) {
  RbtAtomList newAtomList;
  RbtCoordList allCoords;
  GetCoordList(allCoords);
  Rbt::isAtomNearCoordList bInRange(allCoords,maxDist);
  std::copy_if(atomList.begin(),atomList.end(),std::back_inserter(newAtomList),bInRange);
  return newAtomList;
}

//Returns the count of atoms within minDist and maxDist from cavity
RbtUInt RbtDockingSite::GetNumAtoms(const RbtAtomList& atomList, RbtDouble minDist, RbtDouble maxDist) throw (RbtError) {
  if (maxDist > m_border) {
    throw RbtBadArgument(_WHERE_,"maxDist is greater than grid border; recalculate grid");
  }
  if (m_spGrid.Null()) {
    CreateGrid();
  }
  if (m_spGrid.Null()) return 0;
  isAtomInRange bInRange(m_spGrid,minDist,maxDist);
  return std::count_if(atomList.begin(),atomList.end(),bInRange);
}
  
//Create a grid whose values represent the distance to the nearest cavity coord
void RbtDockingSite::CreateGrid() {
  m_spGrid = RbtRealGridPtr();
  if (m_cavityList.empty()) return;

  //Get the grid step to use from the first cavity
  RbtVector gridStep = m_cavityList.front()->GetGridStep();
  RbtCoord minCoord = m_minCoord-m_border;
  RbtCoord maxCoord = m_maxCoord+m_border;
  RbtVector extent = maxCoord-minCoord;
  RbtUInt nX = int(extent.x/gridStep.x)+1;
  RbtUInt nY = int(extent.y/gridStep.y)+1;
  RbtUInt nZ = int(extent.z/gridStep.z)+1;
  m_spGrid = RbtRealGridPtr(new RbtRealGrid(minCoord,gridStep,nX,nY,nZ));
  m_spGrid->SetAllValues(999999.9);

  //Get the total list of cavity coords
  RbtCoordList allCoords;
  for(RbtCavityListConstIter iter = m_cavityList.begin(); iter != m_cavityList.end(); iter++) {
    const RbtCoordList cavCoords = (*iter)->GetCoordList();
    //Reserve enough space for appending the next cavity coord list
    allCoords.reserve(allCoords.size()+cavCoords.size());
    //Perform the append
    std::copy(cavCoords.begin(),cavCoords.end(),std::back_inserter(allCoords));
    //We can save a bit of time by initialising the distance^2 of the grid points nearest to each cavity coord
    //Normally zero, but we allow the possibility that the grid step may be different
    //i.e. the cavity coords may not lie directly on a grid point.
    for(RbtCoordListConstIter cIter = cavCoords.begin(); cIter != cavCoords.end(); cIter++) {
      RbtUInt i = m_spGrid->GetIXYZ(*cIter);//Grid index of nearest grid point
      RbtDouble dist2 = Rbt::Length2(*cIter,m_spGrid->GetCoord(i));
      m_spGrid->SetValue(i,dist2);
    }
    //Sort the coords so we can remove any dups

    std::sort(allCoords.begin(),allCoords.end(),RbtCoordCmp());
    RbtCoordListIter uniqIter = std::unique(allCoords.begin(),allCoords.end());
    allCoords.erase(uniqIter,allCoords.end());
    cout << "Cav = " << cavCoords.size() << "; total = " << allCoords.size() << endl;
  }

  //Loop over all grid points in the distance grid
  //Can terminate when distance^2 is less than or equal to mindist^2 (shortest length of grid interval)
  RbtDouble mindist2 = std::min(gridStep.x,gridStep.y);
  mindist2 = std::min(mindist2,gridStep.z);
  mindist2 *= mindist2;
  for (RbtUInt i = 0; i < m_spGrid->GetN(); i++) {
    const RbtCoord& c = m_spGrid->GetCoord(i);
    //Initialise dist^2 from initial grid value (999999.9 or 0.0 for cavity coords)
    RbtDouble dist2 = m_spGrid->GetValue(i);
    //Determine min distance to any of the cavity coords
    for (RbtCoordListConstIter iter = allCoords.begin(); iter != allCoords.end() && dist2 > mindist2; iter++) {
      dist2 = std::min(dist2,Rbt::Length2(c-(*iter)));
    }
    m_spGrid->SetValue(i,sqrt(dist2));
  }
}
