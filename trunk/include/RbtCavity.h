/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtCavity.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Struct for holding a coordinate list representing an active site cavity
//Properties are read-only once set in the constructor

#ifndef _RBTCAVITY_H_
#define _RBTCAVITY_H_

#include "RbtConfig.h"
#include "RbtCoord.h"
#include "RbtPrincipalAxes.h"
#include "RbtRealGrid.h"

class RbtCavity
{
 public:
  ////////////////////////////////////////
  //Constructors/destructors

  //Constructor
  RbtCavity(const RbtCoordList& coordList, const RbtVector gridStep) :
    m_coordList(coordList),
    m_gridStep(gridStep),
    m_prAxes(Rbt::GetPrincipalAxes(coordList)),
    m_minCoord(Rbt::Min(coordList)),
    m_maxCoord(Rbt::Max(coordList))
  {
  }
  
  //Constructor reading from binary stream
  RbtCavity(istream& istr)
  {
    Read(istr);
  }

  //Destructor
  ~RbtCavity() {}; 

  ///////////////////////////////////////////////
  //Stream functions
  ///////////////////////////////////////////////

  //Insertion operator
  friend ostream& operator<<(ostream& s, const RbtCavity& cavity) {
    cavity.Print(s);
    return s;
  };

  //Virtual function for dumping cavity details to an output stream
  //Derived classes can override if required
  virtual void Print(ostream& s) const {
    s << "Size=" << m_coordList.size() << " points; Vol=" << GetVolume() << " A^3; Min=" << m_minCoord << "; Max=" << m_maxCoord
      << "; Center=" << m_prAxes.com << "; Extent=" << m_maxCoord-m_minCoord;
  }
  
  //Binary output (serialisation)
  virtual void Write(ostream& ostr) const {
    //DM 4 Apr 2002 - Write grid step
    m_gridStep.Write(ostr);
    //Write number of coords
    RbtInt nCoords = m_coordList.size();
    Rbt::WriteWithThrow(ostr, (const char*) &nCoords, sizeof(nCoords));
    for (RbtCoordListConstIter cIter = m_coordList.begin(); cIter != m_coordList.end(); cIter++) {
      (*cIter).Write(ostr);//Write each coord
    }
  }
  
  //Binary input, replaces existing cavity
  virtual void Read(istream& istr) {
    //Clear the existing cavity
    m_coordList.clear();
    //DM 4 Apr 2002 - Read grid step
    m_gridStep.Read(istr);
    //Read number of coords
    RbtInt nCoords;
    Rbt::ReadWithThrow(istr, (char*) &nCoords, sizeof(nCoords));
    m_coordList.reserve(nCoords);
    RbtCoord c;
    //Read each coord and add it to the cavity list
    for (RbtInt i = 0; i < nCoords; i++) {
      c.Read(istr);
      m_coordList.push_back(c);
    }
    //Recalculate the other properties
    m_prAxes = Rbt::GetPrincipalAxes(m_coordList);
    m_minCoord = Rbt::Min(m_coordList);
    m_maxCoord = Rbt::Max(m_coordList);	
  }

  //DM 4 Apr 2002 - return a grid with all cavity points set to 1.0
  RbtRealGridPtr GetGrid() const {
    RbtVector extent = m_maxCoord-m_minCoord;
    RbtUInt nX = int(extent.x/m_gridStep.x)+3;//+3 to allow a small border
    RbtUInt nY = int(extent.y/m_gridStep.y)+3;
    RbtUInt nZ = int(extent.z/m_gridStep.z)+3;
    RbtRealGridPtr spGrid = RbtRealGridPtr(new RbtRealGrid(m_minCoord-m_gridStep,m_gridStep,nX,nY,nZ));  
    for (RbtCoordListConstIter iter = m_coordList.begin(); iter != m_coordList.end(); iter++) {
      spGrid->SetValue(*iter,1.0);
    }
    return spGrid;
  }

  ////////////////////////////////////////
  //Public methods
  ////////////////
  const RbtCoord& GetCenterOfMass() const {return m_prAxes.com;}
  const RbtPrincipalAxes& GetPrincipalAxes() const {return m_prAxes;}
  RbtInt GetNumCoords() const {return m_coordList.size();}
  const RbtCoordList& GetCoordList() const {return m_coordList;}
  const RbtCoord& GetMinCoord() const {return m_minCoord;}
  const RbtCoord& GetMaxCoord() const {return m_maxCoord;}
  const RbtVector& GetGridStep() const {return m_gridStep;}
  RbtDouble GetVolume() const {return m_coordList.size() * m_gridStep.x * m_gridStep.y * m_gridStep.z;}

 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////


 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtCavity(); //Disable default constructor
  RbtCavity(const RbtCavity&);//Copy constructor disabled by default
  RbtCavity& operator=(const RbtCavity&);//Copy assignment disabled by default


 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtCoordList m_coordList;
  RbtVector m_gridStep;
  RbtPrincipalAxes m_prAxes;
  RbtCoord m_minCoord;
  RbtCoord m_maxCoord;
};

//Useful typedefs
typedef SmartPtr<RbtCavity> RbtCavityPtr;//Smart pointer
typedef vector<RbtCavityPtr> RbtCavityList;//Vector of smart pointers
typedef RbtCavityList::iterator RbtCavityListIter;
typedef RbtCavityList::const_iterator RbtCavityListConstIter;

namespace Rbt
{
  ////////////////////////////////////////////////////////
  //Comparison functions for sorting RbtCavityPtr containers
  //For use by STL sort algorithms
  ////////////////////////////////////////////////////////
  //Less than operator for sorting RbtCavityPtrs into ascending order by center of mass distance from a given coord
  class RbtCavityPtrCmp_Distance {
    const RbtCoord& c;
  public:
    RbtCavityPtrCmp_Distance(const RbtCoord& cc) : c(cc) {}
    RbtBool operator()(RbtCavityPtr spCav1, RbtCavityPtr spCav2) const {
      return Rbt::Length2(spCav1->GetCenterOfMass()-c) < Rbt::Length2(spCav2->GetCenterOfMass()-c);
    }
  };
  //Less than operator for sorting RbtCavityPtrs into descending order by volume
  class RbtCavityPtrCmp_Volume {
  public:
    RbtCavityPtrCmp_Volume() {}
    RbtBool operator()(RbtCavityPtr spCav1, RbtCavityPtr spCav2) const {
      return spCav1->GetNumCoords() > spCav2->GetNumCoords();
    }
  };
  ////////////////////////////////////////////////////////
  //Predicate functions for RbtCavityPtr
  //For use by STL algorithms
  ////////////////////////////////////////////////////////
  class isCavityNearCoord : public std::unary_function<RbtCavityPtr,RbtBool> {
    const RbtCoord& c;
    RbtDouble r2;//radius squared (to avoid taking square roots)
  public:
    explicit isCavityNearCoord(const RbtCoord& cc, RbtDouble rr) : c(cc), r2(rr*rr) {}
    RbtBool operator() (RbtCavityPtr spCavity) const {return Rbt::Length2(spCavity->GetCenterOfMass()-c) <= r2;}
  };

  ////////////////////////////////////////////
  //Functions objects for performing actions on cavities
  //For use by STL algorithms
  ////////////////////////////////////////////
  inline const RbtCoord& ExtractCavityMinCoord(RbtCavityPtr spCav) {return spCav->GetMinCoord();}
  inline const RbtCoord& ExtractCavityMaxCoord(RbtCavityPtr spCav) {return spCav->GetMaxCoord();}
}

#endif //_RBTCAVITY_H_
