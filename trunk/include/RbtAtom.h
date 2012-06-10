/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtAtom.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Handles all the properties for a single atom

#ifndef _RBTATOM_H_
#define _RBTATOM_H_

#include <functional> //For STL function objects
using std::ptr_fun;
using std::bind2nd;

#include "RbtConfig.h"
#include "RbtCoord.h"
#include "RbtQuat.h"
#include "RbtPMF.h"
#include "RbtTriposAtomType.h"

class RbtModel; //Forward declaration
class RbtBond; //Forward declaration

namespace Rbt {
	//DM 23 Apr 1999 - provide custom comparison function
	//to bond map, so that map is sorted by bond ID, not by pointer
	class RbtBondPCmp_BondId {
	public:
		RbtBool operator()(RbtBond* pBond1, RbtBond* pBond2) const;
	};
}

//Map of bonds the atom is bonded to
//Key is regular RbtBond*, value is RbtBool
//RbtBool is true if this atom is atom1 in the bond, false if atom2
typedef map<RbtBond*,RbtBool,Rbt::RbtBondPCmp_BondId> RbtBondMap;
typedef RbtBondMap::iterator RbtBondMapIter;
typedef RbtBondMap::const_iterator RbtBondMapConstIter;


class RbtAtom
{
 public:

  ///////////////////////////////////////////////
  //Enums
  ///////////////////////////////////////////////

  //Hybridisation state
  enum eHybridState
  {
    UNDEFINED = 0,
    SP   = 1, //SP hybridised (linear geom)
    SP2  = 2, //SP2 hybridised (planar geom)
    SP3  = 3, //SP3 hybridised (tetrahedral geom)
    AROM = 4, //Aromatic SP2 (planar geom)
    TRI  = 5, //Trigonal nitrogen (planar geom)    
  };

  ///////////////////////////////////////////////
  //Constructors / destructors
  ///////////////////////////////////////////////

  //Default constructor
  RbtAtom();

  //Constructor supplying all 2-D parameters
  RbtAtom(RbtInt nAtomId,
	  RbtInt 		nAtomicNo		= 6,
	  RbtString 	strAtomName		= "C",
	  RbtString 	strSubunitId	= "1",
	  RbtString 	strSubunitName	= "RES",
	  RbtString 	strSegmentName	= "SEG1",
      eHybridState	eState			= UNDEFINED,//DM 8 Dec 1998 Changed from SP3 to UNDEFINED
	  RbtUInt		nHydrogens		= 0,
	  RbtInt		nFormalCharge	= 0 
	  );

  //Default destructor
  virtual ~RbtAtom();

  //Copy constructor
  RbtAtom(const RbtAtom& atom);

  //Copy assignment
  RbtAtom& operator=(const RbtAtom& atom);



  ///////////////////////////////////////////////
  //Stream functions
  ///////////////////////////////////////////////

  //Insertion operator (primarily for debugging)
  //Note: needs to be friend so can access atom private data
  //without using the accessor functions
  friend ostream& operator<<(ostream& s, const RbtAtom& atom);

  //Virtual function for dumping atom details to an output stream
  //Derived classes (e.g. pseudoatom) can override if required
  virtual ostream& Print(ostream& s) const;



  ///////////////////////////////////////////////
  //Public accessor functions
  //2-D attributes
  ///////////////////////////////////////////////

    //DM 7 June 2006
    //atom enabled state is controlled by the parent RbtModel
    //enabled state
    RbtBool GetEnabled() const;
  
  //AtomId
  RbtInt GetAtomId() const {return m_nAtomId;};
  void SetAtomId(const RbtInt nAtomId) {m_nAtomId = nAtomId;};

  //AtomicNo
  RbtInt GetAtomicNo() const {return m_nAtomicNo;};
  void SetAtomicNo(const RbtInt nAtomicNo) {m_nAtomicNo = nAtomicNo;};

  //AtomName
  RbtString GetAtomName() const {return m_strAtomName;};
  void SetAtomName(const RbtString strAtomName) {m_strAtomName = strAtomName;};
  RbtString GetFullAtomName() const;//Returns composite of segment name, subunit id and name, and atom name

  //SubunitId
  RbtString GetSubunitId() const {return m_strSubunitId;}
  void SetSubunitId(const RbtString strSubunitId) {m_strSubunitId = strSubunitId;}

  //SubunitName
  RbtString GetSubunitName() const {return m_strSubunitName;}
  void SetSubunitName(const RbtString strSubunitName) {m_strSubunitName = strSubunitName;}

  //SegmentName
  RbtString GetSegmentName() const {return m_strSegmentName;}
  void SetSegmentName(const RbtString strSegmentName) {m_strSegmentName = strSegmentName;}

  //Hybridisation state
  eHybridState GetHybridState() const {return m_eState;}
  void SetHybridState(const eHybridState eState) {m_eState = eState;}

  //Attached hydrogens
  RbtUInt GetNumImplicitHydrogens() const {return m_nHydrogens;}
  void SetNumImplicitHydrogens(const RbtUInt nHydrogens) {m_nHydrogens = nHydrogens;}

  //FormalCharge (DM 24 Mar 1999 - changed from double to int)
  RbtInt GetFormalCharge() const {return m_nFormalCharge;}
  void SetFormalCharge(const RbtInt nFormalCharge) {m_nFormalCharge = nFormalCharge;}

  //CyclicFlag - flag to indicate atom is in a ring (set by RbtModel::FindRing)
  RbtBool GetCyclicFlag() const {return m_bCyclic;}
  void SetCyclicFlag(RbtBool bCyclic=true) {m_bCyclic = bCyclic;}

  //SelectionFlag - general purpose flag can be set/cleared by various search algorithms (e.g. FindRings)
  RbtBool GetSelectionFlag() const {return m_bSelected;}
  void SetSelectionFlag(RbtBool bSelected=true) {m_bSelected = bSelected;}
  void InvertSelectionFlag() {m_bSelected = !m_bSelected;}

  //DM 21 Jul 1999 - user-defined flag
  RbtBool GetUser1Flag() const {return m_bUser1;}
  void SetUser1Flag(RbtBool bUser1=true) {m_bUser1 = bUser1;}
  
  //DM 29 Jan 2000 - user-defined double
  RbtDouble GetUser1Value() const {return m_dUser1;}
  void SetUser1Value(RbtDouble dUser1=0.0) {m_dUser1 = dUser1;}

  //DM 27 Jul 2000 - user-defined double#2
  RbtDouble GetUser2Value() const {return m_dUser2;}
  void SetUser2Value(RbtDouble dUser2=0.0) {m_dUser2 = dUser2;}

  //ModelPtr - pointer to the parent model
  RbtModel* GetModelPtr() const {return m_pModel;}
  void SetModelPtr(RbtModel* pModel = NULL) {m_pModel = pModel;}

  //DM 04 Dec 1998  Add functions to handle bond map
  //Returns number of bonds in map
  RbtUInt GetNumBonds() const {return m_bondMap.size();}
  //Returns the bond map
  //DM 26 Oct 2000 - return by reference
  //RbtBondMap GetBondMap() const {return m_bondMap;}
  const RbtBondMap& GetBondMap() const {return m_bondMap;}
  //Add a bond to the bond map - returns true if OK, false if this atom is not a member of the bond
  //or if bond already added
  RbtBool AddBond(RbtBond* pBond);
  //Remove a bond from the bond map - returns true if OK, false if bond not present in map
  RbtBool RemoveBond(RbtBond* pBond);
  //Returns number of cyclic bonds in map
  RbtUInt GetNumCyclicBonds() const;
  //Returns bond map with only cyclic bonds included
  RbtBondMap GetCyclicBondMap() const;
  //Returns total formal bond order around atom
  RbtInt GetTotalFormalBondOrder() const;
  //Returns max formal bond order for all bonds to atom
  RbtInt GetMaxFormalBondOrder() const;



  ///////////////////////////////////////////////
  //Public accessor functions
  //3-D attributes
  ///////////////////////////////////////////////

  //Coords
  //DM 8 Dec 1998 - make virtual so we can override in the pseudoatom class
  //DM 11 Jul 2000 - remove overhead of virtual methods
  //Pseudo atoms must now refresh their coords each time the constituent atoms move
  //DM 27 Oct 2000 - return by reference
  const RbtCoord& GetCoords() const {return m_coord;} 
  RbtDouble GetX() const {return m_coord.x;}
  RbtDouble GetY() const {return m_coord.y;}
  RbtDouble GetZ() const {return m_coord.z;}

  void SetCoords(const RbtCoord& coord) {m_coord = coord;}
  void SetCoords(const RbtDouble x, const RbtDouble y, const RbtDouble z) {m_coord.x = x; m_coord.y = y; m_coord.z = z;}
  void SetX(const RbtDouble x) {m_coord.x = x;}
  void SetY(const RbtDouble y) {m_coord.y = y;}
  void SetZ(const RbtDouble z) {m_coord.z = z;}

  //PartialCharge
  RbtDouble GetPartialCharge() const {return m_dPartialCharge;}
  void SetPartialCharge(const RbtDouble dPartialCharge) {m_dPartialCharge = dPartialCharge;}

  //GroupCharge (added DM 24 Mar 1999, for ionic interaction group charges)
  RbtDouble GetGroupCharge() const {return m_dGroupCharge;}
  void SetGroupCharge(const RbtDouble dGroupCharge) {m_dGroupCharge = dGroupCharge;}

  //AtomicMass
  RbtDouble GetAtomicMass() const {return m_dAtomicMass;}
  void SetAtomicMass(const RbtDouble dAtomicMass) {m_dAtomicMass = dAtomicMass;}

  //VdwRadius
  RbtDouble GetVdwRadius() const {return m_dVdwRadius;}
  void SetVdwRadius(const RbtDouble dVdwRadius) {m_dVdwRadius = dVdwRadius;}

  //AtomType
  RbtString GetFFType() const {return m_strFFType;}
  void SetFFType(const RbtString& strFFType) {m_strFFType = strFFType;}
  RbtPMFType GetPMFType() const {return m_nPMFType;} 
  void SetPMFType(RbtPMFType aType) {m_nPMFType = aType;}
  RbtTriposAtomType::eType GetTriposType() const {return m_triposType;}
  void SetTriposType(RbtTriposAtomType::eType aType) {m_triposType = aType;}



  ///////////////////////////////////////////////
  //Other public methods
  ///////////////////////////////////////////////
  //DM 11 Jan 1999 - give atoms the ability to save their current coords
  //and revert to the saved coords. Useful for Monte Carlo simulations.
  //DM 08 Feb 1999 - all saved coords are now saved in a map<RbtUInt,RbtCoord>
  //map key=0 is reserved for the default SaveCoords and RevertCoords
  void SaveCoords(RbtUInt coordNum = 0);
  void RevertCoords(RbtUInt coordNum = 0) throw (RbtError);

  //Translate - translate coordinates by the supplied vector
  void Translate(const RbtVector& vector) {m_coord += vector;}

  void Translate(const RbtDouble vx, const RbtDouble vy, const RbtDouble vz) {m_coord += RbtCoord(vx,vy,vz);}


  //DM 07 Jan 1999 - rotate coordinates using the supplied quaternion
  void RotateUsingQuat(const RbtQuat& q) {m_coord = q.Rotate(m_coord);}


  //DM 04 Dec 1998  Now we have the bond map, we can easily provide coordination numbers
  //This version returns the total number of coordinated atoms (includes implicit hydrogens)
  //Equivalent to GetNumBonds()+GetNumImplicitHydrogens()
  RbtUInt GetCoordinationNumber() const;
  //This version returns the number of coordinated atoms of a given element
  //Note: for hydrogens, the implicit atoms are excluded
  //so GetCoordinationNumber(1) will return the number of explicit hydrogens
  RbtUInt GetCoordinationNumber(RbtInt nAtomicNo) const;
  //This version returns the number of coordinated atoms of a given force field type
  RbtUInt GetCoordinationNumber(const RbtString& strFFType) const;
  //This version returns the number of coordinated atoms of a given hybridisation state
  RbtUInt GetCoordinationNumber(eHybridState e) const;


 protected:

 private:
  //Private methods
  //Clears the bond map - should only need to be called by the copy constructors, hence private
  void ClearBondMap();

 private:
  //Private data

  //These can be considered as 2-D params (i.e. define the chemistry and topology)
  RbtInt 		m_nAtomicNo;		// Atomic number
  RbtInt 		m_nAtomId;			// Atom ID
  RbtString 	m_strAtomName;		// Atom name
  RbtString 	m_strSubunitId;		// Subunit(residue) ID (note: string not int)
  RbtString 	m_strSubunitName;	// Subunit(residue) name
  RbtString 	m_strSegmentName;	// Segment(molecule) name
  eHybridState	m_eState;			// Hybridisation state
  RbtUInt		m_nHydrogens;		// Number of attached (implicit) hydrogens
  RbtInt		m_nFormalCharge;	// Formal charge (DM 24 Mar 1999 - changed from double to int)
  RbtModel*		m_pModel;			// Regular pointer to parent model
  RbtBondMap	m_bondMap;			// Map of bonds this atom is bonded to
  RbtBool		m_bCyclic;			// Is the atom in a ring ?
  RbtBool		m_bSelected;		// Can be set/cleared by various search algorithms (e.g. FindRings)
  RbtBool		m_bUser1;			// User flag 1
  RbtDouble		m_dUser1;			// User value 1
  RbtDouble		m_dUser2;			// User value 2
  RbtPMFType	m_nPMFType;			// PMF atom type: default
  RbtTriposAtomType::eType m_triposType;//DM 14 May 2002: Tripos 5.2 atom type

  //These can be considered as 3-D params (i.e. the extra info required for 3-D calculations)
  RbtCoord m_coord; //cartesian coords
  RbtDouble m_dPartialCharge; //partial charge
  RbtDouble m_dGroupCharge; //interaction group charge (added DM 24 Mar 1999)
  RbtDouble m_dAtomicMass; //atomic mass
  RbtDouble m_dVdwRadius; //atomic mass
  RbtString m_strFFType; //force field atom type

  RbtUIntCoordMap m_savedCoords; //DM 08 Feb 1999 - now store all saved coords in a map<RbtUInt,RbtCoord>
};


//Useful typedefs
typedef SmartPtr<RbtAtom> RbtAtomPtr;//Smart pointer
typedef vector<RbtAtomPtr> RbtAtomList;//Vector of smart pointers
typedef RbtAtomList::iterator RbtAtomListIter;
typedef RbtAtomList::const_iterator RbtAtomListConstIter;

typedef vector<RbtAtom*> RbtAtomRList;//Vector of regular pointers
typedef RbtAtomRList::iterator RbtAtomRListIter;
typedef RbtAtomRList::const_iterator RbtAtomRListConstIter;

typedef vector<RbtAtomList> RbtAtomListList;//A vector of atom vectors (e.g. for storing ring systems)
typedef RbtAtomListList::iterator RbtAtomListListIter;
typedef RbtAtomListList::const_iterator RbtAtomListListConstIter;

typedef vector<RbtAtomRList> RbtAtomRListList;//A vector of atom vectors (regular pointer version)
typedef RbtAtomRListList::iterator RbtAtomRListListIter;
typedef RbtAtomRListList::const_iterator RbtAtomRListListConstIter;

typedef list<RbtAtomPtr> RbtAtomTrueList; //A real list of smart pointers
typedef RbtAtomTrueList::iterator RbtAtomTrueListIter;
typedef RbtAtomTrueList::const_iterator RbtAtomTrueListConstIter;


///////////////////////////////////////////////
// Non-member functions (in Rbt namespace)
//////////////////////////////////////////
namespace Rbt
{
  ////////////////////////////////////////////
  //Calculation functions
  ///////////////////////

  //DM 07 Jul 2000 - converted to regular RbtAtom* versions
  //Returns the "bond" length in Angstroms between 2 atoms
  //Note: atoms do not have to be bonded
  inline RbtDouble BondLength(RbtAtom* pAtom1, RbtAtom* pAtom2) {
    return Rbt::Length(pAtom1->GetCoords(),pAtom2->GetCoords());
  }

  //Returns the "bond" angle in degrees between 3 atoms
  //Note: atoms do not have to be bonded
  inline RbtDouble BondAngle(RbtAtom* pAtom1, RbtAtom* pAtom2, RbtAtom* pAtom3) {
    return Rbt::Angle(pAtom1->GetCoords(),pAtom2->GetCoords(),pAtom3->GetCoords());
  }

  //Returns the "bond" dihedral in degrees between 4 atoms
  //Note: atoms do not have to be bonded
  inline RbtDouble BondDihedral(RbtAtom* pAtom1, RbtAtom* pAtom2, RbtAtom* pAtom3, RbtAtom* pAtom4) {
    return Rbt::Dihedral(pAtom1->GetCoords(),pAtom2->GetCoords(),pAtom3->GetCoords(),pAtom4->GetCoords());
  }


  ////////////////////////////////////////////
  //Conversion functions
  ///////////////////////
  
  //Converts hybridisation state enum to a string
  RbtString ConvertHybridStateToString(RbtAtom::eHybridState eState);
  //Convert formal charge to a string (e.g. +, -, +2, -2 etc)
  RbtString ConvertFormalChargeToString(RbtInt nCharge);


  ////////////////////////////////////////////
  //Predicates (for use by STL algorithms)
  ////////////////////////////////////////////
  //DM 30 Apr 2002 - convert to regular RbtAtom* functions
  //More universal as RbtAtomPtr parameters will be automatically degraded to
  //regular pointers

  //Is atom enabled ?
  class isAtomEnabled : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomEnabled() {}
    RbtBool operator() (RbtAtom* pAtom) const {return pAtom->GetEnabled();}
  };

  //Is atom selected ?
  class isAtomSelected : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomSelected() {}
    RbtBool operator() (RbtAtom* pAtom) const {return pAtom->GetSelectionFlag();}
  };

  //Is atom cyclic ?
  class isAtomCyclic : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomCyclic() {}
    RbtBool operator() (const RbtAtom* pAtom) const {return pAtom->GetCyclicFlag();}
  };

  //DM 9 Feb 2000 - is atom bridgehead? Checks for >2 cyclic bonds
  class isAtomBridgehead : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomBridgehead() {}
    RbtBool operator() (const RbtAtom* pAtom) const {return (pAtom->GetNumCyclicBonds() > 2);}
  };

  //Is atom a H-Bond Acceptor ?
  //Checks atomic number for O(8), N(7) or S(16)
  class isAtomHBondAcceptor : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomHBondAcceptor() {}
    RbtBool operator() (const RbtAtom*) const;
  };

  //Is atom a H-Bond Donor ?
  //Checks 1) is it a hydrogen, 2) does it make exactly one bond, 3) is the bond to O, N or S ?
  class isAtomHBondDonor : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomHBondDonor() {}
    RbtBool operator() (const RbtAtom*) const;
  };

  //Is atom formally charged ?
  //Checks if formal charge is != zero
  class isAtomCharged : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomCharged() {}
    RbtBool operator() (const RbtAtom* pAtom) const {return (pAtom->GetFormalCharge() != 0);}
  };

  //Is atom formally positively charged ?
  //Checks if formal charge is > zero
  class isAtomPosCharged : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomPosCharged() {}
    RbtBool operator() (const RbtAtom* pAtom) const {return (pAtom->GetFormalCharge() > 0);}
  };

  //Is atom formally negatively charged ?
  //Checks if formal charge is < zero
  class isAtomNegCharged : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomNegCharged() {}
    RbtBool operator() (const RbtAtom* pAtom) const {return (pAtom->GetFormalCharge() < 0);}
  };

  //Does atom have implicit hydrogens
  class isAtomExtended : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomExtended() {}
    RbtBool operator() (const RbtAtom* pAtom) const {return (pAtom->GetNumImplicitHydrogens() > 0);}
  };

  //Is atom planar ?
  //Checks if 1) atom makes 2 bonds (in which case must be planar) or
  //2) hybridisation state is SP2, AROM or TRI
  class isAtomPlanar : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomPlanar() {}
    RbtBool operator() (const RbtAtom*) const;
  };

  //Is atom a pi-atom ?
  //SP2,TRI or AROM, or special case of OSP3
  class isPiAtom : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isPiAtom() {}
    RbtBool operator() (const RbtAtom*) const;
  };

  //Is atomic number equal to n ?
  class isAtomicNo_eq : public std::unary_function<RbtAtom*,RbtBool> {
    RbtInt n;
  public:
    explicit isAtomicNo_eq(RbtInt nn) : n(nn) {}
    RbtBool operator() (const RbtAtom* pAtom) const {return pAtom->GetAtomicNo() == n;}
  };

  //Is Force field type equal to s ?
  class isFFType_eq : public std::unary_function<RbtAtom*,RbtBool> {
    RbtString s;
  public:
    explicit isFFType_eq(const RbtString& ss) : s(ss) {}
    RbtBool operator() (const RbtAtom* pAtom) const {return pAtom->GetFFType() == s;}
  };

  //Is Atom name equal to s ?
  class isAtomName_eq : public std::unary_function<RbtAtom*,RbtBool> {
    RbtString s;
  public:
    explicit isAtomName_eq(const RbtString& ss) : s(ss) {}
    RbtBool operator() (const RbtAtom* pAtom) const {return pAtom->GetAtomName() == s;}
  };

  //Is Subunit name equal to s ?
  class isSubunitName_eq : public std::unary_function<RbtAtom*,RbtBool> {
    RbtString s;
  public:
    explicit isSubunitName_eq(const RbtString& ss) : s(ss) {}
    RbtBool operator() (const RbtAtom* pAtom) const {return pAtom->GetSubunitName() == s;}
  };

  //Is Subunit ID equal to s ?
  class isSubunitId_eq : public std::unary_function<RbtAtom*,RbtBool> {
    RbtString s;
  public:
    explicit isSubunitId_eq(const RbtString& ss) : s(ss) {}
    RbtBool operator() (const RbtAtom* pAtom) const {return pAtom->GetSubunitId() == s;}
  };

  //Is Segment name equal to s ?
  class isSegmentName_eq : public std::unary_function<RbtAtom*,RbtBool> {
    RbtString s;
  public:
    explicit isSegmentName_eq(const RbtString& ss) : s(ss) {}
    RbtBool operator() (const RbtAtom* pAtom) const {return pAtom->GetSegmentName() == s;}
  };

  //Is hybridisation state equal to e ? (DM 26 Mar 1999)
  class isHybridState_eq : public std::unary_function<RbtAtom*,RbtBool> {
    RbtAtom::eHybridState e;
  public:
    explicit isHybridState_eq(RbtAtom::eHybridState ee) : e(ee) {}
    RbtBool operator() (const RbtAtom* pAtom) const {return pAtom->GetHybridState() == e;}
  };

  //DM 6 Jan 1999
  //Is atom2 equal to atom1 (checks if underlying regular pointers match)
  //Note: this is a binary rather than unary predicate
  //DM 1 Feb 1999 - renamed from isAtom_eq to isAtomPtr_eq
  class isAtomPtr_eq : public std::binary_function<RbtAtom*,RbtAtom*,RbtBool> {
  public:
    explicit isAtomPtr_eq() {}
    RbtBool operator() (const RbtAtom* pAtom1, const RbtAtom* pAtom2) const {return pAtom1 == pAtom2;}
  };

  //DM 1 Feb 1999
  //Is atom2 equal to atom1 (checks if subunit name, subunit ID and atom name match)
  //Note: this is a binary rather than unary predicate
  class isAtom_eq : public std::binary_function<RbtAtom*,RbtAtom*,RbtBool> {
  public:
    explicit isAtom_eq() {}
    RbtBool operator() (const RbtAtom* pAtom1, const RbtAtom* pAtom2) const {
      return ( (pAtom1->GetSubunitId() == pAtom2->GetSubunitId()) &&
	       (pAtom1->GetSubunitName() == pAtom2->GetSubunitName()) &&
	       (pAtom1->GetAtomName() == pAtom2->GetAtomName()) );
    }
  };

  //DM 22 Apr 2002
  //Is atom2 ID equal to atom1 ID
  //Note: this is a binary rather than unary predicate
  class isAtomId_eq : public std::binary_function<RbtAtom*,RbtAtom*,RbtBool> {
  public:
    explicit isAtomId_eq() {}
    RbtBool operator() (const RbtAtom* pAtom1, const RbtAtom* pAtom2) const {
      return (pAtom1->GetAtomId() == pAtom2->GetAtomId());
    }
  };

  //DM 12 Apr 1999
  //Is atom inside a sphere defined by radius=r, center=c
  class isAtomInsideSphere : public std::unary_function<RbtAtom*,RbtBool> {
    const RbtCoord& c;//center of sphere
    RbtDouble r2;//radius squared (to avoid taking square roots)
  public:
    explicit isAtomInsideSphere(const RbtCoord& cc, RbtDouble rr) : c(cc), r2(rr*rr) {}
    RbtBool operator() (const RbtAtom* pAtom) const {return Rbt::Length2(pAtom->GetCoords()-c) <= r2;}
  };

  //DM 12 Apr 1999
  //Is atom within a cuboid defined by cmin,cmax coords?
  class isAtomInsideCuboid : public std::unary_function<RbtAtom*,RbtBool> {
    const RbtCoord& cmin;
    const RbtCoord& cmax;
  public:
    explicit isAtomInsideCuboid(const RbtCoord& ccmin, const RbtCoord& ccmax) : cmin(ccmin), cmax(ccmax) {}
    RbtBool operator() (const RbtAtom* pAtom) const {
    	//DM 27 Oct 2000 - GetCoords now returns by reference
      const RbtCoord& c = pAtom->GetCoords();
      return (c >= cmin) && (c <= cmax);
    }
  };


  //DM 29 Jul 1999
  //Is atom within a given distance of any coord in the coord list
  class isAtomNearCoordList : public std::unary_function<RbtAtom*,RbtBool> {
    const RbtCoordList& cl;
    RbtDouble r2;//radius squared (to avoid taking square roots)
  public:
    explicit isAtomNearCoordList(const RbtCoordList& ccl, RbtDouble rr) : cl(ccl), r2(rr*rr) {}
    RbtBool operator() (const RbtAtom* pAtom) const;
  };


  //Is atom2 1-2 connected (i.e. bonded) to atom1 ?
  class isAtom_12Connected : public std::unary_function<RbtAtom*,RbtBool> {
    RbtAtom* pAtom1;
    RbtAtomList bondedAtomList1;
  public:
    explicit isAtom_12Connected(RbtAtom* spAtom);
    RbtBool operator() (RbtAtom* spAtom2) const;
  };

  //Is atom2 1-3 connected (i.e. via a bond angle) to atom1 ?
  class isAtom_13Connected : public std::unary_function<RbtAtom*,RbtBool> {
    RbtAtom* pAtom1;
    RbtAtomList bondedAtomList1;
  public:
    explicit isAtom_13Connected(RbtAtom* pAtom);
    RbtBool operator() (RbtAtom* pAtom2) const;
  };


  //Is atom's subunit an RNA-type (i.e. A,G,C or U)
  class isAtomRNA : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomRNA() {}
    RbtBool operator() (const RbtAtom*) const;
  };

  //Is atom classified as an RNA receptor "anionic" atom
  //For our purposes, these are N7,O2,O4,O6,O1P and O2P
  //DM 22 Dec 1998 - center phosphate interactions on the P atom, not on O1P and O2P
  //to avoid double counting of favourable interactions
  //class isAtomRNAAnionic : public std::unary_function<RbtAtomPtr,RbtBool> {
  //  isAtomRNA bIsAtomRNA;
  //public:
  //  explicit isAtomRNAAnionic() {}
  //  RbtBool operator() (RbtAtomPtr) const;
  //};
  
  //DM 22 Dec 1998 - need ligand anionic atoms so we can include repulsive
  //interactions with the RNA
  //
  //Is atom classified as a ligand "anionic" atom
  //Includes COO-, PO3-
  //class isAtomLigandAnionic : public std::unary_function<RbtAtomPtr,RbtBool> {
  //public:
  //  explicit isAtomLigandAnionic() {}
  //  RbtBool operator() (RbtAtomPtr) const;
  //};

  //DM 22 Dec 1998 - renamed from isAtomRNACationic
  //
  //Is atom classified as a ligand "cationic" atom
  //The most important category is charged Nitrogens
  //class isAtomLigandCationic : public std::unary_function<RbtAtomPtr,RbtBool> {
  //public:
  //  explicit isAtomLigandCationic() {}
  //  RbtBool operator() (RbtAtomPtr) const;
  //};

  //////////////////////////////
  //DM 30 Mar 1999
  //Now that the MolecularFileSources define the GroupCharge attribute in a file-format independent manner,
  //we don't need separate predicates for isRNAAnionic, isLigandAnionic etc
  //We just need isIonic, isCationic, isAnionic

  //Is atom defined as an ionic interaction center ?
  //Checks if group charge is != zero
  class isAtomIonic : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomIonic() {}
    RbtBool operator() (const RbtAtom* pAtom) const {return (fabs(pAtom->GetGroupCharge()) > 0.001);}
  };

  //Is atom defined as a cationic interaction center ?
  //Checks if group charge is > zero
  class isAtomCationic : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomCationic() {}
    RbtBool operator() (const RbtAtom* pAtom) const {return (pAtom->GetGroupCharge() > 0.001);}
  };

  //Is atom defined as an anionic interaction center ?
  //Checks if group charge is < zero
  class isAtomAnionic : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomAnionic() {}
    RbtBool operator() (const RbtAtom* pAtom) const {return (pAtom->GetGroupCharge() < -0.001);}
  };

  //Is atom the central carbon in a guanidinium group ?
  //Checks for cationic sp2/arom carbon
  //DM 27 Jan 2000 - also check is acyclic
  class isAtomGuanidiniumCarbon : public std::unary_function<RbtAtom*,RbtBool> {
    Rbt::isAtomicNo_eq bIsCarbon;
    Rbt::isAtomCationic bIsCationic;
    Rbt::isPiAtom bIsPi;
    Rbt::isAtomCyclic bIsCyclic;
  public:
    explicit isAtomGuanidiniumCarbon() : bIsCarbon(6) {}
    RbtBool operator() (const RbtAtom* pAtom) const {return bIsCationic(pAtom) && bIsCarbon(pAtom) && bIsPi(pAtom)
						    && !bIsCyclic(pAtom);}
  };

  //DM 24 Jan 2001
  //Checks for common metal ions by atomic number
  //DM 19 Oct 2001 - extend to all common metals (Na,Mg,K->Zn)
  class isAtomMetal : public std::unary_function<RbtAtom*,RbtBool> {
  public:
    explicit isAtomMetal();
    RbtBool operator() (const RbtAtom* pAtom) const;
  };

  //DM 21 Jul 1999 Is atom lipophilic ?
  //DM 16 May 2003 Total rewrite to be much more comprehensive
  class isAtomLipophilic : public std::unary_function<RbtAtom*,RbtBool> {
    Rbt::isAtomIonic isIonic;
    Rbt::isAtomHBondDonor isHBD;
    Rbt::isAtomHBondAcceptor isHBA;
    Rbt::isAtomMetal isMetal;
    Rbt::isHybridState_eq isSP3;
    Rbt::isHybridState_eq isSP2;
    Rbt::isAtomicNo_eq isO;
    Rbt::isAtomicNo_eq isN;
  public:
    explicit isAtomLipophilic() : isSP3(RbtAtom::SP3), isSP2(RbtAtom::SP2), isO(8), isN(7) {}
    RbtBool operator() (RbtAtom* pAtom) const;
  };

  //////////////////////////////

  //DM 6 April 1999
  //Is coordination number of atom equal to n
  //First constructor checks total coordination number
  //Other constructors check atomic number, force field type and hybridisation state coordination numbers
  class isCoordinationNumber_eq : public std::unary_function<RbtAtom*,RbtBool> {
    enum {TOTAL,ATNO,FFTYPE,HYBRID} eCNType;//Type of coordination number to check
    RbtInt n;//Coordination number value to check
    RbtInt atNo;
    RbtString ffType;
    RbtAtom::eHybridState hybrid;
  public:
    //Total coordination number
    explicit isCoordinationNumber_eq(RbtInt nn) :
      n(nn),eCNType(TOTAL),atNo(0),ffType(""),hybrid(RbtAtom::UNDEFINED) {}
    //Atomic number coordination number
    explicit isCoordinationNumber_eq(RbtInt nn, RbtInt nAt) :
      n(nn),eCNType(ATNO),atNo(nAt),ffType(""),hybrid(RbtAtom::UNDEFINED) {}
    //Force field type coordination number
    explicit isCoordinationNumber_eq(RbtInt nn, const RbtString& strType) :
      n(nn),eCNType(FFTYPE),atNo(0),ffType(strType),hybrid(RbtAtom::UNDEFINED) {}
    //Hybridisation state coordination number
    explicit isCoordinationNumber_eq(RbtInt nn, RbtAtom::eHybridState eState) :
      n(nn),eCNType(HYBRID),atNo(0),ffType(""),hybrid(eState) {}
    RbtBool operator() (RbtAtom* pAtom) const {
      switch (eCNType) {
      case TOTAL:
	return pAtom->GetCoordinationNumber() == n;
      case ATNO:
	return pAtom->GetCoordinationNumber(atNo) == n;
      case FFTYPE:
	return pAtom->GetCoordinationNumber(ffType) == n;
      case HYBRID:
	return pAtom->GetCoordinationNumber(hybrid) == n;      
      default:
	return false;
      }
    }
  };

  ////////////////////////////////////////////////////////
  //Comparison functions for sorting RbtAtom* containers
  //For use by STL sort algorithms
  ////////////////////////////////////////////////////////

  //Less than operator for sorting RbtAtom*s by atom ID
  class RbtAtomPtrCmp_AtomId {
  public:
    RbtBool operator()(RbtAtom* pAtom1, RbtAtom* pAtom2) const {
      return pAtom1->GetAtomId() < pAtom2->GetAtomId();
    }
  };

  //Less than operator for sorting RbtAtom*s by atomic number
  class RbtAtomPtrCmp_AtomicNo {
  public:
    RbtBool operator()(RbtAtom* pAtom1, RbtAtom* pAtom2) const {
      return pAtom1->GetAtomicNo() < pAtom2->GetAtomicNo();
    }
  };

  //Less than operator for sorting RbtAtom*s by atom name
  class RbtAtomPtrCmp_AtomName {
  public:
    RbtBool operator()(RbtAtom* pAtom1, RbtAtom* pAtom2) const {
      return pAtom1->GetAtomName() < pAtom2->GetAtomName();
    }
  };

  //Less than operator for sorting RbtAtomPtrs by pointer value
  class RbtAtomPtrCmp_Ptr {
  public:
    RbtBool operator()(RbtAtom* pAtom1, RbtAtom* pAtom2) const {
      return pAtom1 < pAtom2;
    }
  };

  ////////////////////////////////////////////
  //Functions objects for performing actions on atoms
  //For use by STL algorithms
  ////////////////////////////////////////////

  //DM 28 Jul 1999 - extract atom coord (for use by std::transform)
  //DM 27 Oct 2000 - return by reference
  inline const RbtCoord& ExtractAtomCoord(RbtAtom* pAtom) {return pAtom->GetCoords();}

  //DM 09 Nov 1999 - Accumulate atomic mass (for use by std::accumulate)
  inline RbtDouble AccumAtomicMass(RbtDouble val, RbtAtom* pAtom) {
    return val+pAtom->GetAtomicMass();
  }

  //DM 11 Nov 1999 - Accumulate mass weighted coords (for use by std::accumulate)
  inline RbtCoord AccumMassWeightedCoords(RbtCoord val, RbtAtom* pAtom) {
    return val+pAtom->GetAtomicMass()*pAtom->GetCoords();
  }

  //Translate an atom by the supplied vector
  class TranslateAtom {
    const RbtVector v;
  public:
    explicit TranslateAtom(const RbtVector& vv) : v(vv) {}
    void operator() (RbtAtom* pAtom) {pAtom->Translate(v);}
  };

  //Translate an atom by the supplied vector, but only if the selection flag is true
  class TranslateAtomIfSelected {
    const RbtVector v;
  public:
    explicit TranslateAtomIfSelected(const RbtVector& vv) : v(vv) {}
    void operator() (RbtAtom* pAtom) {if (pAtom->GetSelectionFlag()) pAtom->Translate(v);}
  };

  //Rotate an atom by the supplied quaternion
  class RotateAtomUsingQuat {
    const RbtQuat q;
  public:
    explicit RotateAtomUsingQuat(const RbtQuat& qq) : q(qq) {}
    void operator() (RbtAtom* pAtom) {pAtom->RotateUsingQuat(q);}
  };

  //Rotate an atom by the supplied quaternion, but only if the selection flag is true
  class RotateAtomUsingQuatIfSelected {
    const RbtQuat q;
  public:
    explicit RotateAtomUsingQuatIfSelected(const RbtQuat& qq) : q(qq) {}
    void operator() (RbtAtom* pAtom) {if (pAtom->GetSelectionFlag()) pAtom->RotateUsingQuat(q);}
  };

  //Select/deselect the atom
  class SelectAtom {
    RbtBool b;
  public:
    explicit SelectAtom(RbtBool bb) : b(bb) {}
    void operator() (RbtAtom* pAtom) {pAtom->SetSelectionFlag(b);}
  };

  //Invert the atom selection flag
  class InvertSelectAtom {
  public:
    explicit InvertSelectAtom() {}
    void operator() (RbtAtom* pAtom) {pAtom->InvertSelectionFlag();}
  };

  //Select the intramolecular flexible interactions to this atom (sets SelectionFlag = true)
  class SelectFlexAtoms {
  public:
    explicit SelectFlexAtoms() {}
    void operator() (RbtAtom* pAtom);
  };

  //Set/unset the cyclic flag
  class CyclicAtom {
    RbtBool b;
  public:
    explicit CyclicAtom(RbtBool bb) : b(bb) {}
    void operator() (RbtAtom* pAtom) {pAtom->SetCyclicFlag(b);}
  };

  ////////////////////////////////////////////
  //Atom list functions (implemented as STL algorithms)
  ////////////////////////////////////////////

  //Unary

  //Generic template version of GetNumAtoms, passing in your own predicate
  template<class Predicate> RbtUInt GetNumAtoms(const RbtAtomList& atomList, const Predicate& pred) {
    return std::count_if(atomList.begin(),atomList.end(),pred);
  }

  //Generic template version of GetAtomList, passing in your own predicate
  template<class Predicate> RbtAtomList GetAtomList(const RbtAtomList& atomList, const Predicate& pred) {
    RbtAtomList newAtomList;
    std::copy_if(atomList.begin(),atomList.end(),std::back_inserter(newAtomList),pred);
    return newAtomList;
  }

  //Generic template version of FindAtom, passing in your own predicate
  template<class Predicate> RbtAtomListIter FindAtom(RbtAtomList& atomList, const Predicate& pred) {
    return std::find_if(atomList.begin(),atomList.end(),pred);
  }

  //Selected atoms
  void SetAtomSelectionFlags(RbtAtomList& atomList,RbtBool bSelected=true);
  void InvertAtomSelectionFlags(RbtAtomList& atomList);//DM 08 Jan 1999

  inline RbtUInt GetNumSelectedAtoms(const RbtAtomList& atomList) {
    return Rbt::GetNumAtoms(atomList,Rbt::isAtomSelected());
  }

  inline RbtAtomList GetSelectedAtomList(const RbtAtomList& atomList) {
    return Rbt::GetAtomList(atomList,Rbt::isAtomSelected());
  }

  //Cyclic atoms
  void SetAtomCyclicFlags(RbtAtomList& atomList,RbtBool bCyclic=true);

  inline RbtUInt GetNumCyclicAtoms(const RbtAtomList& atomList) {
    return Rbt::GetNumAtoms(atomList,Rbt::isAtomCyclic());
  }
  inline RbtAtomList GetCyclicAtomList(const RbtAtomList& atomList) {
    return Rbt::GetAtomList(atomList,Rbt::isAtomCyclic());
  }
  
  //Hydrogen bond acceptor atoms
  inline RbtUInt GetNumHBondAcceptorAtoms(const RbtAtomList& atomList) {
    return Rbt::GetNumAtoms(atomList,Rbt::isAtomHBondAcceptor());
  }
  inline RbtAtomList GetHBondAcceptorAtomList(const RbtAtomList& atomList) {
    return Rbt::GetAtomList(atomList,Rbt::isAtomHBondAcceptor());
  }
  
  //Hydrogen bond donor atoms
  inline RbtUInt GetNumHBondDonorAtoms(const RbtAtomList& atomList) {
    return Rbt::GetNumAtoms(atomList,Rbt::isAtomHBondDonor());
  }
  inline RbtAtomList GetHBondDonorAtomList(const RbtAtomList& atomList) {
    return Rbt::GetAtomList(atomList,Rbt::isAtomHBondDonor());
  }

  //(Formally) charged atoms
  inline RbtUInt GetNumChargedAtoms(const RbtAtomList& atomList) {
    return Rbt::GetNumAtoms(atomList,Rbt::isAtomCharged());
  }
  inline RbtAtomList GetChargedAtomList(const RbtAtomList& atomList) {
    return Rbt::GetAtomList(atomList,Rbt::isAtomCharged());
  }

  //Planar atoms
  inline RbtUInt GetNumPlanarAtoms(const RbtAtomList& atomList) {
    return Rbt::GetNumAtoms(atomList,Rbt::isAtomPlanar());
  }
  inline RbtAtomList GetPlanarAtomList(const RbtAtomList& atomList) {
    return Rbt::GetAtomList(atomList,Rbt::isAtomPlanar());
  }

  //RNA anionic atoms
  //inline RbtUInt GetNumRNAAnionicAtoms(const RbtAtomList& atomList) {
  //  return Rbt::GetNumAtoms(atomList,Rbt::isAtomRNAAnionic());
  //}
  //inline RbtAtomList GetRNAAnionicAtomList(const RbtAtomList& atomList) {
  //  return Rbt::GetAtomList(atomList,Rbt::isAtomRNAAnionic());
  //}

  //Ligand anionic atoms
  //DM 22 Dec 1998
  //inline RbtUInt GetNumLigandAnionicAtoms(const RbtAtomList& atomList) {
  //  return Rbt::GetNumAtoms(atomList,Rbt::isAtomLigandAnionic());
  //}
  //inline RbtAtomList GetLigandAnionicAtomList(const RbtAtomList& atomList) {
  //  return Rbt::GetAtomList(atomList,Rbt::isAtomLigandAnionic());
  //}

  //Ligand cationic atoms
  //DM 22 Dec 1998 -renamed from Get..RNACationic..
  //inline RbtUInt GetNumLigandCationicAtoms(const RbtAtomList& atomList) {
  //  return Rbt::GetNumAtoms(atomList,Rbt::isAtomLigandCationic());
  //}
  //inline RbtAtomList GetLigandCationicAtomList(const RbtAtomList& atomList) {
  //  return Rbt::GetAtomList(atomList,Rbt::isAtomLigandCationic());
  //}

  //Guanidinium carbons
  //inline RbtUInt GetNumGuanidiniumCarbonAtoms(const RbtAtomList& atomList) {
  //  return Rbt::GetNumAtoms(atomList,Rbt::isAtomGuanidiniumCarbon());
  //}
  //inline RbtAtomList GetGuanidiniumCarbonAtomList(const RbtAtomList& atomList) {
  //  return Rbt::GetAtomList(atomList,Rbt::isAtomGuanidiniumCarbon());
  //}

  //BondedAtoms
  //Helper function to return the atom pointer for the "other" atom in the bond
  //i.e. return (bondPair.second) ? bondPair.first->GetAtom2Ptr() : bondPair.first->GetAtom1Ptr();
  RbtAtomPtr GetBondedAtomPtr(pair<RbtBond*,RbtBool> bondBoolPair);

  //This pair of functions is a little different to the above
  //They convert a RbtBondMap to an RbtAtomList. The BondMap is assumed to contain all bonds to a particular atom
  //GetBondedAtomPtr is used to get the "other" atom in the pair.
  RbtUInt GetNumBondedAtoms(const RbtBondMap& bondMap);
  RbtAtomList GetBondedAtomList(const RbtBondMap& bondMap);

  //If spAtom is a regular RbtAtom, these two functions just behave like the two above
  //If spAtom can be dynamically_casted to an RbtPseudoAtom, these return the constituent atom list
  //for the pseudoatom
  RbtUInt GetNumBondedAtoms(const RbtAtom* pAtom);
  RbtAtomList GetBondedAtomList(const RbtAtom* pAtom);


  //Binary

  //Atoms with atomic no = nAtomicNo
  inline RbtUInt GetNumAtomsWithAtomicNo_eq(const RbtAtomList& atomList, RbtInt nAtomicNo) {
    return Rbt::GetNumAtoms(atomList,Rbt::isAtomicNo_eq(nAtomicNo));
  }
  inline RbtAtomList GetAtomListWithAtomicNo_eq(const RbtAtomList& atomList, RbtInt nAtomicNo) {
    return Rbt::GetAtomList(atomList,Rbt::isAtomicNo_eq(nAtomicNo));
  }

  //Atoms with FFType = strFFType
  inline RbtUInt GetNumAtomsWithFFType_eq(const RbtAtomList& atomList, RbtString strFFType) {
    return Rbt::GetNumAtoms(atomList,Rbt::isFFType_eq(strFFType));
  }
  inline RbtAtomList GetAtomListWithFFType_eq(const RbtAtomList& atomList, RbtString strFFType) {
    return Rbt::GetAtomList(atomList,Rbt::isFFType_eq(strFFType));
  }

  //Atoms with AtomName = strAtomName
  inline RbtUInt GetNumAtomsWithAtomName_eq(const RbtAtomList& atomList, RbtString strAtomName) {
    return Rbt::GetNumAtoms(atomList,Rbt::isAtomName_eq(strAtomName));
  }
  inline RbtAtomList GetAtomListWithAtomName_eq(const RbtAtomList& atomList, RbtString strAtomName) {
    return Rbt::GetAtomList(atomList,Rbt::isAtomName_eq(strAtomName));
  }

  //DM 1 Feb 1999
  //Attempts to match atoms in atomList2 with those in atomList1, where the match is performed
  //using isAtom_eq (tests subunit ID, subunit name and atom name for equality) and NOT isAtomPtr_eq
  //(which tests the RbtAtom* pointers themselves for equality, i.e. same objects)
  //Returns list of RbtAtom smart pointers to atoms in atomList1 for which a match is found.
  RbtUInt GetNumMatchingAtoms(const RbtAtomList& atomList1, const RbtAtomList& atomList2);
  RbtAtomList GetMatchingAtomList(const RbtAtomList& atomList1, const RbtAtomList& atomList2);

  //DM 30 Mar 1999
  //These versions match atoms against a full atom name specifier of the form given by GetFullAtomName():
  //<segment name>:<subunit name>_<subunit id>:<atom name>  e.g. A:U_23:O4
  //All fields are optional, and if missing, will match all atoms
  //So, for example, these are all valid full names:
  //N7     (equiv. to ::N7) - matches N7's in all subunits, all segments
  //A::O4                   - matches O4's in all subunits in segment A
  //U:O4 (equiv. to :U:O4)  - matches O4's in subunits named U in all segments
  //U_23:O4                 - matches O4's in subunit U23 in all segments
  //_23:                    - matches all atoms in subunit ID=23 in all segments
  RbtUInt GetNumMatchingAtoms(const RbtAtomList& atomList, const RbtString& strFullName);
  RbtAtomList GetMatchingAtomList(const RbtAtomList& atomList, const RbtString& strFullName);

  //DM 15 Apr 1999 - as above, but match against a list of full atom name specifiers
  //Returns total list (i.e. all matches OR'd). Does not remove duplicates.
  RbtUInt GetNumMatchingAtoms(const RbtAtomList& atomList, const RbtStringList& fullNameList);
  RbtAtomList GetMatchingAtomList(const RbtAtomList& atomList, const RbtStringList& fullNameList);


  /////////////////////////
  //DM 07 Jan 1999
  //Actions on atom lists
  /////////////////////////
  //Translate all atoms in the list by the supplied vector
  inline void TranslateAtoms(const RbtAtomList& atomList, const RbtVector& v) {
    std::for_each(atomList.begin(),atomList.end(),Rbt::TranslateAtom(v));
  }

  //Translate all selected atoms in the list by the supplied vector
  inline void TranslateSelectedAtoms(const RbtAtomList& atomList, const RbtVector& v) {
    std::for_each(atomList.begin(),atomList.end(),Rbt::TranslateAtomIfSelected(v));
  }

  //Rotate all atoms in the list by the supplied quaternion
  inline void RotateAtomsUsingQuat(const RbtAtomList& atomList, const RbtQuat& q) {
    std::for_each(atomList.begin(),atomList.end(),Rbt::RotateAtomUsingQuat(q));
  }

  //Rotate all selected atoms in the list by the supplied quaternion
  inline void RotateSelectedAtomsUsingQuat(const RbtAtomList& atomList, const RbtQuat& q) {
    std::for_each(atomList.begin(),atomList.end(),Rbt::RotateAtomUsingQuatIfSelected(q));
  }

  //Save coords by number for all atoms in the list
  void SaveAtomCoords(const RbtAtomList& atomList, RbtUInt coordNum = 0);
  //Revert to numbered coords for all atoms in the list
  void RevertAtomCoords(const RbtAtomList& atomList, RbtUInt coordNum = 0);

  //DM 9 Nov 1999
  //Returns total atomic mass (molecular weight) for all atoms in the list
  RbtDouble GetTotalAtomicMass(const RbtAtomList& atomList);

  //Returns center of mass of atoms in the list
  RbtCoord GetCenterOfMass(const RbtAtomList& atomList);

  //DM 20 May 1999 - returns the coords of all atoms in the list
  RbtCoordList GetCoordList(const RbtAtomList& atomList);
  //DM 09 Aug 2001 - returns coordList via argument
  void GetCoordList(const RbtAtomList& atomList, RbtCoordList& coordList);

  //Streams an atom list in Quanta CSD file format (for easy color coding of selected atoms in Quanta)
  //nFormat  = 0 => Receptor atom format: "zone 1 # pick O5T = col 2"
  //nFormat != 0 => Ligand atom format "pick N1 .and. segm H = col 2"
  //ostream should have been opened before calling PrintQuantaCSDFormat
  void PrintQuantaCSDFormat(const RbtAtomList& atomList, ostream& s, RbtInt nColor = 2, RbtInt nFormat = 0);

  //DM 25 Feb 1999
  //Modified DM 6 Apr 1999 to operate using the GroupCharge attribute
  //Scans an atom list and sets the group charge to zero on any
  //cationic-ionic pairs which are 1-2 or 1-3 connected
  //(e.g. on -OOC.CR.NH3+ would zero the C and N)
  void RemoveZwitterions(RbtAtomList& atomList);
}

#endif //_RBTATOM_H_
