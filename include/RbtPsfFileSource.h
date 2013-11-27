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

//Handles retrieval of molecular info from Charmm PSF files.

#ifndef _RBTPSFFILESOURCE_H_
#define _RBTPSFFILESOURCE_H_

#include "RbtBaseMolecularFileSource.h"
#include "RbtCharmmDataSource.h"
#include "RbtParameterFileSource.h"
#include "RbtElementFileSource.h"

class RbtPsfFileSource : public RbtBaseMolecularFileSource
{
 public:
  //Constructors
  RbtPsfFileSource(const RbtString& fileName,
		   const RbtString& strMassesFile = Rbt::GetRbtFileName("data","masses.rtf"),
		   RbtBool bImplHydrogens=true);

  //Default destructor
  virtual ~RbtPsfFileSource();

  ////////////////////////////////////////
  //Override public methods from RbtBaseMolecularDataSource
  virtual RbtBool isTitleListSupported() {return true;};
  virtual RbtBool isAtomListSupported() {return true;};
  virtual RbtBool isCoordinatesSupported() {return false;};
  virtual RbtBool isBondListSupported() {return true;};

 protected:
  //Pure virtual in RbtBaseFileSource - needs to be defined here
  virtual void Parse() throw (RbtError);

 private:
  //Private methods
  RbtPsfFileSource();//Disable default constructor
  RbtPsfFileSource(const RbtPsfFileSource&);//Copy constructor disabled by default
  RbtPsfFileSource& operator=(const RbtPsfFileSource&);//Copy assignment disabled by default

  //Sets up all the atomic attributes that are not explicitly stored in the PSF file
  void SetupAtomParams() throw (RbtError);
  void SetupVdWRadii() throw (RbtError);
  void SetupPartialIonicGroups();
  void RemoveNonPolarHydrogens();
  //Is atom in same substructure
  class isSS_eq : public std::unary_function<RbtAtom*,RbtBool> {
    const RbtAtom* a1;
  public:
    explicit isSS_eq(const RbtAtom* aa) : a1(aa) {}
    RbtBool operator() (const RbtAtom* a2) const
    {
      return
	(a1->GetSubunitName() == a2->GetSubunitName()) &&
	(a1->GetSubunitId() == a2->GetSubunitId()) &&
	(a1->GetSegmentName() == a2->GetSegmentName());
    }
  };

  //Private data
  RbtElementFileSourcePtr m_spElementData; //Elemental data source for vdW radii etc
  RbtCharmmDataSourcePtr m_spCharmmData;//Charmm data source for fftype lookup
  RbtParameterFileSourcePtr m_spParamSource;//Parameter file source for partial charges etc
  RbtBool m_bImplHydrogens;//If true, non-hydrogen bonding hydrogens are removed
};

//useful typedefs
typedef SmartPtr<RbtPsfFileSource> RbtPsfFileSourcePtr;//Smart pointer

#endif //_RBTPSFFILESOURCE_H_
