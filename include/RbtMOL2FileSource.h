/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtMOL2FileSource.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

#ifndef _RBTMOL2FILESOURCE_H_
#define _RBTMOL2FILESOURCE_H_

#include "RbtBaseMolecularFileSource.h"
#include "RbtElementFileSource.h"
#include "RbtParameterFileSource.h"

//Simple struct for holding MOL2 SUBSTRUCTURE records
class RbtMOL2Substructure {
 public:
  RbtMOL2Substructure() : root_atom(0), chain("UNK"), sub_type("UNK") {};
  RbtMOL2Substructure(const RbtString& name, RbtInt rr, const RbtString& cc, const RbtString& ss) :
    subst_name(name), root_atom(rr), chain(cc), sub_type(ss) {};
  RbtString GetName() const {return subst_name;};
  RbtInt GetRoot() const {return root_atom;};
  RbtString GetChain() const {return chain;};
  RbtString GetType() const {return sub_type;};
 private:
  RbtString subst_name;
  RbtInt root_atom;
  RbtString chain;
  RbtString sub_type;
};

//A MOL2SubstructureMap stores MOL2Substructures keyed by SUBSTRUCTURE ID (subst_id)
typedef map<RbtInt,RbtMOL2Substructure> RbtMOL2SubstructureMap;
typedef RbtMOL2SubstructureMap::iterator RbtMOL2SubstructureMapIter;
typedef RbtMOL2SubstructureMap::const_iterator RbtMOL2SubstructureMapConstIter;

class RbtMOL2FileSource : public RbtBaseMolecularFileSource {
 public:
  static RbtString _CT;
  // record delimiter strings
  static RbtString _TRIPOS_DELIM; 
  static RbtInt _TRIPOS_DELIM_SIZE;
  static RbtString _IDS_MOL2_RECDELIM;

  RbtMOL2FileSource(const RbtString& fileName, RbtBool bImplHydrogens=true);
  ~RbtMOL2FileSource();
  RbtBool isTitleListSupported() {return true;}; 
  RbtBool isAtomListSupported() {return true;}; 
  RbtBool isCoordinatesSupported() {return true;};
  RbtBool isBondListSupported() {return true;};
  RbtBool isDataSupported() {return false;};
  
 protected:
  virtual void Parse() throw (RbtError);
  
 private:
  RbtMOL2FileSource();
  RbtMOL2FileSource( const RbtMOL2FileSource&);
  RbtMOL2FileSource& operator=(const RbtMOL2FileSource&);
  
  //Helper functions to clean up atom attributes after parsing the file
  void SetupAtomParams() throw(RbtError);
  void FixImplicitHydrogenCount();
  void FixHybridState();
  void FixTriposTypes();
  void SetupVdWRadii();
  void RemoveNonPolarHydrogens();

  //Splits subst_name into ID and name components (e.g. ALA27 => ALA and 27)
  void GetSSIDandName(const RbtString& subst_name, RbtInt subst_id, RbtString& sID, RbtString& sName);

  // parser functions (switched by pointers)
  void ParseRecordMOLECULE(const RbtString& aLine);
  void ParseRecordATOM(const RbtString& aLine);
  void ParseRecordBOND(const RbtString& aLine);
  void ParseRecordSUBSTRUCTURE(const RbtString& aLine);
  void ParseRecordUNSUPPORTED(const RbtString& aLine);
  // get delimiter tag
  RbtString GetMOL2Tag(const RbtString& aLine);
  // parse/tokenize second line (number of  atoms, bonds etc) in MOLECULE
  void ParseCountFields(const RbtString& aLine);
  void Tokenize(const RbtString& aString, RbtStringList& aTokensBuf);

  int m_NL;//No. of lines
  // data fields in MOLECULE
  RbtInt nAtoms;// number of atoms
  RbtInt nBonds;//           bonds
  RbtInt nSubstructures;//substructures
  RbtInt nFeatures;//features
  RbtInt nSets;//sets

  // helper data structures used to lookup chain names for each substructure
  //and to sort atom list into a more sensible order (ordered by substructure ID)
  RbtMOL2SubstructureMap m_ssInfo;
  map<RbtInt,RbtAtomList> m_ssAtoms;
  
  RbtElementFileSourcePtr m_spElementData;	//Elemental data source for vdW radii etc
  RbtParameterFileSourcePtr m_spParamSource;	//Parameter file source for partial charges etc
  RbtTriposAtomType m_typer;//Tripos atom typing object
  RbtBool m_bImplHydrogens;
};

// usual smart pointer
typedef SmartPtr<RbtMOL2FileSource> RbtMOL2FileSourcePtr;
// function pointer (callback) for Parse()
typedef void (RbtMOL2FileSource::*fcPtr)(const RbtString&);

#endif // _RBTMOL2FILESOURCE_H_

