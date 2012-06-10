/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtMdlFileSource.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Handles retrieval of molecular info from MDL Mol and SD files.

#ifndef _RBTMDLFILESOURCE_H_
#define _RBTMDLFILESOURCE_H_

#include "RbtBaseMolecularFileSource.h"
#include "RbtElementFileSource.h"

const RbtString IDS_MDL_RECDELIM = "$$$$";

class RbtMdlFileSource : public RbtBaseMolecularFileSource
{
 public:
  //Constructors
  //RbtMdlFileSource(const char* fileName);
  RbtMdlFileSource(const RbtString& fileName, RbtBool bPosIonisable=true, RbtBool bNegIonisable=true, RbtBool bImplHydrogens=true);

  //Default destructor
  virtual ~RbtMdlFileSource();

  ////////////////////////////////////////
  //Override public methods from RbtBaseMolecularDataSource
  virtual RbtBool isTitleListSupported() {return true;};
  virtual RbtBool isAtomListSupported() {return true;};
  virtual RbtBool isCoordinatesSupported() {return true;};
  virtual RbtBool isBondListSupported() {return true;};
  //DM 12 May 1999 - support for data records
  virtual RbtBool isDataSupported() {return true;};

 protected:
  //Pure virtual in RbtBaseFileSource - needs to be defined here
  virtual void Parse() throw (RbtError);

 private:
  //Private methods
  RbtMdlFileSource();//Disable default constructor
  RbtMdlFileSource(const RbtMdlFileSource&);//Copy constructor disabled by default
  RbtMdlFileSource& operator=(const RbtMdlFileSource&);//Copy assignment disabled by default

  //Sets up all the atomic attributes that are not explicitly stored in the MDL file
  void SetupAtomParams() throw (RbtError);

  //Helper functions for SetupAtomParams()

  //Sets hybridisation state and checks for Nsp3 next to sp2 special case
  void SetupHybridState() throw (RbtError);

  //Checks for valency within range
  //Adds implicit hydrogens to C,N,S
  //Defines vdW radius, correcting for extended atoms and H-bond donor hydrogens
  //Defines formal "force-field" type string
  void SetupTheRest() throw (RbtError);
    
  void SetupPosIonisableGroups();
  void SetupNegIonisableGroups();

  void AddHydrogen(RbtAtomPtr spAtom) throw (RbtError);


  //SetupIonicGroups sets the "group charge" attribute of the atoms. The idea is that the "group charge" is
  //file-format independent, and eliminates the need for the scoring function setup to have any knowledge of
  //the input file format. For example, PSF files represent COO- as OC-C-OC (both oxygens charged) whereas
  //SD files represent it as the formal resonance structure O=C-O-. For RiboDock we actually want the negative
  //charge on the central carbon!
  //
  //It is envisaged each molecular file source will have a version of SetupIonicGroups to convert from the native
  //representation to the RiboDock representation. The FormalCharge and PartialCharge attributes should be left
  //untouched to allow the model to be rewritten back in the same format.
  void SetupIonicGroups();
  void SetupNSP3Plus();//Helper function to set up protonated amines
  void SetupNSP2Plus();//Helper function to set up guanidines, imidazoles, amidines
  void SetupOTRIMinus();//Helper function to set up deprotonated acids

  //DM 4 June 1999 - Search for the separate fragments in the record
  //and set distinct segment names for each one
  void SetupSegmentNames();

  //DM 2 Aug 1999 - remove all non-polar hydrogens
  void RemoveNonPolarHydrogens();

  //DM 8 Feb 2000 - setup atom and bond cyclic flags (previously in RbtModel)
  //Set the atom and bond cyclic flags for all atoms and bonds in the source
  void SetCyclicFlags();


  //Private data
  RbtElementFileSourcePtr m_spElementData; //Elemental data source for vdW radii etc
  RbtBool m_bPosIonisable;//If true, protonate Nsp3, and Nsp2 in imidazole and guanidinium
  RbtBool m_bNegIonisable;//If true, deprotonate Otri-Hsp3
  RbtBool m_bImplHydrogens;//If true, non-hydrogen bonding hydrogens are removed
};

//useful typedefs
typedef SmartPtr<RbtMdlFileSource> RbtMdlFileSourcePtr;//Smart pointer

#endif //_RBTMDLFILESOURCE_H_
