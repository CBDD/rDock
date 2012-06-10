/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtMdlFileSink.h#4 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Class for writing RbtModel's to MDL SD and MOL files

#ifndef _RBTMDLFILESINK_H_
#define _RBTMDLFILESINK_H_

#include "RbtBaseMolecularFileSink.h"
#include "RbtElementFileSource.h"

//DM 19 June 2006
//Map to keep track of logical atom IDs as they are rendered to file
//This avoids the need for the atoms in the RbtModel to have consecutive
//atom IDs. The rendering of enabled solvent models is now supported
//correctly.
//Key = RbtAtom* pointer
//Value = logical atom ID as rendered to file
//Logical atom IDs are consecutive for all atoms rendered
//The actual atom ID (RbtAtom::GetAtomId()) is unchanged
typedef map<RbtAtom*,RbtUInt,Rbt::RbtAtomPtrCmp_Ptr> RbtAtomIdMap;

class RbtMdlFileSink : public RbtBaseMolecularFileSink
{
 public:
  ////////////////////////////////////////
  //Constructors/destructors
  RbtMdlFileSink(const RbtString& fileName, RbtModelPtr spModel);

  virtual ~RbtMdlFileSink(); //Default destructor


  ////////////////////////////////////////
  //Public methods
  ////////////////
  //
  //
  ////////////////////////////////////////
  //Override public methods from RbtBaseFileSink
  virtual void Render() throw (RbtError);

 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  void RenderAtomList(const RbtAtomList& atomList);
  void RenderBondList(const RbtBondList& bondList);
  void RenderData(const RbtStringVariantMap& dataMap);
  				
  RbtMdlFileSink(); //Disable default constructor
  RbtMdlFileSink(const RbtMdlFileSink&);//Copy constructor disabled by default
  RbtMdlFileSink& operator=(const RbtMdlFileSink&);//Copy assignment disabled by default
  

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////

 private:
  ////////////////////////////////////////
  //Private data
  //////////////
 RbtElementFileSourcePtr m_spElementData;
 //DM 27 Apr 1999 - default behaviour is for the first Render to overwrite any existing file
 //then subsequent Renders to append.
 RbtBool m_bFirstRender;
 RbtAtomIdMap m_atomIdMap;//Keep track of logical atom IDs as rendered to file
};

//Useful typedefs
typedef SmartPtr<RbtMdlFileSink> RbtMdlFileSinkPtr;//Smart pointer

namespace Rbt
{
  void RenumberScaffold(RbtAtomList& atomList);
}

#endif //_RBTMDLFILESINK_H_
