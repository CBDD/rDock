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

//Abstract base class for molecular file-based output. Provides hooks for
//handling multiple models/conformations being rendered.
//Derived classes must provide a Render method.


#ifndef _RBTBASEMOLECULARFILESINK_H_
#define _RBTBASEMOLECULARFILESINK_H_

#include "RbtBaseFileSink.h"
#include "RbtModel.h"

class RbtBaseMolecularFileSink : public RbtBaseFileSink
{
 public:
  ////////////////////////////////////////
  //Constructors/destructors
  //RbtBaseMolecularFileSink(const char* fileName, RbtModelPtr spModel);
  RbtBaseMolecularFileSink(const RbtString& fileName, RbtModelPtr spModel, RbtBool bUseModelSegmentNames=true);

  virtual ~RbtBaseMolecularFileSink(); //Default destructor


  ////////////////////////////////////////
  //Public methods
  ////////////////
  
  //Derived classes should override
  //Controls whether file type can support the writing of multiple conformations/models to a single file
  virtual RbtBool isMultiConfSupported() {return false;};

  //Get/set the RbtModel that is linked to the sink
  RbtModelPtr GetModel() const;
  //DM 21 Apr 1999 - if bUseModelSegmentNames is true, the model segment names are rendered
  //if bUseModelSegmentNames is false, incremental numeric segment IDs are rendered
  void SetModel(RbtModelPtr spModel, RbtBool bUseModelSegmentNames=false);

  //DM 21 Apr 1999 - change the segment ID counter
  //Set to nSegmentId-1 as the counter gets incremented before each Render
  void SetNextSegmentId(RbtInt nSegmentId) {m_nSegmentId = nSegmentId-1; m_bUseModelSegmentNames=false;}

  //Get/set the multi-conformer status
  //With MultiConf = false, each call to Render should immediately Write the cache to the file
  //With MultiConf = true, each call to Render should continue to fill the cache with atom records
  //The cache is only written by an explicit call to WriteMultiConf, or by the base RbtFileSink destructor
  RbtBool GetMultiConf() const;
  void SetMultiConf(RbtBool bMultiConf);
  //Force writing the file following several multi-conf Renders
  void WriteMultiConf();
  
  RbtModelList GetSolvent() const;
  void SetSolvent(RbtModelList solventList);

 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  void Reset();//Reset the atom, residue and segment numbering

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////

  RbtBaseMolecularFileSink(); //Disable default constructor
  RbtBaseMolecularFileSink(const RbtBaseMolecularFileSink&);//Copy constructor disabled by default
  RbtBaseMolecularFileSink& operator=(const RbtBaseMolecularFileSink&);//Copy assignment disabled by default


 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////
  RbtInt m_nAtomId;//Next available atom number
  RbtInt m_nSubunitId;//Next available subunit number
  RbtInt m_nSegmentId;//Next available segment number
  //DM 21 Apr 1999 - if true, output the atom segment names in the model
  //If false, output incremental segment IDs
  RbtBool m_bUseModelSegmentNames;


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtModelPtr m_spModel;
  RbtModelList m_solventList;
  RbtBool m_bMultiConf;

};

//Useful typedefs
typedef SmartPtr<RbtBaseMolecularFileSink> RbtMolecularFileSinkPtr;//Smart pointer

#endif //_RBTBASEMOLECULARFILESINK_H_
