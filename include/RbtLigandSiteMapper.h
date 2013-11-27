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

//Reference ligand site mapper function class

#ifndef _RBTLIGANDSITEMAPPER_H_
#define _RBTLIGANDSITEMAPPER_H_

#include "RbtSiteMapper.h"

class RbtLigandSiteMapper : public RbtSiteMapper {
 public:
  //Static data member for class type
  static RbtString _CT;
  //Parameter names
  static RbtString _REF_MOL;
  static RbtString _VOL_INCR;
  static RbtString _SMALL_SPHERE;
  static RbtString _GRIDSTEP;
  static RbtString _RADIUS;
  static RbtString _MIN_VOLUME;
  static RbtString _MAX_CAVITIES;
  
  RbtLigandSiteMapper(const RbtString& strName = "LIGAND_MAPPER");
  virtual ~RbtLigandSiteMapper();

  //Override RbtSiteMapper pure virtual
  //This is the function which actually does the mapping
  virtual RbtCavityList operator() ();
};

//Useful typedefs
typedef SmartPtr<RbtLigandSiteMapper> RbtLigandSiteMapperPtr;//Smart pointer

#endif //_RBTLIGANDSITEMAPPER_H_
