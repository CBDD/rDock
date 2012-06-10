/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtLigandSiteMapper.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
