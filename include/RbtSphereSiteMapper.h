/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtSphereSiteMapper.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//2-sphere site mapper function class

#ifndef _RBTSPHERESITEMAPPER_H_
#define _RBTSPHERESITEMAPPER_H_

#include "RbtSiteMapper.h"

class RbtSphereSiteMapper : public RbtSiteMapper {
 public:
  //Static data member for class type
  static RbtString _CT;
  //Parameter names
  static RbtString _VOL_INCR;
  static RbtString _SMALL_SPHERE;
  static RbtString _LARGE_SPHERE;
  static RbtString _GRIDSTEP;
  static RbtString _CENTER;
  static RbtString _RADIUS;
  static RbtString _MIN_VOLUME;
  static RbtString _MAX_CAVITIES;
  
  RbtSphereSiteMapper(const RbtString& strName = "SPHERE_MAPPER");
  virtual ~RbtSphereSiteMapper();

  //Override RbtSiteMapper pure virtual
  //This is the function which actually does the mapping
  virtual RbtCavityList operator() ();
};

//Useful typedefs
typedef SmartPtr<RbtSphereSiteMapper> RbtSphereSiteMapperPtr;//Smart pointer

#endif //_RBTSPHERESITEMAPPER_H_
