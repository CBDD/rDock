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
