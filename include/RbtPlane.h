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

//Simple class for representing 3-D planes in normal form

#ifndef _RBTPLANE_H_
#define _RBTPLANE_H_

#include "RbtCoord.h"

class RbtPlane
{
  ///////////////////////////////////////////////
  // Data members
  ///////////////
  private:
  RbtVector m_vnorm;//The vector normal to the plane (a,b,c)
  RbtDouble m_d;//Distance from plane to origin (in normal form)

  public:
  ///////////////////////////////////////////////
  //Constructors/destructors:
  ///////////////////////////

  //Default constructor (initialise to zero)
  inline RbtPlane() : m_vnorm(0.0,0.0,0.0),m_d(0.0) {}
  
  //Constructor with initial values
  inline RbtPlane(const RbtVector& v, RbtDouble d) : m_vnorm(v), m_d(d) {
    Normalise();//Convert to normal form
  }

  //Constructor taking three coordinates which lie in the plane
  inline RbtPlane(const RbtCoord& c0, const RbtCoord& c1, const RbtCoord& c2) {
    //Determine the vectors from c0 to c1 and c0 to c2
    RbtVector v1 = c1-c0;
    RbtVector v2 = c2-c0;
    m_vnorm = v1.Cross(v2);//Vector normal to the plane
    m_d = -m_vnorm.Dot(c0);
    Normalise();//Convert to normal form
  }

  //Destructor
  virtual ~RbtPlane() {}


  ///////////////////////////////////////////////
  //Operator functions:
  /////////////////////


  ///////////////////////////////////////////////
  //Friend functions:
  ///////////////////

  //Insertion operator
  friend ostream& operator<<(ostream& s, const RbtPlane& plane) {
    return s << plane.m_vnorm.x << "x + "
	     << plane.m_vnorm.y << "y + "
	     << plane.m_vnorm.z << "z + "
	     << plane.m_d << " = 0";
  }

  ///////////////////////////////////////////////
  //Public methods
  ////////////////
  inline RbtDouble D() const {return m_d;}
  inline RbtVector VNorm() const {return m_vnorm;}

  ///////////////////////////////////////////////
  //Private methods
  ////////////////
  inline void Normalise() {
    //Convert to normal form by dividing by +/- vnorm.Length()
    //where sign is opposite the sign of d when d<>0,
    //same as sign of c (vnorm.z) when d==0 and c<>0,
    //and same as sign of b (vnorm.y) otherwise
    RbtDouble l = m_vnorm.Length();
    RbtInt iSign =
      (m_d < 0.0) ? 1 :
      (m_d > 0.0) ? -1 :
      (m_vnorm.z < 0.0) ? -1 :
      (m_vnorm.z > 0.0) ? 1 :
      (m_vnorm.y < 0.0) ? -1 :
      1;
    l *= iSign;
    //Check for divide by zero
    if (l != 0.0) {
      m_vnorm /= l;
      m_d /= l;
    }
  }
};


///////////////////////////////////////////////
// Non-member functions (in Rbt namespace)
//////////////////////////////////////////
namespace Rbt
{
  //Since our planes are always in normal form, we don't need to divide
  //by the magnitude of the normal vector
  //DM 31 Mar 1999 - returned signed distance (not absolute value)
  inline RbtDouble DistanceFromPointToPlane(const RbtCoord& c0, const RbtPlane& p) {
    //return fabs( c0.Dot(p.VNorm()) + p.D() );
    return c0.Dot(p.VNorm()) + p.D();
  }

  //Angle between planes is just the angle between the normal vectors
  inline RbtDouble AngleBetweenPlanes(const RbtPlane& p0, const RbtPlane& p1) {
    return Rbt::Angle(p0.VNorm(),p1.VNorm());
  }  
  
  //Cos(Angle) between planes is just the dot product of the normal vectors
  inline RbtDouble CosAngleBetweenPlanes(const RbtPlane& p0, const RbtPlane& p1) {
    return Rbt::Dot(p0.VNorm(),p1.VNorm());
  }
}

#endif //_RBTPLANE_H_
