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

//Simple class for handling quaternions

#ifndef _RBTQUAT_H_
#define _RBTQUAT_H_

#include <math.h> //for sqrt

#include "RbtCoord.h"

class RbtQuat
{
  ///////////////////////////////////////////////
  // Data members
  ///////////////
  //(leave public so we don't need
  // to write accessor functions)
  public:
  RbtDouble s;//Scalar component
  RbtVector v;//Vector component

  ///////////////////////////////////////////////
  //Constructors/destructors:
  ///////////////////////////
  public:
  
  //Default constructor (corresponds to no rotation (null transformation))
  inline RbtQuat() : s(1.0), v(0.0,0.0,0.0) {}
  
  //Constructor with initial values (vector component passed as 3 RbtDouble's)
  inline RbtQuat(RbtDouble s1, RbtDouble vx, RbtDouble vy, RbtDouble vz) :
    s(s1), v(vx,vy,vz) {}

  //Constructor with initial values (vector component passed as RbtVector)
  //CAUTION: Argument types are the reverse of below
  inline RbtQuat(RbtDouble s1, const RbtVector& v1) :
    s(s1), v(v1) {}
  
  //Constructor accepting a rotation axis and rotation angle (radians).
  //CAUTION: Argument types are the reverse of above
  inline RbtQuat(const RbtVector& axis, RbtDouble phi) {
    RbtDouble halfPhi(0.5*phi);
    s = cos(halfPhi);
    v = sin(halfPhi)*axis.Unit();
  }

  //Destructor
  virtual ~RbtQuat() {}

  //Copy constructor
  inline RbtQuat(const RbtQuat& quat) {
    s = quat.s;
    v = quat.v;
  }


  ///////////////////////////////////////////////
  //Operator functions:
  /////////////////////

  //Copy assignment (*this = coord)
  inline RbtQuat& operator=(const RbtQuat& quat) {
    if (this != &quat) { //beware of self-assignment
      s = quat.s;
      v = quat.v;
    }
    return *this;
  }

  //*this += quat
  inline void operator+=(const RbtQuat& quat) {
    s += quat.s;
    v += quat.v;
  }

  //*this -= quat
  inline void operator-=(const RbtQuat& quat) {
    s -= quat.s;
    v -= quat.v;
  }

  //*this *= const
  inline void operator*=(const RbtDouble& d) {
    s *= d;
    v *= d;
  }

  //*this /= const
  inline void operator/=(const RbtDouble& d) {
    s /= d;
    v /= d;
  }

  ///////////////////////////////////////////////
  //Friend functions:
  ///////////////////

  //Insertion operator
  friend ostream& operator<<(ostream& s, const RbtQuat& quat) {
    return s << quat.s << "," << quat.v;
  }

  //Equality
  inline friend RbtBool operator== (const RbtQuat& quat1, const RbtQuat& quat2) {
    return quat1.s == quat2.s && quat1.v == quat2.v;
  }

  //Non-equality
  inline friend RbtBool operator!= (const RbtQuat& quat1, const RbtQuat& quat2) {
    return quat1.s != quat2.s || quat1.v != quat2.v;
  }

  //Addition
  inline friend RbtQuat operator+(const RbtQuat& quat1, const RbtQuat& quat2) {
    return RbtQuat(quat1.s + quat2.s, quat1.v + quat2.v);
  }

  //Subtraction
  inline friend RbtQuat operator-(const RbtQuat& quat1, const RbtQuat& quat2) {
    return RbtQuat(quat1.s - quat2.s, quat1.v - quat2.v);
  }

  //Negation
  inline friend RbtQuat operator-(const RbtQuat& quat) {
    return RbtQuat(-quat.s, -quat.v);
  }

  //Multiplication (non-commutative)
  inline friend RbtQuat operator*(const RbtQuat& quat1, const RbtQuat& quat2) {
    return RbtQuat(quat1.s*quat2.s - (quat1.v).Dot(quat2.v), quat1.s*quat2.v + quat1.v*quat2.s + (quat1.v).Cross(quat2.v));
  }

  //Scalar product
  inline friend RbtQuat operator*(const RbtQuat& quat, const RbtDouble& d) {
    return RbtQuat(quat.s*d, quat.v*d);
  }

  //Scalar product
  inline friend RbtQuat operator*(const RbtDouble& d, const RbtQuat& quat) {
    return RbtQuat(quat.s*d, quat.v*d);
  }

  //Scalar division
  inline friend RbtQuat operator/(const RbtQuat& quat, const RbtDouble& d) {
    return RbtQuat(quat.s/d, quat.v/d);
  }

  ///////////////////////////////////////////////
  //Public methods
  ////////////////

  //Returns magnitude of quat
  //Member function (L = Q1.Length())
  inline RbtDouble Length() const {
    return sqrt(s*s + v.x*v.x + v.y*v.y + v.z*v.z);
  }

  //Returns unit quat
  //Member function (U1 = Q1.Unit())
  inline RbtQuat Unit() const {
    RbtDouble l = Length();
    return (l>0) ? *this/l : *this;
  }

  //Dot product member function (D = Q1.Dot(Q2))
  inline RbtDouble Dot(const RbtQuat& quat) const {
    return s*quat.s + v.Dot(quat.v);
  }

  //Returns conjugate
  //Member function (Q1* = Q1.Conj())
  inline RbtQuat Conj() const {
    return RbtQuat(s,-v);
  }
  
  //Returns Q W Q* where Q=S+V is a quaternion and W is a 3-D vector
  //Q W Q* = (S^2 - V.V)W + 2S(VxW) + 2V(V.W)
  //where . = Dot product and x = Cross product
  //The end result is a 3-D rotation, assuming q has been obtained from RotationQuat
  inline RbtCoord Rotate(const RbtCoord& w) const {
    return (s*s - v.Dot(v) ) * w
            + 2 * s * v.Cross(w)
            + 2 * v * v.Dot(w);
  }
};

//Useful typedefs
typedef vector<RbtQuat> RbtQuatList;//Vector of quats
typedef RbtQuatList::iterator RbtQuatListIter;
typedef RbtQuatList::const_iterator RbtQuatListConstIter;

///////////////////////////////////////////////
// Non-member functions (in Rbt namespace)
//////////////////////////////////////////
namespace Rbt
{
  //Returns magnitude of quat (L = Length(Q1))
  inline RbtDouble Length(const RbtQuat& quat) {
    return quat.Length();
  }

  //Returns unit quat (U1 = Unit(Q1))
  inline RbtQuat Unit(const RbtQuat& quat) {
    return quat.Unit();
  }

  //Dot product (D = Dot(Q1,Q2))
  inline RbtDouble Dot(const RbtQuat& quat1, const RbtQuat& quat2) {
    return quat1.Dot(quat2);
  }

  //Returns conjugate (Q1* = Conj(Q1))
  inline RbtQuat Conj(const RbtQuat& quat) {
    return quat.Conj();
  }
}

#endif //_RBTQUAT_H_
