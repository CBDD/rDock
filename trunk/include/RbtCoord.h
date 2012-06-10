/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtCoord.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Simple class for handling cartesian coords (points) and 3-d vectors
//The same class is used to represent both points and vectors, although the
//RbtVector typedef can be used in argument lists to distinguish the two cases

#ifndef _RBTCOORD_H_
#define _RBTCOORD_H_

#include <algorithm> //for min,max
#include <numeric> //for accumulate
#include <math.h> //for sqrt

#include "RbtConfig.h"

extern istream& eatSeps(istream& is);

class RbtCoord
{
  ///////////////////////////////////////////////
  // Data members
  ///////////////
  //(leave public so we don't need
  // to write accessor functions)
  public:
  RbtDouble x;
  RbtDouble y;
  RbtDouble z;

  ///////////////////////////////////////////////
  //Constructors/destructors:
  ///////////////////////////
  public:

  //Default constructor (initialise to zero)
  inline RbtCoord() : x(0.0), y(0.0), z(0.0) {}
  
  //Constructor with initial values
  inline RbtCoord(RbtDouble x1, RbtDouble y1, RbtDouble z1) :
    x(x1), y(y1), z(z1) {}

  //Destructor
  virtual ~RbtCoord() {}

  //Copy constructor
  inline RbtCoord(const RbtCoord& coord) {
    x = coord.x;
    y = coord.y;
    z = coord.z;
  }

  //DM 19 Jul 2000 - Read,Write methods to read/write coords to binary streams
  inline ostream& Write(ostream& ostr) const {
    ostr.write( (const char*) &x, sizeof(x));
    ostr.write( (const char*) &y, sizeof(y));
    ostr.write( (const char*) &z, sizeof(z));
    return ostr;
  }
  inline istream& Read(istream& istr) {
    istr.read( (char*) &x, sizeof(x));
    istr.read( (char*) &y, sizeof(y));
    istr.read( (char*) &z, sizeof(z));
    return istr;
  }

  ///////////////////////////////////////////////
  //Operator functions:
  /////////////////////

  //Copy assignment (*this = coord)
  inline RbtCoord& operator=(const RbtCoord& coord) {
    if (this != &coord) { //beware of self-assignment
      x = coord.x;
      y = coord.y;
      z = coord.z;
    }
    return *this;
  }

  //Copy assignment (*this = const)
  inline RbtCoord& operator=(const RbtDouble& d) {
    x = d;
    y = d;
    z = d;
    return *this;
  }

  //*this += coord
  inline void operator+=(const RbtCoord& coord) {
    x += coord.x;
    y += coord.y;
    z += coord.z;
  }

  //*this += const
  inline void operator+=(const RbtDouble& d) {
    x += d;
    y += d;
    z += d;
  }

  //*this -= coord
  inline void operator-=(const RbtCoord& coord) {
    x -= coord.x;
    y -= coord.y;
    z -= coord.z;
  }

  //*this -= const
  inline void operator-=(const RbtDouble& d) {
    x -= d;
    y -= d;
    z -= d;
  }

  //*this *= const
  inline void operator*=(const RbtDouble& d) {
    x *= d;
    y *= d;
    z *= d;
  }

  //*this /= const
  inline void operator/=(const RbtDouble& d) {
    x /= d;
    y /= d;
    z /= d;
  }

  ///////////////////////////////////////////////
  //Friend functions:
  ///////////////////

  //Insertion (output) operator
  friend ostream& operator<<(ostream& s, const RbtCoord& coord) {
    return s << "(" << coord.x << "," << coord.y << "," << coord.z <<")";
  }

  //Input operator
  //DM 15 Apr 1999
  //Valid formats (actually the comma can be any single non-whitespace character):
  // x,y,z
  // (x,y,z)
  // BGD 06 Nov 2002
  // Trying to make it more versatil. The separator now is any
  // amount of white space or commas.
  // x   y,z
  // (x   y, z  )
  // uses the auxiliary function eatSeps
  
  friend istream& operator>>(istream& s, RbtCoord& coord) {
    RbtDouble x=0, y=0, z=0;
    char c = 0;
    s >> c;
    if (c == '(') {
      s >> x;
	  eatSeps(s);
      s >> y ;
      eatSeps(s);
	  s >> z ;
      eatSeps(s);
	  s >> c;
      if (c != ')') s.clear(ios_base::badbit);//Missing closing bracket
    }
    else {
      s.putback(c);
      s >> x;
	  eatSeps(s);
      s >> y ;
      eatSeps(s);
	  s >> z ;
    }
    if (s) coord = RbtCoord(x,y,z);
    return s;
  }

  //Equality
  inline friend RbtBool operator== (const RbtCoord& coord1, const RbtCoord& coord2) {
    return coord1.x == coord2.x && coord1.y == coord2.y && coord1.z == coord2.z;
  }

  //Non-equality
  inline friend RbtBool operator!= (const RbtCoord& coord1, const RbtCoord& coord2) {
    return coord1.x != coord2.x || coord1.y != coord2.y || coord1.z != coord2.z;
  }

  //Greater than
  inline friend RbtBool operator> (const RbtCoord& coord1, const RbtCoord& coord2) {
    return coord1.x > coord2.x && coord1.y > coord2.y && coord1.z > coord2.z;
  }

  //Greater than or equal
  inline friend RbtBool operator>= (const RbtCoord& coord1, const RbtCoord& coord2) {
    return coord1.x >= coord2.x && coord1.y >= coord2.y && coord1.z >= coord2.z;
  }

  //Less than
  inline friend RbtBool operator< (const RbtCoord& coord1, const RbtCoord& coord2) {
    return coord1.x < coord2.x && coord1.y < coord2.y && coord1.z < coord2.z;
  }

  //Less than or equal
  inline friend RbtBool operator<= (const RbtCoord& coord1, const RbtCoord& coord2) {
    return coord1.x <= coord2.x && coord1.y <= coord2.y && coord1.z <= coord2.z;
  }

  //Addition (coord + coord)
  inline friend RbtCoord operator+(const RbtCoord& coord1, const RbtCoord& coord2) {
    return RbtCoord(coord1.x + coord2.x, coord1.y + coord2.y, coord1.z + coord2.z);
  }

  //Addition (coord + const)
  inline friend RbtCoord operator+(const RbtCoord& coord1, RbtDouble d) {
    return RbtCoord(coord1.x + d, coord1.y + d, coord1.z + d);
  }

  //Addition (const + coord)
  inline friend RbtCoord operator+(RbtDouble d, const RbtCoord& coord1) {
    return RbtCoord(d + coord1.x, d+ coord1.y, d + coord1.z);
  }

  //Subtraction (coord - coord)
  inline friend RbtCoord operator-(const RbtCoord& coord1, const RbtCoord& coord2) {
    return RbtCoord(coord1.x - coord2.x, coord1.y - coord2.y, coord1.z - coord2.z);
  }

  //Subtraction (coord - const)
  inline friend RbtCoord operator-(const RbtCoord& coord1, RbtDouble d) {
    return RbtCoord(coord1.x - d, coord1.y - d, coord1.z - d);
  }

  //Subtraction (const - coord)
  inline friend RbtCoord operator-(RbtDouble d, const RbtCoord& coord1) {
    return RbtCoord(d - coord1.x, d - coord1.y, d - coord1.z);
  }

  //Negation
  inline friend RbtCoord operator-(const RbtCoord& coord) {
    return RbtCoord(-coord.x, -coord.y, -coord.z);
  }

  //Dot product
  //DM 7 Jan 1999 - operator* replaced by Dot(coord) member function
  //and by Dot(coord1,coord2) non-member function
  //inline friend RbtDouble operator*(const RbtCoord& coord1, const RbtCoord& coord2) {
  //  return coord1.x*coord2.x + coord1.y*coord2.y + coord1.z*coord2.z;
  //}

  //Scalar product (coord * const)
  inline friend RbtCoord operator*(const RbtCoord& coord, const RbtDouble& d) {
    return RbtCoord(coord.x*d, coord.y*d, coord.z*d);
  }

  //Scalar product (const * coord)
  inline friend RbtCoord operator*(const RbtDouble& d, const RbtCoord& coord) {
    return RbtCoord(coord.x*d, coord.y*d, coord.z*d);
  }

  //Scalar product (coord * coord : component-wise multiplication)
  inline friend RbtCoord operator*(const RbtCoord& coord1, const RbtCoord& coord2) {
    return RbtCoord(coord1.x*coord2.x, coord1.y*coord2.y, coord1.z*coord2.z);
  }

  //Scalar division (coord / const)
  inline friend RbtCoord operator/(const RbtCoord& coord, const RbtDouble& d) {
    return RbtCoord(coord.x/d, coord.y/d, coord.z/d);
  }


  ///////////////////////////////////////////////
  //Public methods
  ////////////////

  //Mostly applicable to vectors (rather than coords)

  //Returns square of magnitude of vector (or distance from origin)
  inline RbtDouble Length2() const {
    return x*x + y*y + z*z;
  }

  //Returns magnitude of vector (or distance from origin)
  inline RbtDouble Length() const {
    return sqrt(Length2());
  }

  //Returns unit vector in same direction
  //Member function (V2 = V1.Unit())
  //Checks for zero length
  inline RbtCoord Unit() const {
    RbtDouble l = Length();
    return (l>0) ? *this/l : *this;
  }

  //Cross product member function (V3 = V1.Cross(V2))
  inline RbtCoord Cross(const RbtCoord& v2) const {
    return RbtCoord( y*v2.z - v2.y*z,
             z*v2.x - v2.z*x,
             x*v2.y - v2.x*y);
  }

  //Dot product member function (D = V1.Dot(V2)
  inline RbtDouble Dot(const RbtCoord& v2) const {
    return x*v2.x + y*v2.y + z*v2.z;
  }
};


//Useful typedefs

//DM 4 Mar 1999 - to make it easier to distinguish arguments that are
//treated as -3D vectors from those are treated as 3-D points, use this typedef
typedef RbtCoord RbtVector;

typedef vector<RbtCoord> RbtCoordList;//Vector of coords
typedef RbtCoordList::iterator RbtCoordListIter;
typedef RbtCoordList::const_iterator RbtCoordListConstIter;

//DM 19/1/99 typedefs for a map of key=RbtInt, value=RbtCoord
//typedef map<RbtString,RbtCoord> RbtStringCoordMap;
//typedef RbtStringCoordMap::iterator RbtStringCoordMapIter;
//typedef RbtStringCoordMap::const_iterator RbtStringCoordMapConstIter;
//DM 08 Feb 1999 - Changed from key=string to key=int
typedef map<RbtUInt,RbtCoord> RbtUIntCoordMap;
typedef RbtUIntCoordMap::iterator RbtUIntCoordMapIter;
typedef RbtUIntCoordMap::const_iterator RbtUIntCoordMapConstIter;


///////////////////////////////////////////////
// Non-member functions (in Rbt namespace)
//////////////////////////////////////////
namespace Rbt
{
  //////////////////
  //VECTOR FUNCTIONS
  //////////////////

  //Returns square of magnitude of vector
  inline RbtDouble Length2(const RbtVector& v1) {
    return v1.Length2();
  }

  //Returns square of distance between two coords
  inline RbtDouble Length2(const RbtCoord& c1, const RbtCoord& c2) {
    return Rbt::Length2(c2-c1);
  }

  //Returns magnitude of vector
  inline RbtDouble Length(const RbtVector& v1) {
    return v1.Length();
  }

  //Returns distance between two coords
  inline RbtDouble Length(const RbtCoord& c1, const RbtCoord& c2) {
    return Rbt::Length(c2-c1);
  }

  //Returns unit vector in same direction
  //V2 = Unit(V1)
  inline RbtVector Unit(const RbtVector& v1) {
    return v1.Unit();
  }

  //Cross product (V3 = Rbt::Cross(V1,V2))
  inline RbtVector Cross(const RbtVector& v11, const RbtVector& v2) {
    return v11.Cross(v2);
  }

  //Dot product (D = Rbt::Dot(V1,V2))
  inline RbtDouble Dot(const RbtVector& v1, const RbtVector& v2) {
    return v1.Dot(v2);
  }

  //Returns angle formed between 2 vectors
  inline RbtDouble Angle(const RbtVector& v1, const RbtVector& v2) {
    //Calculate the product of the lengths
    RbtDouble d1d2 = v1.Length()*v2.Length();

    //Calculate the sin and cos
    RbtDouble cos_theta = v1.Dot(v2)/d1d2;
    RbtDouble sin_theta = (v1.Cross(v2)).Length()/d1d2;

    //Get theta and convert to degrees
    RbtDouble theta = atan2(sin_theta,cos_theta);
    return theta*180.0/M_PI;
  }

  //Returns angle formed between 3 coords coord1..coord2..coord3
  inline RbtDouble Angle(const RbtCoord& coord1, const RbtCoord& coord2, const RbtCoord& coord3) {
    //Determine vectors between coords and call vector version of Angle()
    return Rbt::Angle(coord1-coord2,coord3-coord2);
  }

  //DM 7 June 1999
  //Returns dihedral formed between 3 vectors
  inline RbtDouble Dihedral(const RbtVector& v1, const RbtVector& v2, const RbtVector& v3) {
    //Calculate the cross products
    RbtVector A = v1.Cross(v2);
    RbtVector B = v2.Cross(v3);
    RbtVector C = v2.Cross(A);
    //Calculate the distances
    RbtDouble rA = A.Length();
    RbtDouble rB = B.Length();
    RbtDouble rC = C.Length();
    //Calculate the sin and cos
    RbtDouble cos_phi = A.Dot(B)/(rA*rB);
    RbtDouble sin_phi = C.Dot(B)/(rC*rB);
    //Get phi and convert to degrees
    RbtDouble phi = -atan2(sin_phi,cos_phi);
    return phi*180.0/M_PI;
  }

  //Returns dihedral formed between 4 coords c1..c2..c3..c4
  inline RbtDouble Dihedral(const RbtCoord& c1, const RbtCoord& c2, const RbtCoord& c3, const RbtCoord& c4) {
    //Determine vectors between coords and call vector version of Dihedral()
    return Rbt::Dihedral(c1-c2,c2-c3,c3-c4);
  }

  //////////////////////
  //COORDINATE FUNCTIONS
  //////////////////////

  //DM 28 Jul 1999 - functions for use by STL transform
  inline RbtDouble ExtractXCoord(const RbtCoord& c) {return c.x;}
  inline RbtDouble ExtractYCoord(const RbtCoord& c) {return c.y;}
  inline RbtDouble ExtractZCoord(const RbtCoord& c) {return c.z;}

  //Returns minimum of two coords (component-wise minimum)
  inline RbtCoord Min(const RbtCoord& coord1, const RbtCoord& coord2) {
    return RbtCoord( std::min(coord1.x,coord2.x),
             std::min(coord1.y,coord2.y),
             std::min(coord1.z,coord2.z));
  }

  //Returns maximum of two coords (component-wise maximum)
  inline RbtCoord Max(const RbtCoord& coord1, const RbtCoord& coord2) {
    return RbtCoord( std::max(coord1.x,coord2.x),
             std::max(coord1.y,coord2.y),
             std::max(coord1.z,coord2.z));
  }

  //DM 28 Jul 1999 - returns component-wise min coord for coord list (i.e.min x, min y, min z)
  inline RbtCoord Min(const RbtCoordList& cl) {
    RbtDoubleList cx,cy,cz;
    //Extract arrays of x,y,z components separately
    std::transform(cl.begin(),cl.end(),std::back_inserter(cx),ExtractXCoord);
    std::transform(cl.begin(),cl.end(),std::back_inserter(cy),ExtractYCoord);
    std::transform(cl.begin(),cl.end(),std::back_inserter(cz),ExtractZCoord);
    //Compose the min coord from the min of each of the components
    return RbtCoord(*(std::min_element(cx.begin(),cx.end())),
            *(std::min_element(cy.begin(),cy.end())),
            *(std::min_element(cz.begin(),cz.end())));
  }

  //DM 28 Jul 1999 - returns component-wise max coord for coord list (i.e. max x, max y, max z)
  inline RbtCoord Max(const RbtCoordList& cl) {
    RbtDoubleList cx,cy,cz;
    //Extract arrays of x,y,z components separately
    std::transform(cl.begin(),cl.end(),std::back_inserter(cx),ExtractXCoord);
    std::transform(cl.begin(),cl.end(),std::back_inserter(cy),ExtractYCoord);
    std::transform(cl.begin(),cl.end(),std::back_inserter(cz),ExtractZCoord);
    //Compose the max coord from the max of each of the components
    return RbtCoord(*(std::max_element(cx.begin(),cx.end())),
            *(std::max_element(cy.begin(),cy.end())),
            *(std::max_element(cz.begin(),cz.end())));
  }
  

  //Returns center of "mass" of coords in coord list
  inline RbtCoord GetCenterOfMass(const RbtCoordList& coordList) {
    RbtCoord com;
    com=std::accumulate(coordList.begin(),coordList.end(),com);
    //Check for divide by zero
    return coordList.empty() ? com : com/coordList.size();
  }

}


#endif //_RBTCOORD_H_
