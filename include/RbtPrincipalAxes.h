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

//Principal axes calculation routines (in Rbt namespace)

#ifndef _RBTPRINCIPALAXES_H_
#define _RBTPRINCIPALAXES_H_

#include "RbtConfig.h"
#include "RbtAtom.h"

//Struct for holding principal axes info
//Center of mass, three principal axis vectors, three principal moments
class RbtPrincipalAxes
{
 public:
  //Default is X,Y,Z cartesian axes centered at origin
  RbtPrincipalAxes() :
    com(0.0,0.0,0.0),
    axis1(1.0,0.0,0.0),axis2(0.0,1.0,0.0),axis3(0.0,0.0,1.0),
    moment1(1.0),moment2(1.0),moment3(1.0) {}
  RbtCoord com;
  RbtVector axis1;
  RbtVector axis2;
  RbtVector axis3;
  RbtDouble moment1;
  RbtDouble moment2;
  RbtDouble moment3;
};

typedef vector<RbtPrincipalAxes> RbtPrincipalAxesList;
typedef RbtPrincipalAxesList::iterator RbtPrincipalAxesListIter;
typedef RbtPrincipalAxesList::const_iterator RbtPrincipalAxesListConstIter;

namespace Rbt
{
  //Calculates principal axes and center of mass for the atoms in the atom list
  RbtPrincipalAxes GetPrincipalAxes(const RbtAtomList& atomList);
  //Calculates principal axes and center of mass for the coords in the coord list (assumes all masses=1)
  RbtPrincipalAxes GetPrincipalAxes(const RbtCoordList& coordList);
  //Special case for water
  RbtPrincipalAxes GetSolventPrincipalAxes(const RbtAtomPtr& oAtom, const RbtAtomPtr& h1Atom, const RbtAtomPtr& h2Atom);
  //DM 17 Jul 2001 - returns the quaternion used to effect the transformation
  RbtQuat AlignPrincipalAxes(RbtAtomList& atomList, const RbtPrincipalAxes& alignAxes=RbtPrincipalAxes(), RbtBool bAlignCOM=true) throw (RbtError);
  //Returns the quaternion required to align principal axes with reference axes
  RbtQuat GetQuatFromAlignAxes(const RbtPrincipalAxes& prAxes, const RbtPrincipalAxes& refAxes) throw (RbtError);
  //Returns the quaternion required to effect an alignment of v onto ref vector.
  //Vectors do not have to be unit length.
  //RbtBadArgument exception is thrown if either v or ref is zero length.
  RbtQuat GetQuatFromAlignVectors(const RbtVector& v, const RbtVector& ref) throw (RbtError);
}


#endif //_RBTPRINCIPALAXES_H_
