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
#include <iomanip>
#include <errno.h>

#include "RbtPrincipalAxes.h"
#include "RbtError.h"
#include "RbtPlane.h"

//Template Numerical Toolkit for eigenvalue solver
#include <tnt.h>
#include <jama_eig.h>

//Special case for water
//Allow for symmetry
RbtPrincipalAxes Rbt::GetSolventPrincipalAxes(const RbtAtomPtr& oAtom, const RbtAtomPtr& h1Atom, const RbtAtomPtr& h2Atom) {
  RbtPrincipalAxes retVal;
  RbtCoord oC = oAtom->GetCoords();
  RbtCoord h1C = h1Atom->GetCoords();
  RbtCoord h2C = h2Atom->GetCoords();
  RbtVector v1 = h1C - oC;
  RbtVector v2 = h2C - oC;
  //COM = oxygen
  retVal.com = oC;
  //Axis 1 = Average of O-H bond vectors
  retVal.axis1 = RbtVector(0.5 * (v1 + v2)).Unit();
  //Axis 3 = normal to H-O-H plane
  retVal.axis3 = RbtPlane(oC, h1C, h2C).VNorm();
  //Axis 2 = perpendicular to first two axes
  retVal.axis2 = Rbt::Cross(retVal.axis1,retVal.axis3);
  //Check that H1 is in the positive quadrant of axis1 and axis2
  //Guaranteed to be positive along axis1 so only need to check axis2
  //We do this to ensure canonical axes
  if (v1.Dot(retVal.axis2) < 0.0) {
  	retVal.axis2 = -retVal.axis2;
  	retVal.axis3 = -retVal.axis3;
  }
//  cout << "Solvent: Axis 1 = " << retVal.axis1;
//  cout << ": d1 = " << v1.Dot(retVal.axis1) << "; d2 = " << v2.Dot(retVal.axis1) << endl;
//  cout << "Solvent: Axis 2 = " << retVal.axis2;
//  cout << ": d1 = " << v1.Dot(retVal.axis2) << "; d2 = " << v2.Dot(retVal.axis2) << endl;
//  cout << "Solvent: Axis 3 = " << retVal.axis3;
//  cout << ": d1 = " << v1.Dot(retVal.axis3) << "; d2 = " << v2.Dot(retVal.axis3) << endl;
  return retVal;
}

//Calculates principal axes and center of mass for the atoms in the atom list
//
//
//                             (  Ixx -Ixy -Ixz )
//Moment of Inertia Tensor I = ( -Iyx  Iyy -Iyz )
//                             ( -Izx -Izy  Izz )
//
//             N
//where Ipp = Sum( m(k)*(r(k)*r(k)-p(k)*p(k)) )
//            k=1
//
//                   N
//      Ipq = Iqp = Sum( m(k)*p(k)*q(k) )
//     (p<>q)       k=1
//
//where m(k) is the mass of the k-th atom, r(k) is the vector {x(k),y(k),z(k)}
//from the atom to the center of mass, and p,q = x,y,z.
//
//Principal axes are eigenvectors of I, principal moments are eigenvalues of I
//
RbtPrincipalAxes Rbt::GetPrincipalAxes(const RbtAtomList& atomList)
{
  const RbtUInt N = 3;//Array size

  RbtPrincipalAxes principalAxes;//Return parameter
  if (atomList.empty()) {
    return principalAxes;
  }
  //Special case for water
  else if (atomList.size() == 3) {
    Rbt::isAtomicNo_eq isOxygen(8);
    Rbt::isAtomicNo_eq isHydrogen(1);
    if (isOxygen(atomList[0]) && isHydrogen(atomList[1]) && isHydrogen(atomList[2])) {
      return Rbt::GetSolventPrincipalAxes(atomList[0],atomList[1],atomList[2]);
    }
  }
  principalAxes.com = Rbt::GetCenterOfMass(atomList);//Store center of mass


  //Construct the moment of inertia tensor
  TNT::Array2D<RbtDouble> inertiaTensor(N,N,0.0);
  for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    RbtVector r = (*iter)->GetCoords() - principalAxes.com;//Vector from center of mass to atom
    RbtDouble m = (*iter)->GetAtomicMass();//Atomic mass
    RbtDouble rx2 = r.x*r.x;
    RbtDouble ry2 = r.y*r.y;
    RbtDouble rz2 = r.z*r.z;
    //Diagonal elements (moments of inertia)
    RbtDouble dIxx = m*(ry2+rz2);//=r^2 - x^2
    RbtDouble dIyy = m*(rx2+rz2);//=r^2 - y^2
    RbtDouble dIzz = m*(rx2+ry2);//=r^2 - z^2
    inertiaTensor[0][0] += dIxx;
    inertiaTensor[1][1] += dIyy;
    inertiaTensor[2][2] += dIzz;

    //Off-diagonal elements (products of inertia) - symmetric matrix
    RbtDouble dIxy = m*r.x*r.y;
    RbtDouble dIxz = m*r.x*r.z;
    RbtDouble dIyz = m*r.y*r.z;
    inertiaTensor[0][1] -= dIxy;
    inertiaTensor[1][0] -= dIxy;
    inertiaTensor[0][2] -= dIxz;
    inertiaTensor[2][0] -= dIxz;
    inertiaTensor[1][2] -= dIyz;
    inertiaTensor[2][1] -= dIyz;
  }

  JAMA::Eigenvalue<RbtDouble> eigenSolver(inertiaTensor);
  TNT::Array1D<RbtDouble> eigenValues(N);
  TNT::Array2D<RbtDouble> eigenVectors(N,N);
  eigenSolver.getRealEigenvalues(eigenValues);
  eigenSolver.getV(eigenVectors);

  //Load the principal axes and moments into the return parameter
  //We need to sort these so that axis1 is the first principal axis, axis2 the second, axis3 the third.
  //With only three elements to sort, this is probably as good a way as any:
  RbtUInt idx1=0;
  RbtUInt idx2=1;
  RbtUInt idx3=2;
  if (eigenValues[idx1] > eigenValues[idx2]) std::swap(idx1,idx2);
  if (eigenValues[idx1] > eigenValues[idx3]) std::swap(idx1,idx3);
  if (eigenValues[idx2] > eigenValues[idx3]) std::swap(idx2,idx3);
  principalAxes.axis1 = RbtVector(eigenVectors[0][idx1],eigenVectors[1][idx1],eigenVectors[2][idx1]);
  principalAxes.axis2 = RbtVector(eigenVectors[0][idx2],eigenVectors[1][idx2],eigenVectors[2][idx2]);
  principalAxes.axis3 = RbtVector(eigenVectors[0][idx3],eigenVectors[1][idx3],eigenVectors[2][idx3]);
  principalAxes.moment1 = eigenValues[idx1];
  principalAxes.moment2 = eigenValues[idx2];
  principalAxes.moment3 = eigenValues[idx3];

  //DM 28 Jun 2001 - for GA crossovers in particular we need to ensure the principal axes returned
  //are always consistent for a given ligand conformation. Currently the direction of the axes vectors
  //is not controlled as for e.g. (+1,0,0) is degenerate with (-1,0,0).
  //Method is to arbitrarily check the first atom and invert the axes vectors if necessary
  //to ensure the atom coords lie in the positive quadrant of the coordinate space of principal axes 1 and 2
  //Principal axis 3 is then defined to give a right-handed set of axes
  //
  //LIMITATION: If atom 1 lies exactly on PA#1 or PA#2 this check will fail. Ideally we would like to test an atom
  //on the periphery of the molecule.
  RbtCoord c0 = (atomList.front())->GetCoords() - principalAxes.com;
  RbtDouble d1 = c0.Dot(principalAxes.axis1);
  RbtDouble d2 = c0.Dot(principalAxes.axis2);
  //RbtDouble d3 = c0.Dot(principalAxes.axis3);
  //cout << "Before: d1,d2,d3=" << d1 << "\t" << d2 << "\t" << d3 << endl;
  if (d1 < 0.0)
    principalAxes.axis1 = -principalAxes.axis1;
  if (d2 < 0.0)
    principalAxes.axis2 = -principalAxes.axis2;
  principalAxes.axis3 = Rbt::Cross(principalAxes.axis1,principalAxes.axis2);
  //d1 = c0.Dot(principalAxes.axis1);
  //d2 = c0.Dot(principalAxes.axis2);
  //d3 = c0.Dot(principalAxes.axis3);
  //cout << "After: d1,d2,d3=" << d1 << "\t" << d2 << "\t" << d3 << endl;

  //cout << "(TNT) Axis1=" << principalAxes.axis1 << "; Moment1=" << principalAxes.moment1 << endl;
  //cout << "(TNT) Axis2=" << principalAxes.axis2 << "; Moment2=" << principalAxes.moment2 << endl;
  //cout << "(TNT) Axis3=" << principalAxes.axis3 << "; Moment3=" << principalAxes.moment3 << endl;
  return principalAxes;
}


//Calculates principal axes and center of mass for the coords in the coord list (assumes all masses=1)
RbtPrincipalAxes Rbt::GetPrincipalAxes(const RbtCoordList& coordList)
{
  const RbtUInt N = 3;//Array size

  RbtPrincipalAxes principalAxes;//Return parameter
  principalAxes.com = Rbt::GetCenterOfMass(coordList);//Store center of mass

  //Construct the moment of inertia tensor
  TNT::Array2D<RbtDouble> inertiaTensor(N,N,0.0);
  for (RbtCoordListConstIter iter = coordList.begin(); iter != coordList.end(); iter++) {
    RbtVector r = (*iter) - principalAxes.com;//Vector from center of mass to coord
    RbtDouble rx2 = r.x*r.x;
    RbtDouble ry2 = r.y*r.y;
    RbtDouble rz2 = r.z*r.z;
    //Diagonal elements (moments of inertia)
    RbtDouble dIxx = ry2+rz2;//=r^2 - x^2
    RbtDouble dIyy = rx2+rz2;//=r^2 - y^2
    RbtDouble dIzz = rx2+ry2;//=r^2 - z^2
    inertiaTensor[0][0] += dIxx;
    inertiaTensor[1][1] += dIyy;
    inertiaTensor[2][2] += dIzz;

    //Off-diagonal elements (products of inertia) - symmetric matrix
    RbtDouble dIxy = r.x*r.y;
    RbtDouble dIxz = r.x*r.z;
    RbtDouble dIyz = r.y*r.z;
    inertiaTensor[0][1] -= dIxy;
    inertiaTensor[1][0] -= dIxy;
    inertiaTensor[0][2] -= dIxz;
    inertiaTensor[2][0] -= dIxz;
    inertiaTensor[1][2] -= dIyz;
    inertiaTensor[2][1] -= dIyz;

  }

  JAMA::Eigenvalue<RbtDouble> eigenSolver(inertiaTensor);
  TNT::Array1D<RbtDouble> eigenValues(N);
  TNT::Array2D<RbtDouble> eigenVectors(N,N);
  eigenSolver.getRealEigenvalues(eigenValues);
  eigenSolver.getV(eigenVectors);

  //Load the principal axes and moments into the return parameter
  //We need to sort these so that axis1 is the first principal axis, axis2 the second, axis3 the third.
  //With only three elements to sort, this is probably as good a way as any:
  RbtUInt idx1=0;
  RbtUInt idx2=1;
  RbtUInt idx3=2;
  if (eigenValues[idx1] > eigenValues[idx2]) std::swap(idx1,idx2);
  if (eigenValues[idx1] > eigenValues[idx3]) std::swap(idx1,idx3);
  if (eigenValues[idx2] > eigenValues[idx3]) std::swap(idx2,idx3);
  principalAxes.axis1 = RbtVector(eigenVectors[0][idx1],eigenVectors[1][idx1],eigenVectors[2][idx1]);
  principalAxes.axis2 = RbtVector(eigenVectors[0][idx2],eigenVectors[1][idx2],eigenVectors[2][idx2]);
  principalAxes.axis3 = RbtVector(eigenVectors[0][idx3],eigenVectors[1][idx3],eigenVectors[2][idx3]);
  principalAxes.moment1 = eigenValues[idx1];
  principalAxes.moment2 = eigenValues[idx2];
  principalAxes.moment3 = eigenValues[idx3];

  return principalAxes;
}

//Aligns the principal axes of the atoms in the atom list to lie along refAxes
//If required (bAlignCOM=true), also aligns the center of mass with refAxes.com
//
//LIMITATION: does not check for non-orthogonal alignAxes
RbtQuat Rbt::AlignPrincipalAxes(RbtAtomList& atomList, const RbtPrincipalAxes& refAxes, RbtBool bAlignCOM) throw (RbtError)
{
  RbtPrincipalAxes prAxes = Rbt::GetPrincipalAxes(atomList);
  RbtQuat q = Rbt::GetQuatFromAlignAxes(prAxes,refAxes);

  std::for_each(atomList.begin(),atomList.end(),Rbt::TranslateAtom(-prAxes.com));
  std::for_each(atomList.begin(),atomList.end(),Rbt::RotateAtomUsingQuat(q));
  if (bAlignCOM) {
  	std::for_each(atomList.begin(),atomList.end(),Rbt::TranslateAtom(refAxes.com));
  }
  else {
 	std::for_each(atomList.begin(),atomList.end(),Rbt::TranslateAtom(prAxes.com));
  }
  return q;
}

RbtQuat Rbt::GetQuatFromAlignAxes(const RbtPrincipalAxes& prAxes,
                                    const RbtPrincipalAxes& refAxes) throw (RbtError) {
  //1) Determine the quaternion needed to align axis1 with reference
  RbtQuat q1 = Rbt::GetQuatFromAlignVectors(prAxes.axis1, refAxes.axis1);
  //2) Apply the transformation to axis2 
  RbtVector axis2 = q1.Rotate(prAxes.axis2);
  //3) Determine the quaternion needed to align axis2 with reference
  RbtQuat q2 = Rbt::GetQuatFromAlignVectors(axis2, refAxes.axis2);
  //Return the quaternion product (equivalent to both transformations combined)
  return q2 * q1;
}

RbtQuat Rbt::GetQuatFromAlignVectors(const RbtVector& v, const RbtVector& ref) throw (RbtError) {
  RbtQuat retVal;
  //Unitise the two vectors
  RbtDouble len = Rbt::Length(v);
  RbtDouble refLen = Rbt::Length(ref);
  if ( (len < 0.001) || (refLen < 0.001) ) {
    throw RbtBadArgument(_WHERE_,"zero length vector (v or ref)");
  }
  RbtVector vUnit = v/len;
  RbtVector refUnit = ref/refLen;
  //Determine the rotation axis and angle needed to overlay the two vectors
  RbtVector axis = Rbt::Cross(vUnit,refUnit);
  //DM 15 March 2006: check for zero-length rotation axis
  //This indicates the vectors are already aligned
  RbtDouble axisLen = Rbt::Length(axis);
  if (axisLen > 0.001) {  
    RbtDouble cosPhi = Rbt::Dot(vUnit,refUnit);
    if (cosPhi < -1.0) {
      cosPhi = -1.0;
    }
    else if (cosPhi > 1.0) {
      cosPhi = 1.0;
    }
    errno = 0;
    //Convert rotation axis and angle to a quaternion
    RbtDouble halfPhi = 0.5* acos(cosPhi);
    if ( (halfPhi > 0.001) && (errno != EDOM) ) {
      RbtVector axisUnit = axis/axisLen;
      retVal = RbtQuat(cos(halfPhi),sin(halfPhi)*axisUnit);
    }
  }
  return retVal;
}
