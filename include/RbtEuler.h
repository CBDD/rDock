/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtEuler.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Simple class to manage Euler angles
//Main purpose is to provides conversion methods
//to/from quaternion representation
//
//We adopt a modified "NASA Standard Aeroplane" convention
//in which the y and z axes are swapped relative to the definition in
//http://www.euclideanspace.com/maths/geometry/rotations/euler/index.htm
//(and related links).
//This should be equivalent to the VRML convention.
// 
//An object's orientation is defined in terms of:
//
//1) heading (-PI to +PI) (rotation around z) (y in NASA Aeroplane)
//2) attitude (-PI/2 to +PI/2) (rotation around y) (z in NASA Aeroplane)
//3) bank (-PI to +PI) (rotation around x)
//
//See http://www.euclideanspace.com/maths/geometry/rotations/euler/index.htm
//and related links for definitions and conversion equations.
//
//NOTE: This class is designed to provide fast conversion between Euler angles
//and quaternions. It is the callers responsibility to ensure that:
//a) heading, attitude and bank angles are within range
//b) quaternions passed to the constructor and SetQuat are of unit length
//Run-time exceptions may occur if these conditions are not met.
#ifndef RBTEULER_H_
#define RBTEULER_H_

#include "RbtQuat.h"

class RbtEuler {
    public:
        //Constructor accepting values for heading, attitude and bank (radians).
        //NOTE: For performance reasons, no checks are made that the values
        //are within range
        RbtEuler(RbtDouble heading = 0.0,
                    RbtDouble attitude = 0.0,
                    RbtDouble bank = 0.0)
            : m_heading(heading), m_attitude(attitude), m_bank(bank) {}    
        //Constructor accepting a quaternion.
        //Quaternion is converted to Euler angle representation
        //NOTE: For performance reasons, q is assumed to be of unit length
        //and no further checks are made.
        RbtEuler(const RbtQuat& q) {FromQuat(q);}
        
        //Returns true if all three Euler angles are within their
        //standard ranges
        RbtBool isStandardised() const {
            return (m_heading >= -M_PI) && (m_heading <= M_PI)
            && (m_attitude >= -M_PI/2.0) && (m_attitude <= M_PI/2.0)
            && (m_bank >= -M_PI) && (m_bank <= M_PI);
        }
        //Standardises the Euler angles by conversion to/from quaternion
        void Standardise() {if (!isStandardised()) FromQuat(ToQuat());}
        
        //Gets the current heading angle (radians)
        RbtDouble GetHeading() const {return m_heading;}
        //Gets the current attitude angle (radians)
        RbtDouble GetAttitude() const {return m_attitude;}
        //Gets the current bank angle (radians)
        RbtDouble GetBank() const {return m_bank;}
        //Gets the equivalent quaternion representation
        //NOTE: For performance reasons, no checks are made that the current
        //Euler angles are within range prior to the conversion.
        RbtQuat ToQuat() const;
        //Updates the Euler angles from a quaternion representation
        //NOTE: For performance reasons, q is assumed to be of unit length
        //and no further checks are made.
        void FromQuat(const RbtQuat& q);
        //Convenience method to rotate orientation by a quaternion
        void Rotate(const RbtQuat& q) {FromQuat(q * ToQuat());}
        void Rotate(const RbtVector& axis, RbtDouble theta) {
            Rotate(RbtQuat(axis,theta));
        }
        
        //Insertion operator
        friend ostream& operator<<(ostream& s, const RbtEuler& euler) {
            return s << "Heading = " << euler.m_heading
                    << ", Attitude = " << euler.m_attitude
                    << ", Bank = " << euler.m_bank;
        }
        
        
    private:
        RbtDouble m_heading;//heading angle (radians)
        RbtDouble m_attitude;//attitude angle (radians)
        RbtDouble m_bank;//bank angle (radians)
};
#endif /*RBTEULER_H_*/
