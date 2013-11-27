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

//Base class for pharmacophore constraints

#ifndef _RBT_CONSTRAINT_H_
#define _RBT_CONSTRAINT_H_

#include "RbtTypes.h"
#include "RbtCoord.h"
#include "RbtAtom.h"
#include "RbtModel.h"
#include "RbtError.h"

class RbtConstraint 
{
public:
  static RbtString _CT;
    ///////////////////
    // Constructors
    ///////////////////
  RbtConstraint(RbtCoord, RbtDouble);
  RbtConstraint(const RbtConstraint&);
  RbtConstraint();  
 
    ///////////////////
    // Destructor
    //////////////////
  virtual ~RbtConstraint();
  
  void copy(const RbtConstraint &);
  ostream& Print(ostream&) const;
  friend ostream& operator<<(ostream&, const RbtConstraint &);

  RbtCoord GetCoords() const { return coord; };
  RbtDouble GetTolerance() const { return tolerance; };
  virtual void AddAtomList(RbtModelPtr, RbtBool bCheck=true) = 0;
  RbtDouble Score() const;

    ////////////////////
    // Private methods
    ////////////////////
protected:
  RbtCoord coord;
  RbtDouble tolerance;
  RbtAtomList m_atomList;
  
private:
    /////////////////////
    // Private data
    /////////////////////  
};

//Useful typedefs
typedef SmartPtr<RbtConstraint> RbtConstraintPtr;  //Smart pointer
typedef vector<RbtConstraintPtr> RbtConstraintList;//Vector of smart pointers
typedef RbtConstraintList::iterator RbtConstraintListIter;
typedef RbtConstraintList::const_iterator RbtConstraintListConstIter;

namespace Rbt {
  RbtConstraintPtr CreateConstraint(RbtCoord& c, RbtDouble& t, RbtString& n, RbtBool bCount=true);
  void ZeroCounters();
  void ReadConstraint(istream& ifile, RbtConstraintPtr& cnt, RbtBool bCount=true);
  void ReadConstraintFromMoe(istream& ifile, RbtConstraintPtr& cnt, RbtBool bCount=true);
  void ReadStartMoe(istream& ifile);
  void ReadConstraints(istream& ii, RbtConstraintList& cl, RbtBool bCount=true);
}

//7 Feb 2005 (DM, Enspiral Discovery)
//New constraint type - any heavy atom
//Likely to be inefficient for large ligands, but has the advantage
//of being simple to code and run
class RbtHeavyConstraint : public RbtConstraint
{
    public:
  RbtHeavyConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t) {};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtHBAConstraint : public RbtConstraint
{
    public:
  RbtHBAConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t) {};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtHBDConstraint : public RbtConstraint
{
    public:
  RbtHBDConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtHydroConstraint : public RbtConstraint
{
    public:
  RbtHydroConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtHydroAliphaticConstraint : public RbtConstraint
{
    public:
  RbtHydroAliphaticConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtHydroAromaticConstraint : public RbtConstraint
{
    public:
  RbtHydroAromaticConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};


class RbtNegChargeConstraint : public RbtConstraint
{
    public:
  RbtNegChargeConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtPosChargeConstraint : public RbtConstraint
{
    public:
  RbtPosChargeConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtRingAromaticConstraint : public RbtConstraint
{
    public:
  RbtRingAromaticConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

#endif //_RbtConstraint_H_
