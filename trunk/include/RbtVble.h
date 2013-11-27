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

//Vble class. Stores information about the variables used by the filters.

#ifndef _RBT_VBLE_H_
#define _RBT_VBLE_H_

#include "RbtTypes.h"
#include "RbtGPTypes.h"
#include <strstream>

class RbtVble
{
public:
  static RbtString _CT;
  enum VbleType {CTE, LIG, SCORE, SITE};
  RbtVble() : value(0.0), vt(CTE), name("")
  {
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
  }; // default const.
  RbtVble(RbtString s, RbtReturnType val = 0.0) : value(val), name(s)
  {
    if (s.find("LIG_") == 0)
      vt = LIG;
    else if (s.find("SITE_") == 0)
      vt = SITE;
    else if (s.find("SCORE") == 0)
      vt = SCORE;
    else
      vt = CTE;
   _RBTOBJECTCOUNTER_CONSTR_(_CT);
  }
  RbtVble(const RbtVble& v) : value(v.value), vt(v.vt), name(v.name)
  {
   _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
  }
  void SetValue(RbtReturnType val)
  {
    value = val;
		if (name == "")
		{
			ostrstream s;
			s << value << ends;
			name = s.str();
		}
  }
  RbtReturnType GetValue() const
  {
    return value;
  }
  virtual ~RbtVble()
  {
   _RBTOBJECTCOUNTER_DESTR_(_CT);
  }
  void SetName(RbtString nm)
  {
    name = nm;
  }
  RbtString GetName() const
  {
    return name;
  }

  RbtBool IsLig() {return (vt == LIG);}
  RbtBool IsScore() {return (vt == SCORE);}
  RbtBool IsSite() {return (vt == SITE);}

private:
  VbleType vt;
  RbtReturnType value;
  RbtString name;
};

//Useful typedefs
typedef SmartPtr<RbtVble> RbtVblePtr;  //Smart pointer
typedef vector<RbtVblePtr> RbtVbleList;//Vector of smart pointers
typedef RbtVbleList::iterator RbtVbleListIter;
typedef RbtVbleList::const_iterator RbtVbleListConstIter;

#endif //_RbtVble_H_
