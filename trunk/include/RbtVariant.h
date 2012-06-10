/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtVariant.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Variant class along the lines of the Visual Basic variant data type.
//Handles int, double, const char*, string and vector<string> types
//Underlying value is stored as double and vector<string>

#ifndef _RBTVARIANT_H_
#define _RBTVARIANT_H_

#include <sstream>
using std::istringstream;
using std::ostringstream;
using std::ends;
using std::cout;

#include "RbtTypes.h"
#include "RbtContainers.h"
#include "RbtCoord.h"

//Strings representing TRUE and FALSE
const RbtString _TRUE("TRUE");
const RbtString _FALSE("FALSE");

class RbtVariant
{
 public:	
  ////////////////////////////////////////
  //Constructors/destructors
  RbtVariant() : m_d(0.0) {}
  RbtVariant(RbtInt i) {SetDouble(i);}
  RbtVariant(RbtDouble d) {SetDouble(d);}
  RbtVariant(const RbtString& s) {SetString(s);}
  RbtVariant(const char* c) {SetString(RbtString(c));}
  RbtVariant(const RbtStringList& sl) {SetStringList(sl);}
  RbtVariant(RbtBool b) {SetBool(b);}
  RbtVariant(const RbtCoord& c) {SetCoord(c);}
  //Renders a vector of doubles into comma-separated strings of maxCols in length
  RbtVariant(const RbtDoubleList& dl, RbtInt maxCols, RbtInt precision) {
    SetDoubleList(dl, maxCols, precision);
  }
  virtual ~RbtVariant() {m_sl.clear();}

  ///////////////////////////////////////////////
  //Stream functions
  ///////////////////////////////////////////////

  //Insertion operator
  friend ostream& operator<<(ostream& s, const RbtVariant& v) {
	//Try and guess the most appropriate output format
	//Empty variant
	if (v.isEmpty())
		s << "Undefined variant";
	//String/double
	else if (v.Size() == 1)
		s << v.m_sl.front();
	//String list
	else {
		s << endl;
		for (RbtStringListConstIter iter = v.m_sl.begin(); iter != v.m_sl.end(); iter++)
			s << *iter << endl;
	}
	return s;
  }
  
  ///////////////////////////////////////////////
  //Operator functions:
  /////////////////////

  //Set methods
  RbtVariant& operator=(RbtInt i) {
	SetDouble(i);
	return *this;
  }
  RbtVariant& operator=(RbtDouble d) {
	SetDouble(d);
	return *this;
  }
  RbtVariant& operator=(const RbtString& s) {
	SetString(s);
	return *this;
  }
  RbtVariant& operator=(const char* c) {
	SetString(RbtString(c));
	return *this;
  }
  RbtVariant& operator=(const RbtStringList& sl) {
	SetStringList(sl);
	return *this;
  }
  RbtVariant& operator=(RbtBool b) {
	SetBool(b);
	return *this;
  }  
  RbtVariant& operator=(const RbtCoord& c) {
	SetCoord(c);
	return *this;
  }

  void operator+=(const RbtVariant& v) {
    std::copy(v.m_sl.begin(),v.m_sl.end(),std::back_inserter(m_sl));
    m_d = atof(String().c_str());
  }

  //Conversion operators to convert back to basic types
  operator RbtInt() const {return int(Double());}
  operator RbtDouble() const {return Double();}
  operator RbtString() const {return String();}
  operator RbtStringList() const {return StringList();}
  operator RbtBool() const {return Bool();}
  operator RbtCoord() const {return Coord();}

  ////////////////////////////////////////
  //Public methods
  ////////////////

  //Get methods
  RbtDouble Double() const {return m_d;}
  RbtString String() const {
	return m_sl.empty() ? RbtString() : m_sl.front();
  }
  RbtStringList StringList() const {return m_sl;}
  RbtBool Bool() const {
	return m_d != 0.0 || String() == _TRUE;
  }
  RbtCoord Coord() const {
    RbtCoord c;
    if (m_sl.empty())
      return c;
    else {
      ///istrstream istr(m_sl.front().c_str());
      istringstream istr(m_sl.front().c_str());
      istr >> c;
      return c;
    }
  }      

  RbtUInt Size() const {return m_sl.size();}
  RbtBool isEmpty() const {return m_sl.empty();}

 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  void SetDouble(RbtDouble d) {
	m_d = d;
	m_sl.clear();
	ostringstream ostr;
	//Don't need "ends" with ostringstream apparently
	//(was introducing a non-ASCII \0 char into log files
	//ostr << d << ends;
	ostr << d;
	RbtString s(ostr.str());
	//delete ostr.str();
	m_sl.push_back(s);
  }

  void SetString(const RbtString& s) {
	m_sl.clear();
	if (!s.empty())
		m_sl.push_back(s);
	m_d = atof(String().c_str());
  }

  void SetStringList(const RbtStringList& sl) {
	m_sl = sl;
	m_d = atof(String().c_str()); 
  }

  void SetBool(RbtBool b) {
	if (b) {
		SetString(_TRUE);
	}
	else {
		SetString(_FALSE);
	}
  }

  void SetCoord(const RbtCoord& c) {
    m_d = c.Length();
    m_sl.clear();
    ostringstream ostr;
    //Don't need "ends" with ostringstream apparently
    //(was introducing a non-ASCII \0 char into log files
    //ostr << c << ends;
    ostr << c;
    RbtString s(ostr.str());
    //delete ostr.str();
    m_sl.push_back(s);
  }
  
  void SetDoubleList(const RbtDoubleList& dl, RbtInt maxColumns, RbtInt precision) {
    RbtInt nValues = dl.size();
    m_d = (nValues > 0) ? dl[0] : 0.0;
    m_sl.clear();
    ostringstream ostr;
    ostr.precision(precision);
    ostr.setf(ios_base::fixed,ios_base::floatfield);
    RbtInt lastIndex = nValues - 1;
    for (RbtInt i = 0; i < nValues; ++i) {
        ostr << dl[i];
        if ( (ostr.tellp() >= maxColumns) || (i == lastIndex) ) {
            m_sl.push_back(ostr.str());
            ostr.str("");
        }
        else {
            ostr << ",";
        }
    }
  }
  
 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////

 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtDouble m_d;
  RbtStringList m_sl;
};

//Useful typedefs
//typedef SharedPtr<RbtVariant> RbtVariantPtr;//Smart pointer
typedef vector<RbtVariant> RbtVariantList;//Vector of variants
typedef RbtVariantList::iterator RbtVariantListIter;
typedef RbtVariantList::const_iterator RbtVariantListConstIter;

//Map of (key=string, value=variant)
typedef map<RbtString,RbtVariant> RbtStringVariantMap;
typedef RbtStringVariantMap::iterator RbtStringVariantMapIter;
typedef RbtStringVariantMap::const_iterator RbtStringVariantMapConstIter;

#endif //_RBTVARIANT_H_
