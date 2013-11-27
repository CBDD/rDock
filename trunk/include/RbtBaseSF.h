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

//Base class for all scoring functions
#ifndef _RBTBASESF_H_
#define _RBTBASESF_H_

#include "RbtConfig.h"
#include "RbtBaseObject.h"

class RbtSFAgg;//forward declaration

class RbtBaseSF : public RbtBaseObject
{
 public:
  
  //Class type string
  static RbtString _CT;
  //Parameter names
  static RbtString _WEIGHT;
  //DM 09 Apr 2002 - max distance range over which scoring function operates
  //Subclasses are free not to use if not required
  static RbtString _RANGE;
  //DM 17 Jan 2006 - hardcoded name of scoring function branch under which to save
  //"system" scoring function terms (e.g. intra-receptor, intra-solvent).
  //This is a dirty (and hopefully temporary) solution!
  static RbtString _SYSTEM_SF;
  static RbtString _INTRA_SF;//Ditto, for ligand intramolecular terms
  
  ////////////////////////////////////////
  //Constructors/destructors
  virtual ~RbtBaseSF();
  
  //Give aggregates access to RbtBaseSF private methods/data
  friend class RbtSFAgg;
	
  ////////////////////////////////////////
  //Public methods
  ////////////////
  //Fully qualified name, prefixed by all ancestors (e.g. SCORE.INTER.HBOND)
  RbtString GetFullName() const;
  RbtDouble GetWeight() const;
  void SetWeight(RbtDouble);

  RbtDouble GetRange() const;
  void SetRange(RbtDouble);
  
  //Main public method - returns current weighted score
  RbtDouble Score() const;
  //Returns all child component scores as a string-variant map
  //Key = fully qualified component name, value = weighted score
  //(for saving in a Model's data fields)
  virtual void ScoreMap(RbtStringVariantMap& scoreMap) const;
      
  //Aggregate handling methods
  virtual void Add(RbtBaseSF*) throw (RbtError);
  virtual void Remove(RbtBaseSF*) throw (RbtError);
  virtual RbtBool isAgg() const;
  virtual RbtUInt GetNumSF() const;
  virtual RbtBaseSF* GetSF(RbtUInt iSF) const throw (RbtError);
  void Orphan();//Force removal from the parent aggregate
  RbtBaseSF* GetParentSF() const;
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  RbtBaseSF(const RbtString& strClass, const RbtString& strName);
  RbtBaseSF();
  //PURE VIRTUAL - DERIVED CLASSES MUST OVERRIDE
  virtual RbtDouble RawScore() const = 0;
  //DM 25 Oct 2000 - track changes to parameter values in local data members
  //ParameterUpdated is invoked by RbtParamHandler::SetParameter
  void ParameterUpdated(const RbtString& strName);
  //Helper method for ScoreMap
  void AddToParentMapEntry(RbtStringVariantMap& scoreMap, RbtDouble rs) const;
  
 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtBaseSF(const RbtBaseSF&);//Copy constructor disabled by default      
  RbtBaseSF& operator=(const RbtBaseSF&);//Copy assignment disabled by default
  
 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////
  
  
 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtBaseSF* m_parent;
  RbtDouble m_weight;
  RbtDouble m_range;
};

//Useful typedefs
typedef vector<RbtBaseSF*> RbtBaseSFList;//Vector of smart pointers
typedef RbtBaseSFList::iterator RbtBaseSFListIter;
typedef RbtBaseSFList::const_iterator RbtBaseSFListConstIter;

#endif //_RBTBASESF_H_
