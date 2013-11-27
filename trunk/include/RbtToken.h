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

//Token class. Defines a token of the grammar to be used by the parser.

#ifndef _RBT_TOKEN_H_
#define _RBT_TOKEN_H_

#include "RbtConfig.h"
#include "RbtCommands.h"
#include "RbtVble.h"



class RbtToken 
{
public:
  static RbtString _CT;
    ///////////////////
    // Constructors
    ///////////////////

  RbtToken(const RbtVble&);
  RbtToken(RbtCommands);
  RbtToken(const RbtToken&);
 
    ///////////////////
    // Destructor
    //////////////////
  virtual ~RbtToken();
  const RbtVble& GetVble() const;
//  void SetVbleNumber(RbtInt);
  
//  void copy(const RbtToken &);
  RbtBool IsVble();
  RbtBool IsLog();
  RbtBool IsExp();
  RbtBool IsAdd();
  RbtBool IsSub();
  RbtBool IsMul();
  RbtBool IsDiv();
  RbtBool IsAnd();
  RbtBool IsIf();
  
//  RbtBool operator==(const RbtToken&) const;
//  ostream& Print(ostream&) const;
//  friend ostream& operator<<(ostream& s, const RbtToken &p);

    ////////////////////
    // Private methods
    ////////////////////

private:
  RbtToken();  // Default constructor disabled
  RbtCommands comm;
  const RbtVble& vble;
  RbtBool isvble;
  
};

//Useful typedefs
typedef SmartPtr<RbtToken> RbtTokenPtr;  //Smart pointer
typedef vector<RbtTokenPtr> RbtTokenList;//Vector of smart pointers
typedef RbtTokenList::iterator RbtTokenListIter;
typedef RbtTokenList::const_iterator RbtTokenListConstIter;

  // Compare class to compare different genomes inside a population
  // This is used to sort the genomes in a population depending
  // of the value of their scoring function. 

#endif //_RbtToken_H_
