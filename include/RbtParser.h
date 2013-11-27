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

//Parser class. Takes a list of tokens and parsers them to create a
//Filter Expression

#ifndef _RBT_PARSER_H_
#define _RBT_PARSER_H_

#include "RbtTypes.h"
#include "RbtTokenIter.h"
#include "RbtToken.h"
#include "RbtContext.h"
#include "RbtFilterExpression.h"

class RbtFilterExpression;

class RbtParser 
{
public:
  static RbtString _CT;
    ///////////////////
    // Constructors
    ///////////////////
  RbtParser(const RbtParser&);
  RbtParser();
  RbtFilterExpressionPtr Parse(RbtTokenIterPtr, RbtContextPtr);
 
    ///////////////////
    // Destructor
    //////////////////
  virtual ~RbtParser();
  
    ////////////////////
    // Private methods
    ////////////////////

private:
  
};

//Useful typedefs
typedef SmartPtr<RbtParser> RbtParserPtr;  //Smart pointer
typedef vector<RbtParserPtr> RbtParserList;//Vector of smart pointers
typedef RbtParserList::iterator RbtParserListIter;
typedef RbtParserList::const_iterator RbtParserListConstIter;

  // Compare class to compare different genomes inside a population
  // This is used to sort the genomes in a population depending
  // of the value of their scoring function. 

#endif //_RbtParser_H_
