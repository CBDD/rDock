/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtParser.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
