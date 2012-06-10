/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtCommands.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Defines the interface for a token to be used by the parser

#ifndef _RBT_COMMANDS_H_
#define _RBT_COMMANDS_H_

#include "RbtTypes.h"


class RbtCommands 
{
public:
  enum RbtCommName {ADD=0, SUB=1, MUL=2, DIV=3, IF=4, LOG=5,
                    EXP= 6, RCTE=7, AND=8, NULLC=-1};

    ///////////////////
    // Constructors
    ///////////////////
  RbtCommands(const RbtCommands& c): name(c.name), nargs(c.nargs){};
  inline RbtCommands(RbtInt ncomm)
  {
      name = RbtCommName(ncomm);
      switch (name)
      {
          case ADD:
          case SUB:
          case MUL:
          case DIV:
          case AND:
              nargs = 2;
              break;
          case IF:
              nargs = 3;
              break;
          case LOG:
          case EXP:
              nargs = 1;
              break;
          case RCTE:
              nargs = 0;
              break;
          case NULLC:
              nargs = 0;
              break;
      } 
  };

  inline RbtCommName GetName() const
  {
      return name;
  };

  inline RbtInt GetNArgs() const
  {
      return nargs;
  };

  inline RbtBool IsAdd()const {return (name == ADD);};
  inline RbtBool IsSub()const {return (name == SUB);};
  inline RbtBool IsMul()const {return (name == MUL);};
  inline RbtBool IsDiv()const {return (name == DIV);};
  inline RbtBool IsIf()const {return (name == IF);};
  inline RbtBool IsLog()const {return (name == LOG);};
  inline RbtBool IsExp()const {return (name == EXP);};
  inline RbtBool IsRCte()const {return (name == RCTE);};
  inline RbtBool IsAnd()const {return (name == AND);};

 
    ///////////////////
    // Destructor
    //////////////////
  virtual ~RbtCommands(){};
  
  virtual void copy(const RbtCommands &c) { name = c.name; nargs = c.nargs;};

    ////////////////////
    // Private methods
    ////////////////////

private:
  RbtCommands();  // Default constructor disabled
  RbtCommName name;
  RbtInt nargs;
  
};

//Useful typedefs
//typedef SmartPtr<RbtCommands> RbtCommandsPtr;  //Smart pointer
//typedef vector<RbtCommandsPtr> RbtCommandsList;//Vector of smart pointers
//typedef RbtCommandsList::iterator RbtCommandsListIter;
//typedef RbtCommandsList::const_iterator RbtCommandsListConstIter;

  // Compare class to compare different genomes inside a population
  // This is used to sort the genomes in a population depending
  // of the value of their scoring function. 

#endif //_RbtCommands_H_
