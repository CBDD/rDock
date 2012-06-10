/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtStringTokenIter.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtStringTokenIter.h"
#include "RbtCommands.h"
#include <cstdlib>
#include <errno.h>

RbtString RbtStringTokenIter::_CT("RbtStringTokenIter");

    ///////////////////
    // Constructors
    ///////////////////
RbtStringTokenIter::RbtStringTokenIter(const istreamPtr fn, RbtContextPtr co) 
                    : filep(fn), contextp(co)
{
  (istream&)(*filep) >> strtok;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtStringTokenIter::RbtStringTokenIter(const RbtStringTokenIter& ti) 
  : current(ti.current), filep(ti.filep), 
    strtok(ti.strtok), contextp(ti.contextp)
{
    _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtStringTokenIter::~RbtStringTokenIter()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}
  
void RbtStringTokenIter::Next(RbtContextPtr)
{
  (*filep) >> strtok;
}

RbtTokenPtr RbtStringTokenIter::Current() 
{ 
  current = translate (strtok);
  return current;
}

void RbtStringTokenIter::copy(const RbtStringTokenIter & ti) 
{
  strtok = ti.strtok;
  current = ti.current;
  filep = ti.filep;
  contextp = ti.contextp;
}

RbtTokenPtr RbtStringTokenIter::translate(RbtString s)
{
  if (s.length() == 0)
    throw RbtError(_WHERE_, "Missing token, can't translate the expression");
  if (s == "+")
    return new RbtToken(RbtCommands(RbtCommands::ADD)); 
  if (s == "-")
    return new RbtToken(RbtCommands(RbtCommands::SUB));
  if (s == "*")
    return new RbtToken(RbtCommands(RbtCommands::MUL));
  if (s == "/")
    return new RbtToken(RbtCommands(RbtCommands::DIV));
  if (s == "if")
  {
    return new RbtToken(RbtCommands(RbtCommands::IF));
  }
  if (s == "log")
    return new RbtToken(RbtCommands(RbtCommands::LOG));
  if (s == "exp")
    return new RbtToken(RbtCommands(RbtCommands::EXP));
  if (s == "and")
    return new RbtToken(RbtCommands(RbtCommands::AND));
  if ((string::npos != s.find("SCORE")) ||
    (string::npos != s.find("SITE")) ||
    (string::npos != s.find("LIG")))
  {
		contextp->Assign(s, 0.0);
    return new RbtToken(contextp->GetVble(s));
  }
  else // I assume is a double. If it is not, need to send an RbtError
  {
    char *error;
    errno = 0;
    RbtDouble val = strtod(s.c_str(), &error);
    if (!errno && !*error)  // This checks for errors 
    {
      contextp->Assign(s, val);
      return new RbtToken(contextp->GetVble(s));
    }
  }
  throw RbtError(_WHERE_, "Can't read " + s);
}
