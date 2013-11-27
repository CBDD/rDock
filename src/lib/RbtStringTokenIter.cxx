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
