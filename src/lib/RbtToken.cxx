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

#include "RbtToken.h"
#include "RbtCommands.h"
#include "RbtDebug.h"

RbtString RbtToken::_CT("RbtToken");

    ///////////////////
    // Constructors
    ///////////////////
RbtToken::RbtToken(const RbtVble& v, RbtCommands c, RbtBool isvble) :
    comm(c),
    vble(v),
    isvble(isvble)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtToken::RbtToken(const RbtVble& v) :
    RbtToken(v, -1, true)
{}

RbtToken::RbtToken(RbtCommands c) :
    RbtToken(RbtVble(), c, false) 
{}

RbtToken::RbtToken(const RbtToken& t) :
    RbtToken(t.vble, t.comm, t.isvble)
{}
 
    ///////////////////
    // Destructor
    //////////////////
RbtToken::~RbtToken()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

const RbtVble& RbtToken::GetVble() const
{
    if (!isvble) 
        throw RbtError(_WHERE_, "the token is not a vble");
    return vble;
}

RbtBool RbtToken::IsVble()
{
    return (isvble);
}

RbtBool RbtToken::IsLog()
{
    return (!isvble && comm.IsLog());
}

RbtBool RbtToken::IsExp()
{
    return (!isvble && comm.IsExp());
}

RbtBool RbtToken::IsAdd()
{
    return (!isvble && comm.IsAdd());
}

RbtBool RbtToken::IsSub()
{
    return (!isvble && comm.IsSub());
}

RbtBool RbtToken::IsMul()
{
    return (!isvble && comm.IsMul());
}

RbtBool RbtToken::IsDiv()
{
    return (!isvble && comm.IsDiv());
}

RbtBool RbtToken::IsAnd()
{
    return (!isvble && comm.IsAnd());
}

RbtBool RbtToken::IsIf()
{
    return (!isvble && comm.IsIf());
}
