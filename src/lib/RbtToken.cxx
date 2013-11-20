/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtToken.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtToken.h"
#include "RbtCommands.h"
#include "RbtDebug.h"

RbtString RbtToken::_CT("RbtToken");

    ///////////////////
    // Constructors
    ///////////////////
RbtToken::RbtToken(const RbtVble& v) : isvble(true), vble(v), comm(-1)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtToken::RbtToken(RbtCommands c) : isvble(false), comm(c), vble(RbtVble()) 
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtToken::RbtToken(const RbtToken& t) : isvble(t.isvble), comm(t.comm),
                                        vble(t.vble)
{
    _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}
 
    ///////////////////
    // Destructor
    //////////////////
RbtToken::~RbtToken()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
};
  
/*void RbtToken::copy(const RbtToken &t)
{
    isvble = t.isvble;
    comm = t.comm;
    vble = t.vble;
}*/

const RbtVble& RbtToken::GetVble() const
{
    if (!isvble) 
        throw RbtError(_WHERE_, "the token is not a vble");
    return vble;
}
/* void RbtToken::SetVbleNumber(RbtInt vn)
{
    vble.SetValue(vn);
}*/

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

