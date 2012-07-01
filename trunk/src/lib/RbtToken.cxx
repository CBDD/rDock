/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

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

