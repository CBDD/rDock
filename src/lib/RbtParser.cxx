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

#include "RbtParser.h"
#include "RbtDebug.h"

RbtString RbtParser::_CT("RbtParser");
RbtParser::RbtParser(const RbtParser& p)
{
    _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtParser::RbtParser()
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtParser::~RbtParser()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}
RbtFilterExpressionPtr RbtParser::Parse(RbtTokenIterPtr ti, RbtContextPtr contextp)
{
    RbtTokenPtr t = ti->Current();
    ti->Next(contextp);
    if (t->IsVble())
	{
        return new FilterVbleExp(t->GetVble());
	}
    RbtFilterExpressionPtr op1 = Parse(ti, contextp);
    if (t->IsLog())
        return new FilterLogExp(op1);
    if (t->IsExp())
        return new FilterExpExp(op1);
    RbtFilterExpressionPtr op2 = Parse(ti, contextp);
    if (t->IsAdd())
        return new FilterAddExp(op1,op2);
    if (t->IsSub())
        return new FilterSubExp(op1,op2);
    if (t->IsMul())
        return new FilterMulExp(op1,op2);
    if (t->IsDiv())
        return new FilterDivExp(op1,op2);
    if (t->IsAnd())
        return new FilterAndExp(op1,op2);
    if (t->IsIf())
    {
        RbtFilterExpressionPtr op3 = Parse(ti, contextp);
        return new FilterIfExp(op1, op2, op3);
    }
	delete ti;
    throw RbtError(_WHERE_, "Could not parse " + t->GetVble().GetName());
}
