/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtParser.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

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
