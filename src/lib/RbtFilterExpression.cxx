/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtFilterExpression.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtFilterExpression.h"

RbtString RbtFilterExpression::_CT("FilterExpression");

RbtFilterExpression::RbtFilterExpression()
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtFilterExpression::~RbtFilterExpression()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

FilterVbleExp::FilterVbleExp(const RbtVble& vn) : vble(vn)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterVbleExp::~FilterVbleExp()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void FilterVbleExp::Accept(RbtFilterExpressionVisitor& visitor)
{
    visitor.VisitVbleExp(this);
}

FilterLogExp::FilterLogExp(RbtFilterExpressionPtr op1)
            : operand(op1)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterLogExp::~FilterLogExp()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void FilterLogExp::Accept(RbtFilterExpressionVisitor& visitor)
{
    visitor.VisitLogExp(this);
}


FilterExpExp::FilterExpExp(RbtFilterExpressionPtr op1)
            : operand(op1)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterExpExp::~FilterExpExp()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void FilterExpExp::Accept(RbtFilterExpressionVisitor& visitor)
{
    visitor.VisitExpExp(this);
}

FilterAddExp::FilterAddExp(RbtFilterExpressionPtr op1, 
                           RbtFilterExpressionPtr op2)
            : operand1(op1), operand2(op2)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterAddExp::~FilterAddExp()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void FilterAddExp::Accept(RbtFilterExpressionVisitor& visitor)
{
    visitor.VisitAddExp(this);
}


FilterSubExp::FilterSubExp(RbtFilterExpressionPtr op1, 
                           RbtFilterExpressionPtr op2)
            : operand1(op1), operand2(op2)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterSubExp::~FilterSubExp()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void FilterSubExp::Accept(RbtFilterExpressionVisitor& visitor)
{
    visitor.VisitSubExp(this);
}


FilterMulExp::FilterMulExp(RbtFilterExpressionPtr op1, 
                           RbtFilterExpressionPtr op2)
            : operand1(op1), operand2(op2)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterMulExp::~FilterMulExp()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void FilterMulExp::Accept(RbtFilterExpressionVisitor& visitor)
{
    visitor.VisitMulExp(this);
}


FilterDivExp::FilterDivExp(RbtFilterExpressionPtr op1, 
                           RbtFilterExpressionPtr op2)
            : operand1(op1), operand2(op2)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterDivExp::~FilterDivExp()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void FilterDivExp::Accept(RbtFilterExpressionVisitor& visitor)
{
    visitor.VisitDivExp(this);
}

FilterAndExp::FilterAndExp(RbtFilterExpressionPtr op1, 
                           RbtFilterExpressionPtr op2)
            : operand1(op1), operand2(op2)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterAndExp::~FilterAndExp()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void FilterAndExp::Accept(RbtFilterExpressionVisitor& visitor)
{
    visitor.VisitAndExp(this);
}


FilterIfExp::FilterIfExp(RbtFilterExpressionPtr op1, 
                         RbtFilterExpressionPtr op2, 
                         RbtFilterExpressionPtr op3)
            : operand1(op1), operand2(op2), operand3(op3)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterIfExp::~FilterIfExp()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void FilterIfExp::Accept(RbtFilterExpressionVisitor& visitor)
{
    visitor.VisitIfExp(this);
}


