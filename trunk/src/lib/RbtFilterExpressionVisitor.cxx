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

#include "RbtFilterExpressionVisitor.h"
#include <cmath>
#include <cerrno>
#include <limits.h>

void PrintVisitor::VisitVbleExp(FilterVbleExp* fe)
{
	  RbtString name = fe->GetVble().GetName();
		if (name[0] == 'c') // constant
			cout << fe->GetVble().GetValue() << " ";
		else
      cout << fe->GetVble().GetName() << " ";
}

void PrettyPrintVisitor::VisitVbleExp(FilterVbleExp* fe)
{
//	  cout << "Vble ";
    cout << fe->GetVble().GetName() << " ";
}

void EvaluateVisitor::VisitVbleExp(FilterVbleExp* fe)
{
    fe->SetValue(fe->GetVble().GetValue());
}


void PrintVisitor::VisitAddExp(FilterAddExp* fe)
{
    cout << "+ ";
    for (RbtInt i = 0 ; i < fe->GetNOps() ; i++)
        fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitAddExp(FilterAddExp* fe)
{
    cout << "(";
    fe->GetOp(0)->Accept(*this);
    cout << "+ ";
    fe->GetOp(1)->Accept(*this);
    cout << ")";
}

void EvaluateVisitor::VisitAddExp(FilterAddExp* fe)
{
    fe->GetOp(0)->Accept(*this);
    fe->GetOp(1)->Accept(*this);
    fe->SetValue(fe->GetOp(0)->GetValue() + fe->GetOp(1)->GetValue());
//		cout << "AddExp: " << fe->GetOp(0)->GetValue() << "\t" << fe->GetOp(1)->GetValue() << endl;
}

void PrintVisitor::VisitSubExp(FilterSubExp* fe)
{
    cout << "- ";
    for (RbtInt i = 0 ; i < fe->GetNOps() ; i++)
        fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitSubExp(FilterSubExp* fe)
{
    cout << "(";
    fe->GetOp(0)->Accept(*this);
    cout << "- ";
    fe->GetOp(1)->Accept(*this);
    cout << ")";
}

void EvaluateVisitor::VisitSubExp(FilterSubExp* fe)
{
    fe->GetOp(0)->Accept(*this);
    fe->GetOp(1)->Accept(*this);
    fe->SetValue(fe->GetOp(0)->GetValue() - fe->GetOp(1)->GetValue());
//	cout << "sub " << fe->GetOp(0)->GetValue() <<"\t" << fe->GetOp(1)->GetValue() << endl;
}

void PrintVisitor::VisitMulExp(FilterMulExp* fe)
{
    cout << "* ";
    for (RbtInt i = 0 ; i < fe->GetNOps() ; i++)
        fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitMulExp(FilterMulExp* fe)
{
    cout << "(";
    fe->GetOp(0)->Accept(*this);
    cout << "* ";
    fe->GetOp(1)->Accept(*this);
    cout << ")";
}

void EvaluateVisitor::VisitMulExp(FilterMulExp* fe)
{
    fe->GetOp(0)->Accept(*this);
    fe->GetOp(1)->Accept(*this);
    fe->SetValue(fe->GetOp(0)->GetValue() * fe->GetOp(1)->GetValue());
//	cout << "Mul " << fe->GetOp(0)->GetValue() <<"\t" << fe->GetOp(1)->GetValue() << endl;
}

void PrintVisitor::VisitDivExp(FilterDivExp* fe)
{
    cout << "/ ";
    for (RbtInt i = 0 ; i < fe->GetNOps() ; i++)
        fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitDivExp(FilterDivExp* fe)
{
    cout << "(";
    fe->GetOp(0)->Accept(*this);
    cout << "div ";
    fe->GetOp(1)->Accept(*this);
    cout << ")";
}

void EvaluateVisitor::VisitDivExp(FilterDivExp* fe)
{
    fe->GetOp(1)->Accept(*this);
    fe->GetOp(0)->Accept(*this);
    RbtReturnType v0 = fe->GetOp(0)->GetValue();
    RbtReturnType v1 = fe->GetOp(1)->GetValue();
//	cout << "div " << v0 <<"\t" << v1 << endl;
    if (fabs(v1) < 0.000001)
        fe->SetValue(v0);
    else
        fe->SetValue(v0 / v1);
}

void PrintVisitor::VisitAndExp(FilterAndExp* fe)
{
    cout << "and ";
    for (RbtInt i = 0 ; i < fe->GetNOps() ; i++)
        fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitAndExp(FilterAndExp* fe)
{
    cout << "( (";
    fe->GetOp(0)->Accept(*this);
    cout << " > 0) and (";
    fe->GetOp(1)->Accept(*this);
    cout << " > 0) )";
}

void EvaluateVisitor::VisitAndExp(FilterAndExp* fe)
{
    fe->GetOp(1)->Accept(*this);
    RbtReturnType v0 = fe->GetOp(0)->GetValue();
		if (v0 <= 0.0)
			fe->SetValue(0.0);	//false
		else
		{
      fe->GetOp(0)->Accept(*this);
      RbtReturnType v1 = fe->GetOp(1)->GetValue();
			if (v1 <= 0.0)
				fe->SetValue(0.0); //false
			else
				fe->SetValue(1.0); //true
		}
}

void PrintVisitor::VisitLogExp(FilterLogExp* fe)
{
    cout << "log ";
    for (RbtInt i = 0 ; i < fe->GetNOps() ; i++)
        fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitLogExp(FilterLogExp* fe)
{
    cout << "(log ";
    fe->GetOp(0)->Accept(*this);
    cout << ")";
}

void EvaluateVisitor::VisitLogExp(FilterLogExp* fe)
{
    fe->GetOp(0)->Accept(*this);
    RbtReturnType v0 = fe->GetOp(0)->GetValue();
    if (fabs(v0) < 0.000001)
        fe->SetValue(0.0);
    else
        fe->SetValue(log(fabs(v0)));
//	cout << "log " << v0 << endl;
}

void PrintVisitor::VisitExpExp(FilterExpExp* fe)
{
    cout << "exp ";
    for (RbtInt i = 0 ; i < fe->GetNOps() ; i++)
        fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitExpExp(FilterExpExp* fe)
{
    cout << "(exp ";
    fe->GetOp(0)->Accept(*this);
    cout << ")";
}

void EvaluateVisitor::VisitExpExp(FilterExpExp* fe)
{
    fe->GetOp(0)->Accept(*this);
    RbtReturnType v0 = fe->GetOp(0)->GetValue();
		errno = 0;
    if (v0 > 20)
			
		{
//			float m = FLOAT_MAX; //numeric_limits<float>::max();
//			cout << "greater " << m << "\n";
        fe->SetValue(exp(20));}
    else if (v0 < -200)
        fe->SetValue(0.0);
    else
        fe->SetValue(exp(v0));
//	cout << "exp " << v0 <<  "\t" << fe->GetValue() <<endl;
}

void PrintVisitor::VisitIfExp(FilterIfExp* fe)
{
    cout << "if ";
    for (RbtInt i = 0 ; i < fe->GetNOps() ; i++)
        fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitIfExp(FilterIfExp* fe)
{
    cout << "\nif ";
    fe->GetOp(0)->Accept(*this);
    cout << " > 0 then\n\t";
    fe->GetOp(1)->Accept(*this);
    cout << "\nelse\n\t";
    fe->GetOp(2)->Accept(*this);
    cout << "\nend\n";
}

void EvaluateVisitor::VisitIfExp(FilterIfExp* fe)
{
    fe->GetOp(0)->Accept(*this);
    RbtReturnType v0 = fe->GetOp(0)->GetValue();
    if (v0 > 0.0)
    {
        fe->GetOp(1)->Accept(*this);
        fe->SetValue(fe->GetOp(1)->GetValue());
    }
    else
    {
        fe->GetOp(2)->Accept(*this);
        fe->SetValue(fe->GetOp(2)->GetValue());
    }
}

RbtFilterExpressionVisitor::~RbtFilterExpressionVisitor()
{
}

RbtFilterExpressionVisitor::RbtFilterExpressionVisitor()
{
}

PrintVisitor::PrintVisitor(RbtContextPtr co) : contextp(co)
{
}

PrettyPrintVisitor::PrettyPrintVisitor(RbtContextPtr co) : contextp(co)
{
}

EvaluateVisitor::EvaluateVisitor(RbtContextPtr co) : contextp(co)
{
}
