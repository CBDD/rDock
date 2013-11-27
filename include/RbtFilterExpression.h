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

//Takes a list of tokens and parsers them to create a Filter Expression

#ifndef _RBT_FILTEREXPRESSION_H_
#define _RBT_FILTEREXPRESSION_H_

#include "RbtTypes.h"
#include "RbtGPTypes.h"
#include "RbtTokenIter.h"
#include "RbtFilterExpressionVisitor.h"
#include "RbtError.h"


class RbtFilterExpressionVisitor;

class RbtFilterExpression 
{
public:
  static RbtString _CT;
  virtual void Print()=0;
  virtual void Accept(RbtFilterExpressionVisitor&)=0;
  virtual RbtInt GetNOps()=0;
  virtual SmartPtr<RbtFilterExpression> GetOp(RbtInt)=0;
    ///////////////////
    // Destructor
    //////////////////
  virtual ~RbtFilterExpression();
  RbtReturnType GetValue(){return value;};
  void SetValue(RbtReturnType v){value = v;};
  
    ////////////////////
    // Private methods
    ////////////////////

protected:
  RbtFilterExpression();  // Default constructor disabled
  RbtReturnType value;
};

//Useful typedefs
typedef SmartPtr<RbtFilterExpression> RbtFilterExpressionPtr;  //Smart pointer
typedef vector<RbtFilterExpressionPtr> RbtFilterExpressionList;
                    //Vector of smart pointers
typedef RbtFilterExpressionList::iterator RbtFilterExpressionListIter;
typedef RbtFilterExpressionList::const_iterator 
                                 RbtFilterExpressionListConstIter;

class FilterVbleExp : public RbtFilterExpression
{
public:
    FilterVbleExp(const RbtVble&);
    virtual ~FilterVbleExp();
    void Accept (RbtFilterExpressionVisitor& visitor);
    void Print(){cout << "v" << vble.GetName() << " ";};
    RbtInt GetNOps() { return 0;}
    const RbtVble& GetVble() const { return vble;}
    RbtFilterExpressionPtr GetOp(RbtInt i)
    {
        throw RbtError(_WHERE_, "Vble Expressions don't have operands");
    }
private:
    const RbtVble& vble;
};

class FilterLogExp : public RbtFilterExpression
{
public:
    FilterLogExp(RbtFilterExpressionPtr);
    virtual ~FilterLogExp();
    void Accept (RbtFilterExpressionVisitor& visitor);
    void Print(){cout << "log "; operand->Print();};
    virtual RbtInt GetNOps() { return 1;}
    RbtFilterExpressionPtr GetOp(RbtInt i)
    {
        if (i == 0)
            return operand;
        else
            throw RbtError(_WHERE_, "Log Expressions only have 1 operand");
    }
private:
    RbtFilterExpressionPtr operand;
};

class FilterExpExp : public RbtFilterExpression
{
public:
    FilterExpExp(RbtFilterExpressionPtr);
    virtual ~FilterExpExp();
    void Accept (RbtFilterExpressionVisitor& visitor);
    void Print(){cout << "exp "; operand->Print();};
    virtual RbtInt GetNOps() { return 1;}
    RbtFilterExpressionPtr GetOp(RbtInt i)
    {
        if (i == 0)
            return operand;
        else
            throw RbtError(_WHERE_, "Exp Expressions only have 1 operand");
    }
private:
    RbtFilterExpressionPtr operand;
};

class FilterAddExp : public RbtFilterExpression
{
public:
    FilterAddExp(RbtFilterExpressionPtr, RbtFilterExpressionPtr);
    virtual ~FilterAddExp();
    void Accept (RbtFilterExpressionVisitor& visitor);
    void Print(){cout << "add "; operand1->Print(); operand2->Print();};
    virtual RbtInt GetNOps() { return 2;}
    RbtFilterExpressionPtr GetOp(RbtInt i)
    {
        if (i == 0)
            return operand1;
        if (i == 1)
            return operand2;
        else
            throw RbtError(_WHERE_, "Add Expressions only have 2 operands");
    }
private:
    RbtFilterExpressionPtr operand1, operand2;
};

class FilterSubExp : public RbtFilterExpression
{
public:
    FilterSubExp(RbtFilterExpressionPtr, RbtFilterExpressionPtr);
    virtual ~FilterSubExp();
    void Accept (RbtFilterExpressionVisitor& visitor);
    void Print(){cout << "sub "; operand1->Print(); operand2->Print();};
    virtual RbtInt GetNOps() { return 2;}
    RbtFilterExpressionPtr GetOp(RbtInt i)
    {
        if (i == 0)
            return operand1;
        if (i == 1)
            return operand2;
        else
            throw RbtError(_WHERE_, "Sub Expressions only have 2 operands");
    }
private:
    RbtFilterExpressionPtr operand1, operand2;
};

class FilterMulExp : public RbtFilterExpression
{
public:
    FilterMulExp(RbtFilterExpressionPtr, RbtFilterExpressionPtr);
    virtual ~FilterMulExp();
    void Accept (RbtFilterExpressionVisitor& visitor);
    void Print(){cout << "mul "; operand1->Print(); operand2->Print();};
    virtual RbtInt GetNOps() { return 2;}
    RbtFilterExpressionPtr GetOp(RbtInt i)
    {
        if (i == 0)
            return operand1;
        if (i == 1)
            return operand2;
        else
            throw RbtError(_WHERE_, "Mul Expressions only have 2 operands");
    }
private:
    RbtFilterExpressionPtr operand1, operand2;
};

class FilterDivExp : public RbtFilterExpression
{
public:
    FilterDivExp(RbtFilterExpressionPtr, RbtFilterExpressionPtr);
    virtual ~FilterDivExp();
    void Accept (RbtFilterExpressionVisitor& visitor);
    void Print(){cout << "div "; operand1->Print(); operand2->Print();};
    virtual RbtInt GetNOps() { return 2;}
    RbtFilterExpressionPtr GetOp(RbtInt i)
    {
        if (i == 0)
            return operand1;
        if (i == 1)
            return operand2;
        else
            throw RbtError(_WHERE_, "Div Expressions only have 2 operands");
    }
private:
    RbtFilterExpressionPtr operand1, operand2;
};

class FilterAndExp : public RbtFilterExpression
{
public:
    FilterAndExp(RbtFilterExpressionPtr, RbtFilterExpressionPtr);
    virtual ~FilterAndExp();
    void Accept (RbtFilterExpressionVisitor& visitor);
    void Print(){cout << "and "; operand1->Print(); operand2->Print();};
    virtual RbtInt GetNOps() { return 2;}
    RbtFilterExpressionPtr GetOp(RbtInt i)
    {
        if (i == 0)
            return operand1;
        if (i == 1)
            return operand2;
        else
            throw RbtError(_WHERE_, "And Expressions only have 2 operands");
    }
private:
    RbtFilterExpressionPtr operand1, operand2;
};

class FilterIfExp : public RbtFilterExpression
{
public:
    FilterIfExp(RbtFilterExpressionPtr, RbtFilterExpressionPtr, 
                RbtFilterExpressionPtr);
    virtual ~FilterIfExp();
    void Accept (RbtFilterExpressionVisitor& visitor);
    void Print(){cout << "if "; operand1->Print(); operand2->Print();operand3->Print();};
    virtual RbtInt GetNOps() { return 3;}
    RbtFilterExpressionPtr GetOp(RbtInt i)
    {
        if (i == 0)
            return operand1;
        if (i == 1)
            return operand2;
        if (i == 2)
            return operand3;
        else
            throw RbtError(_WHERE_, "Add Expressions only have 3 operands");
    }
private:
    RbtFilterExpressionPtr operand1, operand2, operand3;
};

#endif //_RbtFilterExpression_H_
