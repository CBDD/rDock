/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/GP/RbtCommand.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

//Command Class. Command represents a function or input variable of the CGP

#ifndef _RBT_COMMAND_H_
#define _RBT_COMMAND_H_

#include "RbtGPTypes.h"

#include <cmath>

class RbtCommand 
{
public:
    virtual RbtReturnType Execute() = 0;
    virtual int GetNArgs() = 0;
    virtual RbtString GetName() = 0;
    virtual void SetArg(int i, RbtReturnType f) = 0;// { *(Arg[i]) = f;};
    virtual void SetNameArg(int i, RbtString& n) = 0;// { NameArg[i] = n;};
    static void Clear() {ntabs = 0;};
 
    ///////////////////
    // Destructor
    //////////////////
    virtual ~RbtCommand(){};
    static RbtInt ntabs; 
    static RbtBool inside; 

protected:
    ///////////////////
    // Constructors
    ///////////////////
    RbtCommand(){};  // Default constructor disabled
    RbtString tabs(RbtInt n)
    {
        RbtString st = "";
        for (RbtInt i = 0 ; i < n ; i++)
            st += "  ";
        return st;
    }
private:
//    friend ostream& operator<<(ostream& s, const RbtCommand &p);
};

//Useful typedefs
typedef SmartPtr<RbtCommand> RbtCommandPtr;  //Smart pointer
typedef vector<RbtCommandPtr> RbtCommandList;//Vector of smart pointers
typedef RbtCommandList::iterator RbtCommandListIter;
typedef RbtCommandList::const_iterator RbtCommandListConstIter;
class AddCommand : public RbtCommand
{
public:
    int GetNArgs() { return 2; };
    RbtString GetName() 
    { return ("(" + NameArg[0] + "+ " + NameArg[1] + ")");};
    AddCommand(){}; //{Arg = RbtReturnType(2);};
    void SetArg(int i, RbtReturnType f) { Arg[i] = f;};
    void SetNameArg(int i, RbtString& n) { NameArg[i] = n;};
    RbtReturnType Execute() 
    {
        return Arg[0] + Arg[1];
    };
private:
    RbtReturnType Arg[2];
    RbtString NameArg[2];
};

class SubCommand : public RbtCommand
{
public:
    int GetNArgs() { return 2; };
    RbtString GetName() 
    { return ("(" + NameArg[0] + "- " + NameArg[1]+")");};
    void SetArg(int i, RbtReturnType f) { Arg[i] = f;};
    void SetNameArg(int i, RbtString& n) { NameArg[i] = n;};
    SubCommand(){};
    RbtReturnType Execute() 
    {
        return Arg[0] - Arg[1];
    };
private:
    RbtReturnType Arg[2];
    RbtString NameArg[2];
};

class MulCommand : public RbtCommand
{
public:
    int GetNArgs() { return 2; };
    RbtString GetName() 
    { return ("(" + NameArg[0] + "* " + NameArg[1] + ")");};
    void SetArg(int i, RbtReturnType f) { Arg[i] = f;};
    void SetNameArg(int i, RbtString& n) { NameArg[i] = n;};
    MulCommand(){};
    RbtReturnType Execute() 
    {
        return Arg[0] * Arg[1];
    };
private:
    RbtReturnType Arg[2];
    RbtString NameArg[2];
};

class DivCommand : public RbtCommand
{
public:
    int GetNArgs() { return 2; };
    RbtString GetName() 
    { return ("(" + NameArg[0] + "div " + NameArg[1]+ ")");};
    DivCommand(){};
    void SetArg(int i, RbtReturnType f) { Arg[i] = f;};
    void SetNameArg(int i, RbtString& n) { NameArg[i] = n;};
    RbtReturnType Execute() 
    {
        if (fabs(Arg[1]) < 0.000001)
            return Arg[0];
        return Arg[0] / Arg[1];
    };
private:
    RbtReturnType Arg[2];
    RbtString NameArg[2];
};

class IfCommand : public RbtCommand
{
public:
    int GetNArgs() { return 3; };
    RbtString GetName() 
    {
	//	return ("\niff(" + NameArg[0] + "," + NameArg[1] + "," + 
	//			NameArg[2] + ")");};
        if (inside)
            return ("iff(" + NameArg[0] + " , " + NameArg[1] + " , " +
                     NameArg[2] + ")");
        return ("\n" + tabs(ntabs - 1) + "if " + NameArg[0] + "> 0 then\n" +
                tabs(ntabs) + NameArg[1] + "\n" + tabs(ntabs - 1)
                + "else \n" + tabs(ntabs) + NameArg[2]+"\n" + 
                tabs(ntabs - 1) + "end\n" );}; 
    IfCommand(){};
    void SetArg(int i, RbtReturnType f) { Arg[i] = f;};
    void SetNameArg(int i, RbtString& n) {NameArg[i] = n;};
    RbtReturnType Execute() 
    {
        if (Arg[0] > 0.0)
            return Arg[1];
        return Arg[2];
    };
private:
    RbtReturnType Arg[3];
    RbtString NameArg[3];
};

class LogCommand : public RbtCommand
{
public:
    int GetNArgs() { return 1; };
    RbtString GetName() 
    { return ("log(" + NameArg[0] + ")");};
    LogCommand(){};
    void SetArg(int i, RbtReturnType f) { Arg[i] = f;};
    void SetNameArg(int i, RbtString& n) { NameArg[i] = n;};
    RbtReturnType Execute() 
    {
        if (fabs(Arg[0]) < 0.000001)
            return 0;
        return (log(fabs(Arg[0])));
    };
private:
    RbtReturnType Arg[1];
    RbtString NameArg[1];
};

class ExpCommand : public RbtCommand
{
public:
    int GetNArgs() { return 1; };
    RbtString GetName() 
    { return ("exp(" + NameArg[0] + ")");};
    ExpCommand(){};
    void SetArg(int i, RbtReturnType f) { Arg[i] = f;};
    void SetNameArg(int i, RbtString& n) { NameArg[i] = n;};
    RbtReturnType Execute() 
    {
        if (Arg[0] > 200)
            return exp(200);
		if (Arg[0] < -200)
			return 0;
        return exp(Arg[0]);
    };
private:
    RbtReturnType Arg[2];
    RbtString NameArg[2];
};
#endif //_RbtCommand_H_
