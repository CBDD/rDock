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

//Stores the information that is global to the interpreter
//in this case, values of ctes, names of variables,...
//It basically contains a map of variables, class RbtVble is also defined here

#ifndef _RBT_CONTEXT_H_
#define _RBT_CONTEXT_H_

#include "RbtGPTypes.h"
#include "RbtBaseSF.h"
#include "RbtModel.h"
#include "RbtDockingSite.h"
#include "RbtVble.h"
#include <fstream>
using std::ifstream;

typedef map < RbtString, RbtVblePtr > RbtStringVbleMap; //Map of Vbles
typedef map < RbtInt, RbtVblePtr > RbtIntVbleMap; //Map of Vbles
typedef RbtStringVbleMap::iterator RbtStringVbleMapIter;
typedef RbtIntVbleMap::iterator RbtIntVbleMapIter;



class RbtContext 
{
public:

    static RbtString _CT;
      ///////////////////
      // Constructors
      ///////////////////
    RbtContext(const RbtContext& c);
    RbtContext(); // default constructor disabled
      ///////////////////
      // Destructors
      ///////////////////
    virtual ~RbtContext();  
    virtual void Assign(RbtString, RbtReturnType) = 0;
    virtual void Assign(RbtInt, RbtReturnType) = 0;
    virtual const RbtVble& GetVble(RbtInt)=0;
    virtual const RbtVble& GetVble(RbtString)=0;
    virtual void SetVble(RbtInt key, const RbtVble& v) = 0;
//    virtual RbtString GetName(RbtInt)=0;
//    virtual RbtReturnType GetValue(RbtInt)=0;
//    virtual RbtString GetName(RbtString)=0;
//    virtual RbtReturnType GetValue(RbtString)=0;

};

class RbtCellContext : public RbtContext
{
public:
    //static RbtString _CT;
    RbtCellContext(ifstream& ifile);
    RbtCellContext();
    RbtCellContext(const RbtCellContext& c);
    virtual ~RbtCellContext();  
    void Assign(RbtInt key, RbtReturnType val)
    {
        RbtIntVbleMapIter it = vm.find(key);
        if (it != vm.end())
          vm[key]->SetValue(val);
        else
        {
            ostrstream s;
            s << val << ends;
            vm[key] = new RbtVble(s.str(), val);
        }
    }
    void Assign (RbtString s, RbtReturnType val)
    {
        throw RbtError (_WHERE_, "This is not a string context");
    }
//    RbtString GetName(RbtInt key) { return vm[key].GetName();}
//    RbtReturnType GetValue(RbtInt key) {return vm[key].GetValue();}
//    RbtString GetName(RbtString){return "";}
//    RbtReturnType GetValue(RbtString){return 0.0;}
    const RbtVble& GetVble(RbtInt key){ return *(vm[key]);};
    void SetVble(RbtInt key, const RbtVble& v){ *(vm[key]) = v;};
    const RbtVble& GetVble(RbtString key)
    { 
        throw RbtError (_WHERE_, "This is not a string context");
    }
//    void Clear();


private:
    RbtIntVbleMap vm;
    RbtInt ninputs;
};

class RbtStringContext : public RbtContext
{
    
public:
    RbtStringContext();
    RbtStringContext(SmartPtr<ifstream> ifile);  
    RbtStringContext(const RbtStringContext& c);
    virtual ~RbtStringContext();  
    void Assign(RbtString key, RbtReturnType val)
    {
      RbtStringVbleMapIter it = vm.find(key);
      if (it != vm.end())
        vm[key]->SetValue(val);
      else 
      {
        vm[key] = new RbtVble(key, val);
      }
    }
    void Assign (RbtInt i, RbtReturnType val)
    {
      throw RbtError (_WHERE_, "This is not a cell context");
    }

//    RbtString GetName(RbtInt){return "";}
//    RbtReturnType GetValue(RbtInt){return 0.0;}
//    RbtString GetName(RbtString key) { return vm[key].GetName();}
//    RbtReturnType GetValue(RbtString key) {return vm[key].GetValue();}
    RbtDouble Get(RbtModelPtr lig, RbtString name);
    RbtDouble Get(RbtModelPtr rec, RbtDockingSitePtr site, RbtString name);
    RbtDouble Get(RbtBaseSF* spSF, RbtString name, RbtModelPtr lig);
    const RbtVble& GetVble(RbtString key)
    { 
      return *(vm[key]);
    }
    const RbtVble& GetVble(RbtInt key)
    { 
        throw RbtError (_WHERE_, "This is not a cell context");
    }
    void SetVble(RbtInt key, const RbtVble& v){ *(vm[""]) = v;};
    void UpdateLigs(RbtModelPtr lig);
    void UpdateSite(RbtModelPtr rec, RbtDockingSitePtr site);
    void UpdateScores(RbtBaseSF* spSF, RbtModelPtr lig);
private:
    RbtStringVbleMap vm;

};

//Useful typedefs
typedef SmartPtr<RbtCellContext> RbtCellContextPtr;  //Smart pointer
typedef vector<RbtCellContextPtr> RbtCellContextList;//Vector of smart pointers
typedef RbtCellContextList::iterator RbtCellContextListIter;
typedef RbtCellContextList::const_iterator RbtCellContextListConstIter;
typedef SmartPtr<RbtStringContext> RbtStringContextPtr;  //Smart pointer
typedef vector<RbtStringContextPtr> RbtStringContextList;//Vector of smart pointers
typedef RbtStringContextList::iterator RbtStringContextListIter;
typedef RbtStringContextList::const_iterator RbtStringContextListConstIter;
typedef SmartPtr<RbtContext> RbtContextPtr;  //Smart pointer
typedef vector<RbtContextPtr> RbtContextList;//Vector of smart pointers
typedef RbtContextList::iterator RbtContextListIter;
typedef RbtContextList::const_iterator RbtContextListConstIter;

#endif //_RbtContext_H_
