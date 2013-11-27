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

//Cell Class. Cell represents a function or input variable of the CGP

#ifndef _RBT_CELL_H_
#define _RBT_CELL_H_

#include "RbtGPTypes.h"

class RbtCell 
{
public:
    ///////////////////
    // Constructors
    ///////////////////
    RbtCell(){evaluated = false; name = "";constant = false;};
    RbtCell(const RbtCell& c) 
        {evaluated = c.evaluated; name = c.name; constant = c.constant;
         result = c.result;};

    ///////////////////
    // Destructor
    //////////////////
    virtual ~RbtCell(){};
 
    RbtBool Evaluated()const {return evaluated;};
    RbtBool Named()const {return (name != "");};
    RbtBool Constant()const {return constant ;};
    void Clear()
    {
        if (!constant)
        {
            evaluated = false; 
            name = "";
        }
    };
    RbtReturnType GetResult()const {return result;};
    RbtString GetName()const {return name;};
    void SetName(RbtString s)
        {name = s; };
    void SetResult(RbtReturnType r)
        {result = r; evaluated = true;};
    void SetConstant(RbtReturnType r)
        {result = r; evaluated = true; constant = true;
        ostrstream nm; nm << r << ends; name = nm.str();};
    void ResetConstant() 
        { evaluated = false; name = ""; 
            if (constant) constant = false; else cout <<"ERROR\n";};

private:
    RbtBool evaluated, constant;
    RbtString name;
    RbtReturnType result;
//    friend ostream& operator<<(ostream& s, const RbtCell &p);
};

//Useful typedefs
typedef SmartPtr<RbtCell> RbtCellPtr;  //Smart pointer
typedef vector<RbtCellPtr> RbtCellList;//Vector of smart pointers
typedef RbtCellList::iterator RbtCellListIter;
typedef RbtCellList::const_iterator RbtCellListConstIter;

#endif //_RbtCell_H_
