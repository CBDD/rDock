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

#include "RbtDebug.h"
#include "RbtObserver.h"

////////////////////////////////////////
//Constructors/destructors
RbtObserver::RbtObserver() {
  _RBTOBJECTCOUNTER_CONSTR_("RbtObserver");
}

RbtObserver::~RbtObserver() {
  _RBTOBJECTCOUNTER_DESTR_("RbtObserver");
}

