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

#ifndef _RBTMACROS_H_
#define _RBTMACROS_H_

/*
    definition for __FALLTHROUGH__ macro.
    if std=c++17 or higher use new attribute notation [[fallthrough]]
    else go for compiler specific options
        if gcc/g++ and std=c++14 or lower use __attribute__ ((fallthrough))
    to be defined for other compilers, like clang
*/

#if __cplusplus > 201402L
    #define __FALLTHROUGH__ [[fallthrough]]
#else
    #ifdef __GNUG__
        #define __FALLTHROUGH__ __attribute__ ((fallthrough))
    #endif
#endif

#endif  // _RBTMACROS_H_
