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

//Wrapper around Randint class
//Function provided to return reference to single instance (singleton) of
//RbtRand

#ifndef _RBTRAND_H_
#define _RBTRAND_H_

#include "RandInt.h"
#include "RbtTypes.h"
#include "RbtCoord.h"

class RbtRand
{
  /////////////
  //Constructor
 public:
  RbtRand();
  /////////////
  //Destructor
  ~RbtRand();

  ////////////////
  //Public methods

  //Seed the random number generator
  void Seed(RbtInt seed=0);
  //Seed the random number generator from the system clock
  void SeedFromClock();
  //Returns current seed
  RbtInt GetSeed();
  //Get a random double between 0 and 1 (inlined)
  RbtDouble GetRandom01() {return m_rand.fdraw();};
  //Get a random integer between 0 and nMax-1
  RbtInt GetRandomInt(RbtInt nMax);
  //Get a random unit vector distributed evenly over the surface of a sphere
  RbtVector GetRandomUnitVector();
  RbtDouble GetGaussianRandom(RbtDouble, RbtDouble);
  RbtDouble GetCauchyRandom(RbtDouble, RbtDouble);

 private:
  Randint m_rand;//Random number generator
};

///////////////////////////////////////
//Non-member functions in Rbt namespace

namespace Rbt
{
  //Returns reference to single instance of RbtRand class (singleton)
  RbtRand& GetRbtRand();
}
#endif //_RBTRAND_H_
