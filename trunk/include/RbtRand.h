/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtRand.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
