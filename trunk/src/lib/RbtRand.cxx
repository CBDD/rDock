/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

#include <time.h> //Time functions for initialising the random number generator from the system clock

#include "RbtRand.h"
#include "Singleton.h"

/////////////
//Constructor
RbtRand::RbtRand()
{
  //Seed the random number generator
  //Fixed seed in debug mode
  //Seed from system clock in release mode
#ifdef _DEBUG
  Seed();
#else //_DEBUG
  SeedFromClock();
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_("RbtRand");
}

/////////////
//Destructor
RbtRand::~RbtRand()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtRand");
}


////////////////
//Public methods

//Seed the random number generator
void RbtRand::Seed(RbtInt seed)
{
  m_rand.seed(seed);
}

//Seed the random number generator from the system clock
void RbtRand::SeedFromClock()
{
  m_rand.seed(::time(NULL));
}

//Returns current seed
RbtInt RbtRand::GetSeed()
{
  return m_rand.GetSeed();
}

//Get a random integer between 0 and nMax-1
RbtInt RbtRand::GetRandomInt(RbtInt nMax)
{
  RbtInt r = nMax*m_rand.fdraw();
  return (r==nMax) ? nMax-1 : r;
}

//Get a random unit vector distributed evenly over the surface of a sphere
RbtVector RbtRand::GetRandomUnitVector()
{
  RbtDouble z = 2.0 * GetRandom01() - 1.0;
  RbtDouble t = 2.0 * M_PI * GetRandom01();
  RbtDouble w = sqrt( 1.0 - z*z );
  RbtDouble x = w * cos( t );
  RbtDouble y = w * sin( t );
  return RbtVector(x,y,z);
}

//Get a random number from the Normal distribution (mean, variance)
RbtDouble RbtRand::GetGaussianRandom(RbtDouble mean, RbtDouble variance)
{
  for(;;) 
  {
    RbtDouble u1 = GetRandom01();
    RbtDouble u2 = GetRandom01();
    RbtDouble v1 = 2 * u1 - 1;
    RbtDouble v2 = 2 * u2 - 1;
    RbtDouble w = (v1 * v1) + (v2 * v2);
    if (w <= 1) 
    {
      RbtDouble y = sqrt( (-2 * log(w)) / w);
      RbtDouble x1 = v1 * y;
      return(x1 * sqrt(variance) + mean);
    }
  }
}

//Get a random number from the Cauchy distribution (mean, variance)
RbtDouble RbtRand::GetCauchyRandom(RbtDouble mean, RbtDouble variance)
{
  RbtDouble v1 = GetGaussianRandom(mean, variance);
  RbtDouble v2 = GetGaussianRandom(mean, variance);
  if (fabs(v2) > 0.001)
    return (v1 / v2);
  else
    return 0.0;
}
  

///////////////////////////////////////
//Non-member functions in Rbt namespace

//Returns reference to single instance of RbtRand class (singleton)
RbtRand& Rbt::GetRbtRand()
{
  return Singleton<RbtRand>::instance();
}
