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

#ifdef _RBTOBJECTCOUNTER
#include "RbtDebug.h"

RbtObjectCounter::RbtObjectCounter() :
  nConstr(0),nCopyConstr(0),nDestr(0)
{
}

RbtObjectCounter::~RbtObjectCounter()
{
}

ostream& operator<<(ostream& s, const RbtObjectCounter& counter)
{
  return s << "Constr=" << counter.nConstr
	   << ", Copy constr=" << counter.nCopyConstr
	   << ", Destr=" << counter.nDestr
	   << ", Net total=" << counter.nConstr + counter.nCopyConstr - counter.nDestr;
}

//The one and only RbtObjectCounterMap definition
RbtObjectCounterMap Rbt::theObjectCounterMap;

#endif //_RBTOBJECTCOUNTER


