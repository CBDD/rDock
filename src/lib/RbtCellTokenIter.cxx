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

#include "RbtCellTokenIter.h"

RbtString RbtCellTokenIter::_CT("RbtCellTokenIter");

  ///////////////////
  // Constructors
  ///////////////////
RbtCellTokenIter::RbtCellTokenIter(const RbtGPChromosomePtr c,
                                   RbtContextPtr contextp) : chrom(c)
{
  RbtInt startCell = chrom->GetStartingCell();
  if (chrom->IsProgramInput(startCell))
  {
    current = new RbtToken(contextp->GetVble(startCell));
  }
  else
  {
    RbtCommands comm = chrom->GetCommand(startCell);
    if (comm.IsRCte())
    {
      RbtDouble val = chrom->GetRCte(startCell);
      contextp->Assign(startCell, val);
      current = new RbtToken(contextp->GetVble(startCell));
    }
    else
    {
      current = new RbtToken(comm);
      RbtInt nargs = comm.GetNArgs();
      for (RbtInt i = nargs - 1 ; i >= 0 ; i--)
        cells.push(chrom->GetArgument(startCell, i));
    }
  }
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtCellTokenIter::RbtCellTokenIter(const RbtCellTokenIter& ti) 
  : chrom(ti.chrom), cells(ti.cells), current(ti.current)
{
  _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtCellTokenIter::~RbtCellTokenIter()
{
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}
  
void RbtCellTokenIter::Next(RbtContextPtr contextp)
{
  if (cells.size() > 0)
  {
    RbtInt cell = cells.top();
    cells.pop();
    if (chrom->IsProgramInput(cell))
    {
      current = new RbtToken(contextp->GetVble(cell));
    }
    else
    {
      RbtCommands comm = chrom->GetCommand(cell);
      if (comm.IsRCte())
      {
        RbtDouble val = chrom->GetRCte(cell);
        contextp->Assign(cell, val);
        current = new RbtToken(contextp->GetVble(cell));
      }
      else
      {
        current = new RbtToken(comm);
        RbtInt nargs = comm.GetNArgs();
        for (RbtInt i = nargs-1 ; i >= 0 ; i--)
          cells.push(chrom->GetArgument(cell, i));
      }
    }
  }
}

RbtTokenPtr RbtCellTokenIter::Current() { return current;}

void RbtCellTokenIter::copy(const RbtCellTokenIter & ti) 
{
  chrom = ti.chrom;
  current = ti.current;
  cells = ti.cells;
}

