/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtCellTokenIter.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

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

