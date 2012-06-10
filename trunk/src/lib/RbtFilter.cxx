/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtFilter.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include <sstream>
using std::ostringstream;
using std::ends;

#include "RbtFilter.h"
#include "RbtParser.h"
#include "RbtStringTokenIter.h"
#include "RbtWorkSpace.h"


const int STOP = 0;
const int CONT = -1;
const int NEXT = 1;

//Static data members
RbtString RbtFilter::_CT("RbtFilter");
////////////////////////////////////////
//Constructors/destructors

// filter is used to determining if strfilter actually contains the
// filter (filter = true) or if strfilter is the name of a file
// that contains the filter (filter = false) This is the most
// common, so filter's value by default is false
RbtFilter::RbtFilter(RbtString strfilter, RbtBool filter) 
                     : RbtBaseObject(_CT, "Filter")
{
#ifdef _DEBUG
  cout << _CT << " default constructor" << endl;
#endif //_DEBUG
//  RbtString filterfilen = GetParameter("_FILTER_FILE");
  SmartPtr <istream> filterfile;
  if (filter)  // filterfilen is the filter
    filterfile = new istrstream(strfilter.c_str());
  else         // assume the strfilter is the name of the file 
               // where the filter is
    filterfile = new ifstream(strfilter.c_str(), ios_base::in);
  (*filterfile) >> nTermFilters;
  filteridx = 0;
  contextp = RbtContextPtr(new RbtStringContext(filterfile));
  RbtParser p;
  for (RbtInt i = 0 ; i < nTermFilters ; i++)
  {
    cout << "\n------------- Terminate filter " << i << "------------" << endl;
    RbtString s;
    getline(*filterfile, s, ',');
    SmartPtr <istream> istrp (new istrstream(s.c_str()));
    RbtTokenIterPtr ti(new RbtStringTokenIter(istrp,contextp));
    RbtFilterExpressionPtr filter = p.Parse(ti, contextp);
    PrettyPrintVisitor visitor1(contextp);
    filter->Accept(visitor1);
    terminationFilters.push_back(filter);
  }
  (*filterfile) >> nWriteFilters;
  for (RbtInt i = 0 ; i < nWriteFilters ; i++)
  {
    cout << "\n------------- Write filter -----------------" << endl;
    RbtString s;
    getline(*filterfile, s, ',');
    SmartPtr <istream> istrp (new istrstream(s.c_str()));
    RbtTokenIterPtr ti(new RbtStringTokenIter(istrp,contextp));
    RbtFilterExpressionPtr filter = p.Parse(ti, contextp);
    PrettyPrintVisitor visitor1(contextp);
    filter->Accept(visitor1);
    writtingFilter.push_back(filter);
  }
  maxnruns = 1000;
  cout << endl;
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtFilter::~RbtFilter()
{
#ifdef _DEBUG
  cout << _CT << "  destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
//Public methods
////////////////

RbtModelPtr RbtFilter::GetReceptor() const {return m_spReceptor;}
RbtModelPtr RbtFilter::GetLigand() const {return m_spLigand;}
    
//Override RbtObserver pure virtual
//Notify observer that subject has changed
void RbtFilter::Update(RbtSubject* theChangedSubject) 
{
  RbtWorkSpace* pWorkSpace = GetWorkSpace();
  if (theChangedSubject == pWorkSpace) 
  {
    //Check if receptor has been updated (model #0)
    if (pWorkSpace->GetNumModels() >= 1) 
    {
      RbtModelPtr spReceptor = GetWorkSpace()->GetModel(0);
      if (spReceptor != m_spReceptor) 
      {
#ifdef _DEBUG
        cout << "RbtBaseInterSF::Update(): Receptor has been updated" << endl;
#endif //_DEBUG
        m_spReceptor = spReceptor;
        SetupReceptor();
      }
    }
    //Check if ligand has been updated (model #1)
    if (pWorkSpace->GetNumModels() >= 2) 
    {        
      RbtModelPtr spLigand = GetWorkSpace()->GetModel(1);
      if (spLigand != m_spLigand) 
      {
#ifdef _DEBUG
        cout << "RbtBaseInterSF::Update(): Ligand has been updated" << endl;
#endif //_DEBUG
        m_spLigand = spLigand;
        SetupLigand();
      }
    }
  }
}

//Called by Update when receptor is changed
void RbtFilter::SetupReceptor()
{
  // read filter file
  ((RbtStringContextPtr)contextp)->UpdateSite(m_spReceptor, 
                                              GetWorkSpace()->GetDockingSite());
}

//Called by Update when ligand is changed
void RbtFilter::SetupLigand()
{
  ((RbtStringContextPtr)contextp)->UpdateLigs(m_spLigand);
  filteridx = 0;
  nruns = 1;
  contextp->Assign("SCORE.NRUNS", nruns);
}

//Called by Update when either model has changed
void RbtFilter::SetupScore()
{
  ((RbtStringContextPtr)contextp)->UpdateScores(GetWorkSpace()->GetSF(), GetWorkSpace()->GetModel(1));
    // write down SCORE.NRUNS to make sure is getting the
    // right value
//  exit(1);
}

//Finished with ligand?
RbtBool RbtFilter::Terminate()
{
  SetupScore();
  RbtBool bTerm;
  if (nTermFilters > 0) {
    EvaluateVisitor visitor2(contextp);
    terminationFilters[filteridx]->Accept(visitor2);
    //    cout << filteridx << "\t" 
    //         << terminationFilters[filteridx]->GetValue() 
    //         << "\tnruns: " << nruns << endl;
    RbtDouble val = terminationFilters[filteridx]->GetValue();
    if (val == STOP) {
      if (GetTrace() > 1) {
	cout << "Terminate with this ligand\n";
      }
      bTerm = true;
    }
    else if (val == NEXT) {
      // keep checking filters until no more or one of them is not passed
      // by this conformation. I do not check if it return stop
      // or cont, just next filter. 
      while ((val == NEXT) && (filteridx < nTermFilters)) {
        filteridx++;
	if (filteridx < nTermFilters) {
          nruns = 0; // it should not stop because of NRUNS
          contextp->Assign("SCORE.NRUNS", nruns);
          terminationFilters[filteridx]->Accept(visitor2);
          val = terminationFilters[filteridx]->GetValue();
	  if (GetTrace() > 1) {
	    cout << "Go to next phase\n";
	  }
	}
      }
      if (filteridx == nTermFilters) {
	if (GetTrace() > 1) {
	  cout << "Terminate with this ligand\n";
	}
        bTerm = true; // output of filter is next phase, but there are no
                      // more phases
      }
      else {
        nruns = 1;
        contextp->Assign("SCORE.NRUNS", nruns);
        bTerm = false;
      }
    }
    else if (val == CONT) {
      // finish because it's reached limit runs 
      if (nruns >= maxnruns) {
	if (GetTrace() > 1) {
	  cout << "Terminate with this ligand\n";
	}
        bTerm = true;
      }
      else
      {
        nruns++;
        contextp->Assign("SCORE.NRUNS", nruns);
	if (GetTrace() > 1) {
	  cout << "Continue in this phase\n";
	}
        bTerm = false;
      }
    }
    else
    {
      ostringstream error;
      error << "Wrong output: " << val;
      error << "; Termination filters should return 0,1 or -1" << ends;
      throw RbtError(_WHERE_, error.str());
    }
  }
  else
    bTerm = true; // if there are no Term filters, no need to write anything
                  // Simply terminate.
  return bTerm;
}

//Output conformation?
RbtBool RbtFilter::Write()
{

  RbtBool bWrite = true;
  for (RbtInt i = 0 ; i < nWriteFilters ; i++)
  {
    EvaluateVisitor visitor2(contextp);
    writtingFilter[i]->Accept(visitor2);
//    cout << writtingFilter[i]->GetValue() << endl;
    RbtDouble val = writtingFilter[i]->GetValue();
    if (val >= 0.0)
    {
      bWrite = false;
      break;
    }
    // else it evaluates next filter. If there are no more, all
    // have returned a negative value and it returns true
  }
  //if (bWrite)
  //  cout << "write this ligand\n";
  //else
  //  cout << "Do not write this ligand\n";
  return bWrite;
}


