/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtConstraint.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include <sstream>
#include "RbtConstraint.h"
#include "RbtLigandError.h"

using std::ostringstream;
using std::ends;

// initialization of the static data of RbtConstraint
RbtString RbtConstraint::_CT("RbtConstraint");
RbtInt RbtHeavyConstraint::counter = 0;//7 Feb 2005 (DM) new constraint type
RbtInt RbtHBAConstraint::counter = 0;
RbtInt RbtHBDConstraint::counter = 0;
RbtInt RbtHydroConstraint::counter = 0;
RbtInt RbtHydroAliphaticConstraint::counter = 0;
RbtInt RbtHydroAromaticConstraint::counter = 0;
RbtInt RbtNegChargeConstraint::counter = 0;
RbtInt RbtPosChargeConstraint::counter = 0;
RbtInt RbtRingAromaticConstraint::counter = 0;

RbtConstraint::RbtConstraint(RbtCoord c, RbtDouble t)
{
    m_atomList.clear();
    coord = c;
    tolerance = t;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtConstraint::RbtConstraint(const RbtConstraint& c)
{
    coord = c.coord;
    tolerance = c.tolerance;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}
 
RbtConstraint::RbtConstraint()
{
    coord = RbtCoord(0,0,0);
    tolerance = 0.0;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}
 
RbtConstraint::~RbtConstraint()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

  
void RbtConstraint::copy(const RbtConstraint & c)
{
    coord = c.coord;
    tolerance = c.tolerance;
}

ostream& RbtConstraint::Print(ostream& o) const
{
    o << coord << "\t" << tolerance;
    return o;
}

ostream& operator<<(ostream& s, const RbtConstraint &c)
{
    return (c.Print(s));
}

//DM 17 Oct 2003 - rewrite to avoid step increase in score at the tolerance distance
//Allow for fact that atomlist may be empty
RbtDouble RbtConstraint::Score() const
{
  if (m_atomList.empty()) {
    //Return artifically high score if the ligand has none of the desired features
    //This will only apply to the optional constraints (mandatory constraints are prefiltered)
    return 99.0;
  }

  //Compile a list of all the distance**2 to the constraint center
  RbtDoubleList d2list;
  for (RbtAtomListConstIter iter = m_atomList.begin(); iter != m_atomList.end() ; iter++) {
    d2list.push_back(Rbt::Length2((*iter)->GetCoords(), coord));
  }
  //Min distance**2
  RbtDouble d2min = *(std::min_element(d2list.begin(),d2list.end()));
  //Score is zero up to the tolerance range, then quadratic beyond.
  RbtDouble dr = sqrt(d2min) - tolerance;
  return (dr > 0.0) ? dr*dr : 0.0;
}


RbtConstraintPtr Rbt::CreateConstraint(RbtCoord& c, RbtDouble& t, RbtString& n, RbtBool bCount)
{
  if (n == "Any") {
    if (bCount) RbtHeavyConstraint::counter++;
    return new RbtHeavyConstraint(c,t);
  }
  else if (n == "Don") {
    if (bCount) RbtHBDConstraint::counter++;
    return new RbtHBDConstraint(c,t);
  }
  else if (n == "Acc") {
    if (bCount) RbtHBAConstraint::counter++;
    return new RbtHBAConstraint(c,t);
  }
  else if (n == "Hyd") {
    if (bCount) RbtHydroConstraint::counter++;
    return new RbtHydroConstraint(c,t);
  }
  else if (n == "Ani") {
    if (bCount) RbtNegChargeConstraint::counter++;
    return new RbtNegChargeConstraint(c,t);
  }
  else if (n == "Cat") {
    if (bCount) RbtPosChargeConstraint::counter++;
    return new RbtPosChargeConstraint(c,t);
  }
  else if (n == "Aro") {
    if (bCount) RbtRingAromaticConstraint::counter++;
    return new RbtRingAromaticConstraint(c,t);
  }
  else if (n == "Hal") {
    if (bCount) RbtHydroAliphaticConstraint::counter++;
    return new RbtHydroAliphaticConstraint(c,t);
  }
  else if (n == "Har") {
    if (bCount) RbtHydroAromaticConstraint::counter++;
    return new RbtHydroAromaticConstraint(c,t);
  }
  else
    throw RbtError(_WHERE_,"Constraint " + n + " not recognized");
}

void Rbt::ZeroCounters()
{
  RbtHeavyConstraint::counter = 0;
  RbtHBAConstraint::counter = 0;
  RbtHBDConstraint::counter = 0;
  RbtHydroConstraint::counter = 0;
  RbtHydroAliphaticConstraint::counter = 0;
  RbtHydroAromaticConstraint::counter = 0;
  RbtNegChargeConstraint::counter = 0;
  RbtPosChargeConstraint::counter = 0;
  RbtRingAromaticConstraint::counter = 0;
}

void Rbt::ReadConstraint(istream& ifile, RbtConstraintPtr& cnt, RbtBool bCount)
{
  RbtString n;
  RbtCoord c;
  RbtDouble t;
  ifile >> c;
  ifile >> t;
  ifile >> n;
  cout << c << "\t" << t << "\t" << n << endl;
  if (ifile.fail())
    throw RbtError(_WHERE_,"Problems reading constraint " + n );
  cnt = Rbt::CreateConstraint(c,t,n,bCount);
}

void Rbt::ReadConstraintFromMoe(istream& ifile, RbtConstraintPtr& cnt, RbtBool bCount)
{
    string n, extra;
    RbtDouble x, y, z, t;
    ifile >> n;
    cout << n << endl;
    ifile >> extra;
    ifile >> x >> y >> z >> t;
    ifile >> extra >> extra;
    if (ifile.fail())
        throw RbtError(_WHERE_,"Problems reading constraint " + n );
    RbtCoord c(x,y,z);
    cout << c << "\t" << t << "\t" << n << endl;
    cnt = Rbt::CreateConstraint(c,t,n,bCount);
}

void Rbt::ReadStartMoe(istream& ifile)
{
    const int MAXLINELENGTH = 100; // I'm quite sure the max length for any 
                                   // line in moe is 72, but
                                   // prefer to be on the safe side 
    char* line = new char[MAXLINELENGTH+1];
    RbtConstraintPtr cons;
    char ch;
    ifile.getline(line, MAXLINELENGTH);  // get rid of first four lines
    ifile.getline(line, MAXLINELENGTH);
    ifile.getline(line, MAXLINELENGTH);  
        // This line has some data that could be neccesary later on, like
        // partial match information also, it contains the
        // values for the fields comment and title, that can be as
        // long as necessary. Need to make sure I read all the lines
        // before the next commentted line (starting with a #)
    while ((ch = ifile.peek()) != '#')
        ifile.getline(line, MAXLINELENGTH);
    ifile.getline(line, MAXLINELENGTH);  // next comment line
}

void Rbt::ReadConstraints(istream& ifile, RbtConstraintList& cl, RbtBool bCount)
{
    char ch;
    RbtBool moe = false;
    cl.clear();
    RbtConstraintPtr cons;
    ch = ifile.peek();
    if (ch == '#')  // moe file
    {
        moe = true;
        ReadStartMoe(ifile);
    }
    while ((!moe && !ifile.eof()) || 
           (moe  && ((ch = ifile.peek()) != '#')))
    {
        if (moe)
            ReadConstraintFromMoe(ifile, cons, bCount);
        else
            ReadConstraint(ifile, cons, bCount);
        cl.push_back(cons);
        while (ifile.get(ch))
        {
            if (!isspace(ch))
            {
                ifile.putback(ch);
                break;
            }
        }
    }
    if (cl.size() == 0)
        cout << RbtConstraint::_CT 
             << "** WARNING: number of constraints read is 0\n";
}

//07 Feb 2005 (DM) - new constraint type, any heavy atom
void RbtHeavyConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), std::not1(Rbt::isAtomicNo_eq(1)));
  if (bCheck && (m_atomList.size() < counter)) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size()
         << " heavy atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

//16 May 2003 (DM) - limit to neutral H-bond donor hydrogens
void RbtHBDConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomHBondDonor());
  m_atomList = Rbt::GetAtomList(m_atomList, std::not1(Rbt::isAtomCationic()));
  if (bCheck && (m_atomList.size() < counter)) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " neutral H-bond donor hydrogens(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

//16 May 2003 (DM) - limit to neutral H-bond acceptors
void RbtHBAConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomHBondAcceptor());
  m_atomList = Rbt::GetAtomList(m_atomList, std::not1(Rbt::isAtomAnionic()));
  if (bCheck && (m_atomList.size() < counter)) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " neutral H-bond acceptor(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

void RbtHydroConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomLipophilic());
  if (bCheck && (m_atomList.size() < counter)) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " hydrophobic atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

void RbtHydroAliphaticConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomLipophilic());
  m_atomList = Rbt::GetAtomList(m_atomList, 
				Rbt::isHybridState_eq(RbtAtom::SP3));
  if (bCheck && (m_atomList.size() < counter)) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " hydrophobic aliphatic atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

void RbtHydroAromaticConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomLipophilic());
  m_atomList = Rbt::GetAtomList(m_atomList, 
				Rbt::isHybridState_eq(RbtAtom::AROM));
  if (bCheck && (m_atomList.size() < counter)) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " hydrophobic aromatic atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

void RbtNegChargeConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomAnionic());
  if (bCheck && (m_atomList.size() < counter)) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " negatively charged atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

void RbtPosChargeConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomCationic());
  if (bCheck && (m_atomList.size() < counter)) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " positively charged atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}


void RbtRingAromaticConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
    m_atomList.clear();
    RbtAtomList at = lig->GetAtomList();
    RbtAtomListList ligRingLists = lig->GetRingAtomLists();
    for (RbtAtomListListConstIter rIter = ligRingLists.begin() ;
         rIter != ligRingLists.end(); rIter++)
    {
        if (Rbt::GetNumAtoms(*rIter,Rbt::isPiAtom()) == 
                 (*rIter).size())
        {
            RbtPseudoAtomPtr spPseudoAtom = 
                   lig->AddPseudoAtom(*rIter);
            m_atomList.push_back(spPseudoAtom);
        }
    }
    if (bCheck && (m_atomList.size() < counter))
      {
	ostringstream ostr;
	ostr << "The ligand has only " << m_atomList.size() 
	     << " aromatic ring(s) (" << counter << " required)" << ends;
        throw RbtLigandError(_WHERE_, ostr.str());
      }
}
