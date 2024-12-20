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

#include "RbtAlignTransform.h"

#include "RbtDockingSite.h"
#include "RbtWorkSpace.h"

// Static data member for class type
RbtString RbtAlignTransform::_CT("RbtAlignTransform");
// Parameter names
RbtString RbtAlignTransform::_COM("COM");
RbtString RbtAlignTransform::_AXES("AXES");

const RbtAlignTransform::Config
    RbtAlignTransform::DEFAULT_CONFIG{};  // Empty initializer to fall back to default values

RbtAlignTransform::RbtAlignTransform(const RbtString& strName, const Config& config):
    RbtBaseBiMolTransform(_CT, strName),
    m_rand(Rbt::GetRbtRand()),
    m_totalSize(0),
    config{config} {
    DEBUG_ERR(_CT << " parameterised constructor" << endl);
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtAlignTransform::~RbtAlignTransform() {
    DEBUG_ERR(_CT << " destructor" << endl);
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Protected methods
///////////////////
void RbtAlignTransform::SetupReceptor() {
    m_cavities.clear();
    m_cumulSize.clear();
    m_totalSize = 0;
    // Return now if docking site not registered with workspace
    RbtDockingSitePtr spDS = GetWorkSpace()->GetDockingSite();
    if (spDS.Null()) return;

    m_cavities = spDS->GetCavityList();
    if (m_cavities.empty()) return;

    // Determine cumulative sizes, for selecting cavities with a probability proportional to their size
    for (RbtCavityListConstIter cIter = m_cavities.begin(); cIter != m_cavities.end(); cIter++) {
        m_totalSize += (*cIter)->GetNumCoords();
        m_cumulSize.push_back(m_totalSize);
    }
}

void RbtAlignTransform::SetupLigand() {}     // No ligand setup required
void RbtAlignTransform::SetupTransform() {}  // No further setup required

////////////////////////////////////////
// Private methods
///////////////////
// Pure virtual in RbtBaseTransform
// Actually apply the transform
void RbtAlignTransform::Execute() {
    RbtModelPtr spLigand(GetLigand());
    if (spLigand.Null() || m_cavities.empty()) return;

    RbtInt iTrace = GetTrace();

    // Select a cavity at random, weighted by each cavity size
    RbtInt iRnd = m_rand.GetRandomInt(m_totalSize);
    RbtInt iCavity(0);
    for (iCavity = 0; iRnd >= m_cumulSize[iCavity]; iCavity++)
        ;
    RbtCavityPtr spCavity = m_cavities[iCavity];

    // Get the coord list and principal axes for the cavity
    const RbtCoordList& coordList = spCavity->GetCoordList();
    const RbtPrincipalAxes& prAxes = spCavity->GetPrincipalAxes();

    switch (config.center_of_mass_placement_strategy) {
        case COM_RANDOM:
            if (!coordList.empty()) {
                RbtInt iRand = m_rand.GetRandomInt(coordList.size());
                RbtCoord asCavityCoord = coordList[iRand];
                if (iTrace > 1) {
                    cout << "Translating ligand COM to active site coord #" << iRand << ": " << asCavityCoord << endl;
                }
                // Translate the ligand center of mass to the selected coord
                spLigand->SetCenterOfMass(asCavityCoord);
            }
            break;
        case COM_ALIGN:
            if (iTrace > 1) cout << "Translating ligand COM to active site COM: " << prAxes.com << endl;
            spLigand->SetCenterOfMass(prAxes.com);
            break;
        default:
            throw RbtBadArgument(_WHERE_, "Bad ligand center of mass placement strategy");
            break;
    }

    switch (config.axes_alignment_strategy) {
        case AXES_RANDOM:
            {  // Braces to prevent compiler complains due to variables initialized inside this case

                RbtDouble thetaDeg = 180.0 * m_rand.GetRandom01();
                RbtCoord axis = m_rand.GetRandomUnitVector();
                if (iTrace > 1)
                    cout << "Rotating ligand by " << thetaDeg << " deg around axis=" << axis << " through COM" << endl;
                spLigand->Rotate(axis, thetaDeg);
                break;
            }
        case AXES_ALIGN:
            spLigand->AlignPrincipalAxes(prAxes, false);  // false = don't translate COM as we've already done it above
            if (iTrace > 1) cout << "Aligning ligand principal axes with active site principal axes" << endl;
            // Make random 180 deg rotations around each of the principal axes
            if (m_rand.GetRandom01() < 0.5) {
                spLigand->Rotate(prAxes.axis1, 180.0, prAxes.com);
                if (iTrace > 1) cout << "180 deg rotation around PA#1" << endl;
            }
            if (m_rand.GetRandom01() < 0.5) {
                spLigand->Rotate(prAxes.axis2, 180.0, prAxes.com);
                if (iTrace > 1) cout << "180 deg rotation around PA#2" << endl;
            }
            if (m_rand.GetRandom01() < 0.5) {
                spLigand->Rotate(prAxes.axis3, 180.0, prAxes.com);
                if (iTrace > 1) cout << "180 deg rotation around PA#3" << endl;
            }
            break;
        default:
            throw RbtBadArgument(_WHERE_, "Bad ligand axes alignment strategy");
            break;
    }
}
