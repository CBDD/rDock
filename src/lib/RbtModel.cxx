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

#include "RbtModel.h"

#include <iomanip>

#include "RbtAtomFuncs.h"
#include "RbtBaseMolecularFileSource.h"
#include "RbtChromElement.h"
#include "RbtChromFactory.h"
#include "RbtFlexData.h"
#include "RbtModelError.h"
#include "RbtQuat.h"

RbtModel::RbtModel(RbtBaseMolecularFileSource* pMolSource): m_occupancy(1.0), m_enabled(true) {
    Create(pMolSource);
    _RBTOBJECTCOUNTER_CONSTR_("RbtModel");
}

//(Fairly) temporary constructor taking arbitrary atom and bond lists
// Use with caution
RbtModel::RbtModel(RbtAtomList& atomList, RbtBondList& bondList):
    m_pFlexData(NULL), m_pChrom(NULL), m_occupancy(1.0), m_enabled(true) {
    AddAtoms(atomList);  // Register atoms with model
    m_bondList = bondList;
    // Rbt::FindRings(m_atomList,m_bondList,m_ringList);

    // DM 8 Feb 1999 Add the unnamed coord to the coord map
    m_coordNames[""] = 0;
    // DM 15 June 2006. Define the model name as the fully qualified
    // name of the first atom (mainly used for solvent)
    if (!atomList.empty()) {
        m_strName = (atomList.front())->GetFullAtomName();
    }
    _RBTOBJECTCOUNTER_CONSTR_("RbtModel");
}

// Default destructor
RbtModel::~RbtModel() {
#ifdef _DEBUG
    cout << "~RbtModel: deleting " << m_strName << " (" << m_atomList.size() << " atoms, " << m_bondList.size()
         << " bonds)" << endl;
#endif        //_DEBUG
    Clear();  // clear the current model
    _RBTOBJECTCOUNTER_DESTR_("RbtModel");
}

//////////////////////
// Public methods
//////////////////////

// DM 12 May 1999 - support for associated model data (e.g. SD file, or generated by rbdock etc)
// Get list of field names as string list
RbtStringList RbtModel::GetDataFieldList() const {
    RbtStringList dataFieldList;
    dataFieldList.reserve(m_dataMap.size());
    for (RbtStringVariantMapConstIter iter = m_dataMap.begin(); iter != m_dataMap.end(); iter++)
        dataFieldList.push_back((*iter).first);
    return dataFieldList;
}

// Query as to whether a particular data field name is present
RbtBool RbtModel::isDataFieldPresent(const RbtString& strDataField) const {
    return m_dataMap.find(strDataField) != m_dataMap.end();
}

// Get a particular data value
// Note: unlike the MolecularFileSource method, the RbtModel version
// doesn't throw an error if the field name is not found.
// Instead, an empty variant is returned.
RbtVariant RbtModel::GetDataValue(const RbtString& strDataField) const {
    RbtStringVariantMapConstIter iter = m_dataMap.find(strDataField);
    if (iter != m_dataMap.end())
        return (*iter).second;
    else
        return RbtVariant();
}

// Set a data value (replaces existing value if field name already exists)
void RbtModel::SetDataValue(const RbtString& strDataField, const RbtVariant& dataValue) {
    m_dataMap[strDataField] = dataValue;
}

// Removes a data field completely from the data map
void RbtModel::ClearDataField(const RbtString& strDataField) { m_dataMap.erase(strDataField); }

// Removes all data fields starting with a given prefix from the data map
void RbtModel::ClearAllDataFields(const RbtString& strDataFieldPrefix) {
    RbtStringList dataFields = GetDataFieldList();
    for (RbtStringListConstIter iter = dataFields.begin(); iter != dataFields.end(); iter++) {
        if ((*iter).find(strDataFieldPrefix) == 0) ClearDataField(*iter);
    }
}

// Removes all data fields from the data map
void RbtModel::ClearAllDataFields() { m_dataMap.clear(); }

// DM 11 Jul 2000 - pseudoatom handling
// DM 17 Oct 2001 - now checks for whether pseudo atom has previously been created
// If so, returns the original pseudo atom pointer (does not create a new one)
// Reduces overhead of updating pseudo atom coordinates for flexible molecules
// Argument list has changed from RbtPseudoAtomPtr to const RbtAtomList&
// i.e. pseudo atom is created inside the method, and RbtPseudoAtomPtr returned
RbtPseudoAtomPtr RbtModel::AddPseudoAtom(const RbtAtomList& atomList) {
    // cout << "AddPseudoAtom";
    // for (RbtAtomListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    //   cout << "\t" << (*iter)->GetFullAtomName();
    // }
    // cout << endl;
    // Check if we have this pseudo atom already
    for (RbtPseudoAtomListConstIter pIter = m_pseudoAtomList.begin(); pIter != m_pseudoAtomList.end(); pIter++) {
        // cout << "Checking pseudo atom #" << (*pIter)->GetAtomId() << "\t" << (*pIter)->GetNumAtoms() << " atoms" <<
        // endl;
        RbtAtomList atomList2 = (*pIter)->GetAtomList();
        if (atomList.size() != atomList2.size()) {
            // cout << "No match: Unequal number of atoms" << endl;
            continue;
        }
        RbtBool bMatch = true;
        for (RbtAtomListConstIter aIter = atomList2.begin(); aIter != atomList2.end() && bMatch; aIter++) {
            // cout << "Checking " << (*aIter)->GetFullAtomName() << endl;
            bMatch = (Rbt::GetNumAtoms(atomList, std::bind2nd(Rbt::isAtomPtr_eq(), *aIter)) == 1);
        }
        if (bMatch) {
            // cout << "Match found" << endl;
            return (*pIter);
        }
    }
    RbtInt nPseudoAtomId = -1 - m_pseudoAtomList.size();
    RbtPseudoAtomPtr spPseudoAtom(new RbtPseudoAtom(atomList, nPseudoAtomId));
    m_pseudoAtomList.push_back(spPseudoAtom);
    // cout << "No match: creating new pseudo atom #" << nPseudoAtomId << endl;
    return spPseudoAtom;
}

void RbtModel::ClearPseudoAtoms() { m_pseudoAtomList.clear(); }

// Updates pseudoatom coords
void RbtModel::UpdatePseudoAtoms() {
    for (RbtPseudoAtomListIter iter = m_pseudoAtomList.begin(); iter != m_pseudoAtomList.end(); iter++) {
        (*iter)->UpdateCoords();
    }
}

RbtUInt RbtModel::GetNumPseudoAtoms() const { return m_pseudoAtomList.size(); }
RbtPseudoAtomList RbtModel::GetPseudoAtomList() const { return m_pseudoAtomList; }

// DM 1 Jul 2002 - tethered atom handling
RbtUInt RbtModel::GetNumTetheredAtoms() const throw(RbtError) { return GetTetheredAtomList().size(); }

// Parses the ">  <TETHERED ATOMS>" data field
//(comma-separated list of atom IDs)
// and returns an atomlist of tethered atoms in the same order
RbtAtomList RbtModel::GetTetheredAtomList() const throw(RbtError) {
    RbtAtomList tetheredAtomList;
    // Can be multi-line so process as a vector of strings (RbtStringList)
    RbtStringList dataField = GetDataValue("TETHERED ATOMS");
    for (RbtStringListConstIter iter = dataField.begin(); iter != dataField.end(); iter++) {
        // Each line is comma-separated list of atom IDs
        RbtStringList strTetheredAtoms = Rbt::ConvertDelimitedStringToList(*iter);
        for (RbtStringListConstIter iter2 = strTetheredAtoms.begin(); iter2 != strTetheredAtoms.end(); iter2++) {
            // Remember to subtract 1 from atom ID to convert to atom list index
            RbtInt i = atoi((*iter2).c_str()) - 1;
            if ((i >= 0) && (i < m_atomList.size())) {
                tetheredAtomList.push_back(m_atomList[i]);
            } else {
                throw RbtModelError(_WHERE_, "Tethered atom ID out of range (" + (*iter) + ")");
            }
        }
    }
    return tetheredAtomList;
}

void RbtModel::SetOccupancy(RbtDouble occupancy, RbtDouble threshold) {
    m_occupancy = occupancy;
    m_enabled = occupancy >= threshold;
}

// Update coords from a data source
void RbtModel::UpdateCoords(RbtBaseMolecularFileSource* pMolSource) throw(RbtError) {
    try {
        // Read coord atom list
        // DM 16 Feb 2000 - attempt to match atom names in the source with those already in the model
        // Allows updating the coords of an implicit-hydrogen model from an all-atom CRD file for e.g.
        // RbtAtomList crdAtomList(Rbt::GetMatchingAtomList(pMolSource->GetAtomList(),m_atomList));

        // DM 31 Oct 2000 - GetMatchingAtomList is too slow for large files
        // Make the assumption that atoms in the CRD file will be in the same order as the PSF,
        // but that there may be gaps
        RbtAtomList crdAtomList(pMolSource->GetAtomList());

        // Copy coords from crd file to model
        // NOTE: no checking of matching atom types at present
        RbtAtomListIter modelIter = m_atomList.begin();
        RbtAtomListIter crdIter = crdAtomList.begin();
        Rbt::isAtom_eq bIsAtomEq;

        // DM 30/11/98 Recompute the segment map in case the segment names in the CRD file are different
        m_segmentMap.clear();
        RbtInt nUpdated(0);
        while (modelIter != m_atomList.end()) {
            // DM 31 Oct 2000 - hunt for next matching atom in CRD file
            // cout << "Comparing (CRD) " << (*crdIter)->GetFullAtomName() << " with (MODEL) " <<
            // (*modelIter)->GetFullAtomName() << endl;
            while (!bIsAtomEq(*modelIter, *crdIter)) {
                //#ifdef _DEBUG
                cout << (*crdIter)->GetFullAtomName() << " does not match " << (*modelIter)->GetFullAtomName() << endl;
                //#endif //_DEBUG
                crdIter++;
            }
            (*modelIter)->SetCoords((*crdIter)->GetCoords());
            (*modelIter)->SetSegmentName((*crdIter)->GetSegmentName());
            m_segmentMap[(*modelIter)->GetSegmentName()]++;  // Increment the segment map atom counter
            modelIter++;
            crdIter++;
            nUpdated++;
        }
        UpdatePseudoAtoms();  // DM 11 Jul 2000 - need to update pseudoatom coords by hand

        // DM 1/12/98 Update the model name to reflect the coordinate file name
        m_strName = pMolSource->GetFileName();  // Filename will do as a model name for now
        // Add list of segment names in segment filter to model name
        if (pMolSource->isSegmentFilterMapDefined()) {
            RbtSegmentMap segmentFilterMap = pMolSource->GetSegmentFilterMap();
            m_strName += "::";
            m_strName += Rbt::ConvertSegmentMapToString(segmentFilterMap);
        }
        // cout << "RbtModel::UpdateCoords: " << nUpdated << " atoms in " << pMolSource->GetFileName()
        //  << " found that match atoms in model" << endl;
    }

    catch (RbtError& error) {
        // Clear();//Clear the model so we don't leave incomplete data structures hanging around
        pMolSource->Reset();  // Reset the source as we don't know what state we've left it in
        throw;                // Rethrow the RbtError
    }
}

// DM 07 Jan 1999
// Translate molecule by the given vector
void RbtModel::Translate(const RbtVector& vector) {
    Rbt::TranslateAtoms(m_atomList, vector);
    UpdatePseudoAtoms();  // DM 11 Jul 2000 - need to update pseudoatom coords by hand
}

// Rotate molecule around the given axis (through the center of mass) by theta degrees
// DM 09 Feb 1999 - Rotate now calls generic Rotate method
void RbtModel::Rotate(const RbtVector& axis, RbtDouble thetaDeg) {
    Rotate(axis, thetaDeg, GetCenterOfMass());
    UpdatePseudoAtoms();  // DM 11 Jul 2000 - need to update pseudoatom coords by hand
}

// DM 09 Feb 1999
// Rotate molecule around the given axis (through the given coordinate) by theta degrees
void RbtModel::Rotate(const RbtVector& axis, RbtDouble thetaDeg, const RbtCoord& center) {
    // Translate all atoms so that the center of rotation lies at the origin
    Rbt::TranslateAtoms(m_atomList, -center);
    // Apply a rotation through theta degrees to all atoms
    RbtQuat quat(axis, thetaDeg * M_PI / 180.0);
    Rbt::RotateAtomsUsingQuat(m_atomList, quat);
    // Translate all atoms back again so that the center of rotation is back where it started
    Rbt::TranslateAtoms(m_atomList, center);
    UpdatePseudoAtoms();  // DM 11 Jul 2000 - need to update pseudoatom coords by hand
}

// Rotate around a given bond by theta degrees, keeping the center of mass invariant
void RbtModel::RotateBond(RbtBondPtr spBond, RbtDouble thetaDeg) {
    RbtAtomPtr spAtom1 = spBond->GetAtom1Ptr();
    RbtAtomPtr spAtom2 = spBond->GetAtom2Ptr();

    // Check if both atoms are actually in the model
    if ((spAtom1->GetModelPtr() != this) || (spAtom2->GetModelPtr() != this)) return;

    // Store the center of mass
    RbtCoord com(GetCenterOfMass());

    // Coords of atom 1
    RbtCoord coord1(spAtom1->GetCoords());
    // Vector along the bond between atom 1 and atom 2 (rotation axis, doesn't have to be unit length)
    RbtVector bondVector(spAtom2->GetCoords() - coord1);

    // Translate all atoms so that atom 1 lies at the origin
    Rbt::TranslateAtoms(m_atomList, -coord1);

    // Select the atoms on one side of the bond
    // DM 30 Oct 2000 - call standalone version of ToSpin
    Rbt::ToSpin(spBond, m_atomList, m_bondList);

    // Apply a rotation through theta/2 to the selected atoms
    RbtQuat quat(bondVector, 0.5 * thetaDeg * M_PI / 180.0);
    Rbt::RotateSelectedAtomsUsingQuat(m_atomList, quat);

    // Invert the atom selection so that atoms on the other end of the bond are selected
    Rbt::InvertAtomSelectionFlags(m_atomList);
    spAtom1->SetSelectionFlag(false);
    spAtom2->SetSelectionFlag(false);
    bondVector = -bondVector;

    // Apply a rotation through -theta/2 to the selected atoms
    quat = RbtQuat(bondVector, 0.5 * thetaDeg * M_PI / 180.0);
    Rbt::RotateSelectedAtomsUsingQuat(m_atomList, quat);

    // Finally, translate the atoms back so that the center of mass is where it started
    Rbt::TranslateAtoms(m_atomList, com - GetCenterOfMass());
    UpdatePseudoAtoms();  // DM 11 Jul 2000 - need to update pseudoatom coords by hand
}

// Rotate around a given bond by theta degrees, keeping the given atom fixed (only spins the other end of the bond)
void RbtModel::RotateBond(RbtBondPtr spBond, RbtDouble thetaDeg, RbtAtomPtr spFixedAtom) {
    RbtAtomPtr spAtom1 = spBond->GetAtom1Ptr();
    RbtAtomPtr spAtom2 = spBond->GetAtom2Ptr();

    // Check if both atoms in the bond and the fixed atom are actually in the model
    if ((spAtom1->GetModelPtr() != this) || (spAtom2->GetModelPtr() != this) || (spFixedAtom->GetModelPtr() != this))
        return;

    // Coords of atom 1
    RbtCoord coord1(spAtom1->GetCoords());
    // Vector along the bond between atom 1 and atom 2 (rotation axis, doesn't have to be unit length)
    RbtVector bondVector(spAtom2->GetCoords() - coord1);

    // Translate all atoms so that atom 1 lies at the origin
    Rbt::TranslateAtoms(m_atomList, -coord1);

    // Select the atoms on one side of the bond
    // DM 30 Oct 2000 - call standalone version of ToSpin
    Rbt::ToSpin(spBond, m_atomList, m_bondList);

    // If the fixed atom is selected for rotation then we need to invert the atom selection
    // so that the other end of the bond is spun instead
    if (spFixedAtom->GetSelectionFlag()) {
        Rbt::InvertAtomSelectionFlags(m_atomList);
        spAtom1->SetSelectionFlag(false);
        spAtom2->SetSelectionFlag(false);
        bondVector = -bondVector;  // DM 25 Feb 1999 - Need to invert bond vector as well
    }

    // Apply a rotation through theta to the selected atoms
    RbtQuat quat(bondVector, thetaDeg * M_PI / 180.0);
    Rbt::RotateSelectedAtomsUsingQuat(m_atomList, quat);

    // Finally, translate the atoms back so that atom 1 is back where it started
    Rbt::TranslateAtoms(m_atomList, coord1);
    UpdatePseudoAtoms();  // DM 11 Jul 2000 - need to update pseudoatom coords by hand
}

// DM 25 Feb 1999 - Rotate around a given bond by theta degrees, only spinning one end of the bond
// If bSwap is false, spins the end bonded to atom2 in the bond
// If bSwap is true, spins the end bonded to atom1 in the bond
void RbtModel::RotateBond(RbtBondPtr spBond, RbtDouble thetaDeg, RbtBool bSwap) {
    RbtAtomPtr spAtom1 = spBond->GetAtom1Ptr();
    RbtAtomPtr spAtom2 = spBond->GetAtom2Ptr();

    // Check if both atoms in the bond and the fixed atom are actually in the model
    if ((spAtom1->GetModelPtr() != this) || (spAtom2->GetModelPtr() != this)) return;

    // Coords of atom 1
    RbtCoord coord1(spAtom1->GetCoords());
    // Vector along the bond between atom 1 and atom 2 (rotation axis, doesn't have to be unit length)
    RbtVector bondVector(spAtom2->GetCoords() - coord1);

    // Translate all atoms so that atom 1 lies at the origin
    Rbt::TranslateAtoms(m_atomList, -coord1);

    // Select the atoms on one side of the bond
    // DM 30 Oct 2000 - call standalone version of ToSpin
    Rbt::ToSpin(spBond, m_atomList, m_bondList);

    // If bSwap is true then we need to invert the atom selection
    // so that the other end of the bond is spun instead
    if (bSwap) {
        Rbt::InvertAtomSelectionFlags(m_atomList);
        spAtom1->SetSelectionFlag(false);
        spAtom2->SetSelectionFlag(false);
        bondVector = -bondVector;
    }

    // Apply a rotation through theta to the selected atoms
    RbtQuat quat(bondVector, thetaDeg * M_PI / 180.0);
    Rbt::RotateSelectedAtomsUsingQuat(m_atomList, quat);

    // Finally, translate the atoms back so that atom 1 is back where it started
    Rbt::TranslateAtoms(m_atomList, coord1);
    UpdatePseudoAtoms();  // DM 11 Jul 2000 - need to update pseudoatom coords by hand
}

void RbtModel::SaveCoords(const RbtString& coordName) {
    // Look up the coord name in the map
    RbtStringIntMapConstIter iter = m_coordNames.find(coordName);
    if (iter != m_coordNames.end()) {
        // cout << "Saving coords under name=" << iter->first << ",index=" << iter->second << endl;
        // If we find the name, reuse the existing index
        Rbt::SaveAtomCoords(m_atomList, (*iter).second);
        m_currentCoord = (*iter).second;
    } else {
        // Add a new index to the map and save the coords using this index
        RbtUInt newIdx = m_coordNames.size();
        m_coordNames[coordName] = newIdx;
        // cout << "Saving coords under name=" << coordName << ",new index=" << newIdx << endl;
        Rbt::SaveAtomCoords(m_atomList, newIdx);
        m_currentCoord = newIdx;
    }
}

void RbtModel::RevertCoords(const RbtString& coordName) throw(RbtError) {
    // Look up the coord name in the map
    RbtStringIntMapConstIter iter = m_coordNames.find(coordName);
    if (iter != m_coordNames.end()) {
        // If we find the name, revert the coords
        // cout << "Reverting coords under name=" << iter->first << ",index=" << iter->second << endl;
        Rbt::RevertAtomCoords(m_atomList, (*iter).second);
        UpdatePseudoAtoms();  // DM 11 Jul 2000 - need to update pseudoatom coords by hand
        m_currentCoord = (*iter).second;
    } else {
        // Coord name not found, don't try and revert the coords
        // cout << "Error reverting coords, name=" << coordName << " not found" << endl;
        throw RbtInvalidRequest(_WHERE_, "RevertCoords failed on model " + GetName() + " for coord name=" + coordName);
    }
}

void RbtModel::RevertCoords(RbtInt i) {
    if (i != m_currentCoord) {
        // cout << "Model: Reverting to coords #" << i << endl;
        Rbt::RevertAtomCoords(m_atomList, i);
        UpdatePseudoAtoms();
        m_currentCoord = i;
    }
}

// Returns center of mass of model
RbtCoord RbtModel::GetCenterOfMass() const { return Rbt::GetCenterOfMass(m_atomList); }

// DM 9 Nov 1999
// Returns total atomic mass (molecular weight) for the model
RbtDouble RbtModel::GetTotalAtomicMass() const { return Rbt::GetTotalAtomicMass(m_atomList); }

// DM 14 Apr 1999 - principal axes methods
// Return principal axes and center of mass for the model
RbtPrincipalAxes RbtModel::GetPrincipalAxes() const { return Rbt::GetPrincipalAxes(m_atomList); }

// Aligns the principal axes of the model to lie along alignAxes
// If required (bAlignCOM=true), also aligns the center of mass with alignAxes.com
// Default is to align with X,Y,Z Cartesian axes centered at origin
void RbtModel::AlignPrincipalAxes(const RbtPrincipalAxes& alignAxes, RbtBool bAlignCOM) {
    Rbt::AlignPrincipalAxes(m_atomList, alignAxes, bAlignCOM);
    UpdatePseudoAtoms();  // DM 11 Jul 2000 - need to update pseudoatom coords by hand
}

// DM 19 Oct 2005 - new chromosome handling
void RbtModel::SetFlexData(RbtFlexData* pFlexData) {
    m_spMutator.SetNull();
    // Manual mem management
    // DM 12 June 2006 - check that we are not passing in
    // the same RbtFlexData object before deleting the old one
    if (m_pFlexData && (m_pFlexData != pFlexData)) {
        delete m_pFlexData;
        m_pFlexData = NULL;
    }
    if (m_pChrom) {
        delete m_pChrom;
        m_pChrom = NULL;
    }
    m_pFlexData = pFlexData;
    if (m_pFlexData) {
        m_pFlexData->SetModel(this);
        RbtChromFactory chromFactory;
        m_pFlexData->Accept(chromFactory);
        m_pChrom = chromFactory.GetChrom();
        m_spMutator = chromFactory.GetModelMutator();
    }
}

RbtFlexData* RbtModel::GetFlexData() const { return m_pFlexData; }

// Returns a clone of the current chromosome for this model
// The caller has the responsibility for mem management of the clone
RbtChromElement* RbtModel::GetChrom() const { return (m_pChrom) ? m_pChrom->clone() : NULL; }

RbtBool RbtModel::isFlexible() const { return !m_spMutator.Null(); }

const RbtAtomRList& RbtModel::GetFlexIntns(RbtAtom* pAtom) const throw(RbtError) {
    if (isFlexible()) {
        // Check if atom is actually in the model
        if (pAtom->GetModelPtr() != this) {
            throw RbtBadArgument(_WHERE_,
                                 "GetFlexIntns: " + pAtom->GetFullAtomName() + " is not in model " + GetName());
        }
        const RbtAtomRListList& flexIntns(m_spMutator->GetFlexIntns());
        RbtInt id = pAtom->GetAtomId() - 1;
        // Assertion - check id is within range
        Assert<RbtAssert>(!MUT_CHECK || (id >= 0 && id < flexIntns.size()));
        RbtAtomRListListConstIter lIter = flexIntns.begin() + id;
        return *lIter;
    } else {
        throw RbtInvalidRequest(_WHERE_, "GetFlexIntns invalid for rigid models");
    }
}
RbtBondList RbtModel::GetFlexBonds() const throw(RbtError) {
    if (isFlexible()) {
        return m_spMutator->GetFlexBonds();
    } else {
        throw RbtInvalidRequest(_WHERE_, "GetFlexBonds invalid for rigid models");
    }
}

// Select all flexible interactions to the specified atom
void RbtModel::SelectFlexAtoms(RbtAtom* pAtom) {
    // First deselect all atoms in the model
    // std::for_each(m_atomList.begin(),m_atomList.end(), Rbt::SelectAtom(false));

    if (isFlexible()) {
        // Check if atom is actually in the model
        if (pAtom->GetModelPtr() != this) {
            return;
        }
        const RbtAtomRListList& flexIntns(m_spMutator->GetFlexIntns());
        RbtInt id = pAtom->GetAtomId() - 1;
        // Assertion - check id is within range
        Assert<RbtAssert>(!MUT_CHECK || (id >= 0 && id < flexIntns.size()));
        RbtAtomRListListConstIter lIter = flexIntns.begin() + id;
        std::for_each((*lIter).begin(), (*lIter).end(), Rbt::SelectAtom(true));
    }
}

// Selects all atoms that are rotated by at least one rotable bond
void RbtModel::SelectFlexAtoms() {
    // First deselect all atoms in the model
    // std::for_each(m_atomList.begin(),m_atomList.end(),Rbt::SelectAtom(false));
    if (isFlexible()) {
        const RbtAtomRListList& flexAtoms(m_spMutator->GetFlexAtoms());
        // flexAtoms is a vector of the atom lists to rotate across each rotable bond
        // If we select all of these then we create the total list of all flexible atoms
        for (RbtAtomRListListConstIter iter1 = flexAtoms.begin(); iter1 != flexAtoms.end(); iter1++) {
            std::for_each((*iter1).begin(), (*iter1).end(), Rbt::SelectAtom(true));
        }
    }
}

////////////////////////////////////////////
// Atom list functions (provided for convenience, as user could just as well
// call the Rbt:: functions with RbtModel::GetAtomList)
// e.g. RbtAtomList atomList = Rbt::GetSelectedAtomList(spModel->GetAtomList);
// is equivalent to
// RbtAtomList atomList = spModel->GetSelectedAtomList();
////////////////////////////////////////////

// Unary

// Generic template version of GetNumAtoms, passing in your own predicate
// template<class Predicate> RbtUInt RbtModel::GetNumAtoms(const Predicate& pred)
//{
//   return Rbt::GetNumAtoms(m_atomList,pred);
// }

// Generic template version of GetAtomList, passing in your own predicate
// template<class Predicate> RbtAtomList RbtModel::GetAtomList(const Predicate& pred)
//{
//   return Rbt::GetAtomList(m_atomList,pred);
// }

// Selected atoms
void RbtModel::SetAtomSelectionFlags(RbtBool bSelected) { Rbt::SetAtomSelectionFlags(m_atomList, bSelected); }

RbtUInt RbtModel::GetNumSelectedAtoms() { return Rbt::GetNumSelectedAtoms(m_atomList); }

RbtAtomList RbtModel::GetSelectedAtomList() { return Rbt::GetSelectedAtomList(m_atomList); }

// Cyclic atoms
void RbtModel::SetAtomCyclicFlags(RbtBool bCyclic) { Rbt::SetAtomCyclicFlags(m_atomList, bCyclic); }

RbtUInt RbtModel::GetNumCyclicAtoms() { return Rbt::GetNumCyclicAtoms(m_atomList); }

RbtAtomList RbtModel::GetCyclicAtomList() { return Rbt::GetCyclicAtomList(m_atomList); }

// DM 21 Jul 1999 User1 flag
void RbtModel::SetAtomUser1Flags(RbtBool bUser1) {
    for (RbtAtomListIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) (*iter)->SetUser1Flag(bUser1);
}

// DM 29 Jan 2000 User1 value
void RbtModel::SetAtomUser1Values(RbtDouble dUser1) {
    for (RbtAtomListIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) (*iter)->SetUser1Value(dUser1);
}

// DM 27 Jul 2000 User2 value
void RbtModel::SetAtomUser2Values(RbtDouble dUser2) {
    for (RbtAtomListIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) (*iter)->SetUser2Value(dUser2);
}

// Hydrogen bond acceptor atoms
RbtUInt RbtModel::GetNumHBondAcceptorAtoms() { return Rbt::GetNumHBondAcceptorAtoms(m_atomList); }

RbtAtomList RbtModel::GetHBondAcceptorAtomList() { return Rbt::GetHBondAcceptorAtomList(m_atomList); }

// Hydrogen bond donor atoms
RbtUInt RbtModel::GetNumHBondDonorAtoms() { return Rbt::GetNumHBondDonorAtoms(m_atomList); }

RbtAtomList RbtModel::GetHBondDonorAtomList() { return Rbt::GetHBondDonorAtomList(m_atomList); }

//(Formally) charged atoms
RbtUInt RbtModel::GetNumChargedAtoms() { return Rbt::GetNumChargedAtoms(m_atomList); }

RbtAtomList RbtModel::GetChargedAtomList() { return Rbt::GetChargedAtomList(m_atomList); }

// Planar atoms
RbtUInt RbtModel::GetNumPlanarAtoms() { return Rbt::GetNumPlanarAtoms(m_atomList); }

RbtAtomList RbtModel::GetPlanarAtomList() { return Rbt::GetPlanarAtomList(m_atomList); }

// Binary

// Atoms with atomic no = nAtomicNo
RbtUInt RbtModel::GetNumAtomsWithAtomicNo_eq(RbtInt nAtomicNo) {
    return Rbt::GetNumAtomsWithAtomicNo_eq(m_atomList, nAtomicNo);
}

RbtAtomList RbtModel::GetAtomListWithAtomicNo_eq(RbtInt nAtomicNo) {
    return Rbt::GetAtomListWithAtomicNo_eq(m_atomList, nAtomicNo);
}

// Atoms with FFType = strFFType
RbtUInt RbtModel::GetNumAtomsWithFFType_eq(RbtString strFFType) {
    return Rbt::GetNumAtomsWithFFType_eq(m_atomList, strFFType);
}

RbtAtomList RbtModel::GetAtomListWithFFType_eq(RbtString strFFType) {
    return Rbt::GetAtomListWithFFType_eq(m_atomList, strFFType);
}

////////////////////////////////////////////
// Bond list functions (provided for convenience, as user could just as well
// call the Rbt:: functions with RbtModel::GetBondList)
// e.g. RbtBondList bondList = Rbt::GetSelectedBondList(spModel->GetBondList);
// is equivalent to
// RbtBondList bondList = spModel->GetSelectedBondList();
////////////////////////////////////////////

// Unary

// Generic template version of GetNumBonds, passing in your own predicate
template <class Predicate>
RbtUInt RbtModel::GetNumBonds(const Predicate& pred) {
    return Rbt::GetNumBonds(m_bondList, pred);
}
// Generic template version of GetBondList, passing in your own predicate
template <class Predicate>
RbtBondList RbtModel::GetBondList(const Predicate& pred) {
    return Rbt::GetBondList(m_bondList, pred);
}

// Selected bonds
void RbtModel::SetBondSelectionFlags(RbtBool bSelected) { Rbt::SetBondSelectionFlags(m_bondList, bSelected); }

RbtUInt RbtModel::GetNumSelectedBonds() { return Rbt::GetNumSelectedBonds(m_bondList); }

RbtBondList RbtModel::GetSelectedBondList() { return Rbt::GetSelectedBondList(m_bondList); }

// Cyclic bonds
void RbtModel::SetBondCyclicFlags(RbtBool bCyclic) { Rbt::SetBondCyclicFlags(m_bondList, bCyclic); }

RbtUInt RbtModel::GetNumCyclicBonds() { return Rbt::GetNumCyclicBonds(m_bondList); }

RbtBondList RbtModel::GetCyclicBondList() { return Rbt::GetCyclicBondList(m_bondList); }

// DM 10 Dec 1998 - at some point these functions should be implemented as generic Rbt:: functions
// operating on arbitrary atom and bond lists

// D Morley, 8 Dec 1998 - modified to now call RbtAtom::isHBondDonor (rather than RbtBond::isHBondDonor
// D Morley, 2 Dec 1998 - go back to the old way, it's more convenient to get all donors in the same list
// and separate them later
// void RbtModel::GetHBondDonorLists(RbtAtomList& donorList, RbtAtomList& donorHList)
//{
// Clear the lists
//  donorList.clear();
//  donorHList.clear();

// DM 10 Dec 1998 - this returns the list of hydrogen-bonding hydrogen atoms
//  donorHList = RbtModel::GetHBondDonorAtomList();

// Now traverse the hydrogen atom list to find the matching list of heavy atoms
//  for (RbtAtomListIter donorHIter = donorHList.begin(); donorHIter != donorHList.end(); donorHIter++) {
// Get the bonded atom list for this atom (should be exactly one)
//    RbtAtomList bondedAtomList = GetBondedAtomList((*donorHIter)->GetBondMap());
//    donorList.push_back(bondedAtomList.front());//Push the heavy atom
//  }
//}

// Get min and max x,y,z coords
// DM 28 Jul 1999 - use new RbtCoordList Min,Max functions. bInit is ignored
void RbtModel::GetMinMaxCoords(RbtCoord& minCoord, RbtCoord& maxCoord, RbtBool bInit /*=true*/) {
    RbtCoordList coordList = Rbt::GetCoordList(m_atomList);
    minCoord = Rbt::Min(coordList);
    maxCoord = Rbt::Max(coordList);
}

// Get map of (key=force field atom type string, value=no. of occurrences)
RbtStringIntMap RbtModel::GetAtomTypeMap() {
    RbtStringIntMap atomTypeMap;
    for (RbtAtomListIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++)
        atomTypeMap[(**iter).GetFFType()]++;
    return atomTypeMap;
}

// Get map of (key=force field bond type (atom type pair) string, value=no. of occurrences)
RbtStringIntMap RbtModel::GetBondTypeMap() {
    RbtStringIntMap bondTypeMap;
    for (RbtBondListIter iter = m_bondList.begin(); iter != m_bondList.end(); iter++) {
        RbtString atom1Type = (*iter)->GetAtom1Ptr()->GetFFType();
        RbtString atom2Type = (*iter)->GetAtom2Ptr()->GetFFType();
        if (atom1Type > atom2Type) std::swap(atom1Type, atom2Type);
        RbtString bondPair = atom1Type + " " + atom2Type;
        bondTypeMap[bondPair]++;
    }
    return bondTypeMap;
}

//////////////////////
// Private methods
//////////////////////

// Clear the current model
void RbtModel::Clear() {
    m_strName = "";
    m_titleList.clear();

    // Set the parent model pointer to NULL for each atom, before
    // clearing the atom list, in case someone has copies of the atom list
    RbtAtomListIter iter;
    for (iter = m_atomList.begin(); iter != m_atomList.end(); iter++) (*iter)->SetModelPtr(NULL);

    m_atomList.clear();
    m_bondList.clear();
    m_segmentMap.clear();
    // Clear each ring atom list in the list of lists
    // std::for_each(m_ringList.begin(),m_ringList.end(),std::mem_fun_ref(&RbtAtomList::clear));
    for (RbtAtomListListIter liter = m_ringList.begin(); liter != m_ringList.end(); liter++) (*liter).clear();
    m_ringList.clear();    // Now clear the list of lists
    m_coordNames.clear();  // Clear map of named coords
    m_dataMap.clear();     // DM 12 May 1999 - clear associated data
    ClearPseudoAtoms();
    SetFlexData(NULL);
}

// Helper function for the constructor
// Create a new model from a data source
void RbtModel::Create(RbtBaseMolecularFileSource* pMolSource) throw(RbtError) {
    m_pFlexData = NULL;
    m_pChrom = NULL;
    Clear();  // Clear previous model, if any

    try {
        // Title list isn't exactly crucial so don't throw an error if
        // it's not supported by the source
        if (pMolSource->isTitleListSupported()) m_titleList = pMolSource->GetTitleList();

        // Atom and bond lists are mandatory for creating a model
        // so allow an error to be thrown if not supported by the source
        RbtAtomList atomList = pMolSource->GetAtomList();
        AddAtoms(atomList);  // Register atoms with model
        m_bondList = pMolSource->GetBondList();
        m_strName = pMolSource->GetFileName();  // Filename will do as a model name for now
        // Add list of segment names in segment filter to model name
        if (pMolSource->isSegmentFilterMapDefined()) {
            RbtSegmentMap segmentFilterMap = pMolSource->GetSegmentFilterMap();
            m_strName += "::";
            m_strName += Rbt::ConvertSegmentMapToString(segmentFilterMap);
        }

        // DM 12 May 1999 Read any associated data if supported by the source
        // DM 18 May 1999 Set the model name to the SD Reg Number or Name field if present
        if (pMolSource->isDataSupported()) {
            m_dataMap = pMolSource->GetDataMap();
            RbtString strRegNum = GetDataValue("REG_Number");
            RbtString strName = GetDataValue("Name");
            // In general, we want to use the Reg Number for the model name
            // except in the case of RBT compounds, where it's nice to keep the RBT prefix
            if (!strRegNum.empty()) {
                // Case 1. Reg num defined, name begins with RBT => use name
                if (strName.substr(0, 3) == "RBT") m_strName = strName;
                // Case 2. Reg num defined, name does not begin with RBT => use reg num
                else
                    m_strName = strRegNum;
            }
            // Case 3. Reg num undefined, name defined => use name
            else if (!strName.empty())
                m_strName = strName;
            // Case 4. Reg num undefined, name undefined => keep name = filename (null op)
        }

        // We've registered some of the source's atoms as our own, so reset the source
        // This will force it to recreate the atom objects next time it is used
        pMolSource->Reset();

        // 31 Oct 2000 (DM) Hack to disable ring detection
        // if $RBT_NORINGS is defined
        char* szRbtNoRings = getenv("RBT_NORINGS");
        if (szRbtNoRings == (char*)NULL) {
            Rbt::FindRings(m_atomList, m_bondList, m_ringList);
            // than set aromatic type for pi atoms. m_ringList is RbtAtomListList
            for (RbtAtomListListIter rIter = m_ringList.begin(); rIter != m_ringList.end(); rIter++) {
                int nCy = (*rIter).size();
                // DM 19 Jun 2003 - more limited definition of AROM hybrid state
                // Specifically only used for 6-membered rings where all hybrid states are initially SP2
                // DM 23 Oct 2003 - this does not work for fused rings, where the 1st ring has already
                // been assigned as AROM. Subsequent fused rings were left as SP2
                // DM 07 Nov 2003 - looser definition. Check for 6-membered rings that are all
                // pi-atoms (SP2, AROM or TRI). Change all SP2 atoms to AROM. This will change
                // C_SP2 and N_SP2 to AROM, but will not change N_TRI. Will work with thymine,
                // which under the previous definition was not changed to AROM.
                RbtInt nPi = Rbt::GetNumAtoms(*rIter, Rbt::isPiAtom());
                if ((nCy == 6) && (nCy == nPi)) {
                    for (RbtAtomListIter aIter = (*rIter).begin(); aIter != (*rIter).end(); aIter++) {
                        if ((*aIter)->GetHybridState() == RbtAtom::SP2) {
                            (*aIter)->SetHybridState(RbtAtom::AROM);
                        }
                    }
                }
            }
        }

        // DM 14 May 2002 - set the Tripos atom type property for each atom
        RbtTriposAtomType triposType;
        for (RbtAtomListIter iter = m_atomList.begin(); iter != m_atomList.end(); iter++) {
            // only if it was not set yet
            if ((*iter)->GetTriposType() == RbtTriposAtomType::UNDEFINED)
                (*iter)->SetTriposType(triposType(*iter, true));  // true = use extended types
        }

        // DM 8 Feb 1999 Add the unnamed coord to the coord map
        m_coordNames[""] = 0;
    }

    catch (RbtError& error) {
        Clear();              // Clear the model so we don't leave incomplete data structures hanging around
        pMolSource->Reset();  // Reset the source as we don't know what state we've left it in
        throw;                // Rethrow the RbtError
    }
}

void RbtModel::AddAtoms(RbtAtomList& atomList) {
    for (RbtAtomListIter iter = atomList.begin(); iter != atomList.end(); iter++) {
        // Tell the atom it belongs to this model
        (*iter)->SetModelPtr(this);
        // Add atom smart pointer to the model's atom list
        m_atomList.push_back(*iter);
        // Increment the segment map atom counter
        m_segmentMap[(*iter)->GetSegmentName()]++;
    }
}
