//Ancillary Software from Vernalis
//$Id: //depot/dev/client3/rdock/2006.1/src/exe/rbmoegrid.cxx#4 $

#include <fstream>
#include <iomanip>
#include <unistd.h>

#include "RbtBiMolWorkSpace.h"
#include "RbtParameterFileSource.h"
#include "RbtPRMFactory.h"
#include "RbtSFFactory.h"
#include "RbtRealGrid.h"
#include "RbtTriposAtomType.h"
#include "RbtMOEGrid.h"

const RbtString EXEVERSION = " ($Id: //depot/dev/client3/rdock/2006.1/src/exe/rbmoegrid.cxx#4 $)";
const RbtString _ROOT_SF = "SCORE";

//Creates list of probe models
//NOTE: MUST BE IN ORDER OF ASCENDING RADII
RbtModelList CreateProbes(RbtString anAtomTypeStr) {
	RbtModelList			probes;
	RbtTriposAtomTypeList	atomTypes;
	RbtTriposAtomType		tAtomType;
	atomTypes.push_back(RbtTriposAtomType::UNDEFINED);
	//	keep is as a vector though there is only one element in there
	for (RbtTriposAtomTypeListConstIter iter = atomTypes.begin(); iter != atomTypes.end(); iter++) {
		RbtAtomList atomList;
		RbtBondList bondList;
		RbtAtomPtr spAtom(new RbtAtom(1));
		spAtom->SetTriposType(tAtomType.Str2Type(anAtomTypeStr) );
		if(spAtom->GetTriposType() == RbtTriposAtomType::UNDEFINED)
			throw RbtBadArgument(_WHERE_,"Undefined atom type (try standard Tripos types)");
		spAtom->SetAtomicMass(12.0);//Mass is irrelevant as we are not rotating models around COM
		atomList.push_back(spAtom);
		probes.push_back(new RbtModel(atomList,bondList));
	}
	return probes;  	
}

void PrintUsage(void)
{
	cout << endl << "rbmoegrid - calculates grids for a given atom type" << endl;
	cout << endl << "Usage:\trbmoegrid -o <OutputRoot> -r <ReceptorPrmFile> -p <SFPrmFile> [-g <GridStep> -b <border> -t <tripos_type>]" << endl;
	cout << endl << "Options:\t-o <OutFileName> (.grd is suffixed)" << endl;
	cout << "\t\t-r <ReceptorPrmFile> - receptor param file (contains active site params)" << endl;
	cout << "\t\t-p <SFPrmFile> - scoring function param file (default calcgrid_vdw.prm)" << endl;
	cout << "\t\t-g <GridStep> - grid step (default=0.5A)" << endl;
	cout << "\t\t-b <Border> - grid border around docking site (default=1.0A)" << endl;
	cout << "\t\t-t <AtomType> - Tripos atom type (default is C.3)"<< endl;
}

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{
	cout.setf(ios_base::left,ios_base::adjustfield);

	//Strip off the path to the executable, leaving just the file name
	RbtString strExeName(argv[0]);
	RbtString::size_type i = strExeName.rfind("/");
	if (i != RbtString::npos)
		strExeName.erase(0,i+1);

	//Print a standard header
	Rbt::PrintStdHeader(cout,strExeName+EXEVERSION);

	//Command line arguments and default values
	RbtString strSuffix(".grd");					// file out suffix
	RbtString strOutfName("moegrid");				// default output filename
	vector<RbtString> strReceptorPrmFiles;			// Receptor param files (you can have more!)
	RbtString strSFFile("calcgrid_vdw.prm");		// Scoring function file
	RbtDouble gs(0.5);								// grid step
	RbtDouble border(1.0);							// grid border around docking site

	//Brief help message
	if (argc == 1) {
		PrintUsage();
		exit(1);
	}
	// parsing command-line options
	opterr 	= 0;
	int			c;
	string		optStr("");
	RbtString	strGridStep;	// due to switch scope better to declare here
	RbtString	strBorder;
	RbtString	strAtomType("C.3");
	while( (c=getopt(argc,argv,"o:r:p:g:b:t:")) != -1) {
		switch (c) {
			case 'o':
				cout << "\t -o "<<optarg<<endl;
				strOutfName = optarg;
				break;
			case 'r':
				cout << "\t -r "<<optarg<<endl;	
				strReceptorPrmFiles.push_back(optarg);
				break;
			case 'p':
				cout << "\t -p "<<optarg<<endl;
				strSFFile = optarg;
				break;
			case 'g':
				cout << "\t -g "<<optarg<<endl;
				strGridStep	= optarg;
				gs 			= atof(strGridStep.c_str());
				break;
			case 'b':
				cout << "\t -b "<<optarg<<endl;	
				strBorder	= optarg;
				border		= atof(strBorder.c_str());
				break;
			case 't':
				cout << "\t -t "<<optarg<<endl;
				strAtomType	= optarg;
				break;
			case '?':
				if (isprint(optopt)) {
					optStr+=optopt;
					std::cerr << "Unknown / missing option `-"<< optStr <<"'"<<endl;
				}
				else
					std::cerr << "Garbage option character. "<<endl;
				return 1;
			default:
				PrintUsage();
				exit(1);
		}
	}
	cout << endl;

	// start processing receptors
	try {
		// generate output file name
		RbtString	strOutputFile(strOutfName+strSuffix);
		RbtMOEGrid	theMOEGrid;
		theMOEGrid.SetOutputFileName(strOutputFile);
		// read each receptors one-by-one to get a common set of grid extents and origin
		theMOEGrid.CalculateCommonExtents(strReceptorPrmFiles);
		RbtCoord	minCoord	= theMOEGrid.GetMinExtents()-border;
		RbtCoord	maxCoord	= theMOEGrid.GetMaxExtents()+border;
		RbtVector	gridStep(gs,gs,gs);
		RbtVector	recepExtent = maxCoord-minCoord;
		RbtUInt		nX			= int(recepExtent.x/gridStep.x)+1;
		RbtUInt		nY			= int(recepExtent.y/gridStep.y)+1;
		RbtUInt		nZ			= int(recepExtent.z/gridStep.z)+1;
		cout << "Constructing grid of size " << nX << " x " << nY << " x " << nZ << endl;
		vector<RbtString>::iterator	strReceptorPrmFilesIter;
		for(	strReceptorPrmFilesIter = strReceptorPrmFiles.begin();
				strReceptorPrmFilesIter != strReceptorPrmFiles.end();
				++strReceptorPrmFilesIter) {
			//Create a bimolecular workspace
			RbtBiMolWorkSpacePtr spWS(new RbtBiMolWorkSpace());
			//Set the workspace name to the root of the receptor .prm filename
			RbtStringList componentList = Rbt::ConvertDelimitedStringToList((*strReceptorPrmFilesIter),".");
			RbtString wsName = componentList.front();
			spWS->SetName(wsName);

			//Read the receptor parameter file
			RbtParameterFileSourcePtr spRecepPrmSource(new RbtParameterFileSource(Rbt::GetRbtFileName("data/receptors",(*strReceptorPrmFilesIter))));
			cout << endl << "RECEPTOR:" << endl << spRecepPrmSource->GetFileName() << endl << spRecepPrmSource->GetTitle() << endl;

			//Read the scoring function file
			RbtParameterFileSourcePtr spSFSource(new RbtParameterFileSource(Rbt::GetRbtFileName("data/sf",strSFFile)));
			RbtSFFactoryPtr spSFFactory(new RbtSFFactory());//Factory class for scoring functions
			RbtSFAggPtr spSF(spSFFactory->CreateAggFromFile(spSFSource,_ROOT_SF));//Root SF aggregate

			//Register the scoring function with the workspace
			//Dump details to cout
			spWS->SetSF(spSF);
			cout << endl << "SCORING FUNCTION DETAILS:" << endl << *spSF << endl;

			//Create the receptor model from the file names in the receptor parameter file
			spRecepPrmSource->SetSection();
			RbtPRMFactory prmFactory(spRecepPrmSource);
			RbtModelPtr spReceptor = prmFactory.CreateReceptor();

			//Read docking site from file and register with workspace
			RbtString strASFile = spWS->GetName()+".as";
			RbtString strInputFile = Rbt::GetRbtFileName("data/grids",strASFile);
			//DM 26 Sep 2000 - ios_base::binary is invalid with IRIX CC
#if defined(__sgi) && !defined(__GNUC__)
			ifstream istr(strInputFile.c_str(),ios_base::in);
#else
			ifstream istr(strInputFile.c_str(),ios_base::in|ios_base::binary);
#endif
			RbtDockingSitePtr spDS(new RbtDockingSite(istr));
			istr.close();
			spWS->SetDockingSite(spDS);

			//Register receptor with workspace
			spWS->SetReceptor(spReceptor);
			cout << endl << "DOCKING SITE" << endl << (*spDS) << endl;

			//Create a grid covering the docking site, plus user-defined border
			//RbtCoord minCoord = spDS->GetMinCoord()-border;
			//RbtCoord maxCoord = spDS->GetMaxCoord()+border;
			//RbtVector recepExtent = maxCoord-minCoord;
			//RbtVector gridStep(gs,gs,gs);
			//RbtUInt nX = int(recepExtent.x/gridStep.x)+1;
			//RbtUInt nY = int(recepExtent.y/gridStep.y)+1;
			//RbtUInt nZ = int(recepExtent.z/gridStep.z)+1;
			//cout << "Constructing grid of size " << nX << " x " << nY << " x " << nZ << endl;
			RbtRealGridPtr spGrid(new RbtRealGrid(minCoord,gridStep,nX,nY,nZ));
			float* gridData = spGrid->GetGridData();

			//Create probes
			RbtModelList probes = CreateProbes(strAtomType);

			// write the grid into a MOE grid file
			cout << "=================================="<<endl;
			cout << "Generating MOE grid "<<endl;
			int dim		=	3;	// its a 3D grid
			// grid extents are in minCoord and maxCoord
			vector<double>	grid_origin;
			grid_origin.push_back(minCoord.x);
			grid_origin.push_back(minCoord.y);
			grid_origin.push_back(minCoord.z);
			vector<double>	grid_extents;
			grid_extents.push_back(maxCoord.x);
			grid_extents.push_back(maxCoord.y);
			grid_extents.push_back(maxCoord.z);
			// create MOE grid shape vector
			RbtMOEGridShape	theShape(grid_origin, grid_extents, gs);
			// MOE grid data vector
			RbtMOEGridData	theData;
			// reserve memory for data
			theData.reserve(theShape.GetDataSize());
			//Store regular pointers to avoid smart pointer dereferencing overheads
			RbtRealGrid* pGrid(spGrid);
			RbtSFAgg* pSF(spSF);
			RbtTriposAtomType triposType;
			//Main loop over each probe model
			for (RbtModelListConstIter mIter = probes.begin(); mIter != probes.end(); mIter++) {
				RbtModelPtr spLigand(*mIter);
				RbtAtom* pAtom = spLigand->GetAtomList().front();
				RbtTriposAtomType::eType atomType = pAtom->GetTriposType();
				RbtString strType = triposType.Type2Str(atomType);
				cout << "Atom type=" << strType << endl;
				//Register ligand with workspace
				spWS->SetLigand(spLigand);
				//Loop over all grid coords and calculate the score at each position		
				for (RbtUInt i = 0; i < spGrid->GetN(); i++ ) {
					pAtom->SetCoords(pGrid->GetCoord(i));
					RbtMOEGridPoint p;
					p.SetValue(pSF->Score());
					theData.push_back(p);
				}
			}
			// create grid object from shape and data 
			theMOEGrid.SetShape(theShape);
			theMOEGrid.SetData(theData);
			cout << "writing grid... ";
			// if it is the very first, then override
			if(strReceptorPrmFilesIter == strReceptorPrmFiles.begin()) {
				theMOEGrid.WriteGrid(ios_base::trunc);
			// otherwise append
			} else {
				theMOEGrid.WriteGrid(ios_base::app);
			}
			cout << "done."<<endl;
		}
	}
	catch (RbtError& e) {
		cout << e << endl;
	}
	catch (...) {
		cout << "Unknown exception" << endl;
	}

	_RBTOBJECTCOUNTER_DUMP_(cout)

		return 0;
}
