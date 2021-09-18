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

#include <fstream>
#include <cstdlib>

#include "RbtMOEGrid.h"

#include "RbtParameterFileSource.h"
#include "RbtRealGrid.h"

const	string	GRID_SIGNATURE	= "SVLgrid\n";
const	int		VERSION			= 100;
const	int		SKIP			= 0;

// utility SVL-related functions that are not really OO and are not strictly 
// grid-related only. still the main purpose of these to export numbers into 
// format that is readable by MOE. If doubt, use 
// $ hexdump -dc test.grd |head -30 
// to make sure the generated grid is the same as the "official" MOE version

namespace svl {
	void write_int4m(std::ofstream& ofstr, int an_int){
		char *i_string = (char*)&an_int;
		ofstr.put(i_string[3]);
		ofstr.put(i_string[2]);
		ofstr.put(i_string[1]);
		ofstr.put(i_string[0]);
	}

	void write_IEEE8m(std::ofstream& ofstr, double a_double){
		char *d_string = (char*)&a_double;
		ofstr.put(d_string[7]);
		ofstr.put(d_string[6]);
		ofstr.put(d_string[5]);
		ofstr.put(d_string[4]);
		ofstr.put(d_string[3]);
		ofstr.put(d_string[2]);
		ofstr.put(d_string[1]);
		ofstr.put(d_string[0]);
	}
}


//////////////////////////////////////////////////////////////////////////////////////
// RbtMOEGridPoint 
// singe grid point related methods
//////////////////////////////////////////////////////////////////////////////////////
//
// default ctor
RbtMOEGridPoint::RbtMOEGridPoint()
{
	value	=	1.0;
	coords.push_back(0.0);	// it is a 1D single point
}

// parameterised ctor receives an n-dimensional vector and a value as
// input
RbtMOEGridPoint::RbtMOEGridPoint(vector<double> a_coord, double a_value)
{
	// copy everything
	vector<double>::iterator crditer;
	for(crditer=a_coord.begin(); crditer != a_coord.end(); ++crditer) 
		coords.push_back(*crditer);
	value	= a_value;
}

// convenience method for RbtCoord
RbtMOEGridPoint::RbtMOEGridPoint(RbtCoord& a_coord, double a_value)
{
	coords.push_back(Rbt::ExtractXCoord(a_coord));
	coords.push_back(Rbt::ExtractYCoord(a_coord));
	coords.push_back(Rbt::ExtractZCoord(a_coord));
	value	= a_value;
}
//////////////////////////////////////////////////////////////////////////////////////
// RbtMOEGridShape
// SVL grid shape vector related methods
//////////////////////////////////////////////////////////////////////////////////////
// default ctor allocates a 1D grid starting at 0.0 to 1.0 
// with 1.0 stepsize (2 elements at least needed for MOE)
RbtMOEGridShape::RbtMOEGridShape() 
{
	vector<double>	grid_origin(1,0.0);
	vector<double>	grid_extents(1,1.0);
	double			grid_step	= 1.0;
	RbtMOEGridShape(grid_origin,grid_extents,grid_step);
}
	
// parametrized ctor receives grid as the origin vector, the extent
// vector and stepsize
RbtMOEGridShape::RbtMOEGridShape(vector<double> grid_origin, vector<double> grid_extents, double grid_steps) 
{
	// grid origin should to be a N-dimensional point with the 
	// same dimensions as the extends
	if(grid_origin.size() != grid_extents.size())
		cout <<"Warning: grid origin and extents dimension mismatch." <<endl;
	// copying ctor parameter vectors (it is not guaranteed they will
	// live forever during the code)
	vector<double>::iterator o_iter,e_iter;	// for the O_rigin and the E_xtents
	for(o_iter=grid_origin.begin(); o_iter!=grid_origin.end(); ++o_iter) 
		origin.push_back(*o_iter);
	for(e_iter=grid_extents.begin(); e_iter!=grid_extents.end(); ++e_iter) 
		extents.push_back(*e_iter);
	step_size	= grid_steps;

	// calculate sizes for all vertices of the grid cube
	for(	o_iter=origin.begin(), e_iter=extents.begin(); 
			o_iter != origin.end() && e_iter != extents.end(); 
			++o_iter, ++e_iter) {
		i_extents.push_back( (int) (((*e_iter)-(*o_iter) )/step_size)+1 );
	}
	// calculate the total number of points (initialize multiplication by 1)
	data_size	=	1;
	for( 	vector<int>::const_iterator ext_iter	=  i_extents.begin();
			ext_iter								!= i_extents.end();
			++ext_iter)
		data_size	*= (*ext_iter);
}

//////////////////////////////////////////////////////////////////////////////////////
// RbtMOEGrid
// grid IO methods involving both the data and the shape vectors
//////////////////////////////////////////////////////////////////////////////////////

// default ctor
RbtMOEGrid::RbtMOEGrid()
{
	// default stream name (should be overridden by user in command-line)
	stream_name = "moegrid.grd";
}

// parameterised ctor needs the shape and the data vector
RbtMOEGrid::RbtMOEGrid(RbtMOEGridShape& a_shape, RbtMOEGridData& a_data, RbtString a_strfName="moegrid.grd")
{
	stream_name = a_strfName;
	myShape		= a_shape;
	myData		= a_data;
}

// read a single rDock .as /.prm file to get the docking site extents
void RbtMOEGrid::GetDockingSiteExtents(RbtString& a_strPrmFile)
{
	//Create a bimolecular workspace
	RbtBiMolWorkSpacePtr spWS(new RbtBiMolWorkSpace());
	//Set the workspace name to the root of the receptor .prm filename
	RbtStringList componentList = Rbt::ConvertDelimitedStringToList(a_strPrmFile,".");
	RbtString wsName = componentList.front();
	spWS->SetName(wsName);
	//Read docking site from file
	RbtString strASFile = spWS->GetName()+".as";
	RbtString strInputFile = Rbt::GetRbtFileName("data/grids",strASFile);
#if defined(__sgi) && !defined(__GNUC__)
	ifstream istr(strInputFile.c_str(),ios_base::in);
#else
	ifstream istr(strInputFile.c_str(),ios_base::in|ios_base::binary);
#endif
	RbtDockingSitePtr dSite(new RbtDockingSite(istr));
	istr.close();
	// we need only these two values
	c_min	= dSite->GetMinCoord();
	c_max	= dSite->GetMaxCoord();
}
// reads all the .prm files and calculates the common extents for all 
// cavities to get a grid that includes them all
void RbtMOEGrid::CalculateCommonExtents(vector<RbtString> strPrmFiles)
{
	
	cout << "Receptors: "<<endl;
	vector<RbtString>::iterator	strPrmFilesIter;
	for(	strPrmFilesIter = strPrmFiles.begin();
			strPrmFilesIter != strPrmFiles.end();
			++strPrmFilesIter) {
		cout <<"\t"<<(*strPrmFilesIter)<<endl;
		// initialize values with the very first
		if(strPrmFilesIter == strPrmFiles.begin()) {
			GetDockingSiteExtents( (*strPrmFilesIter) );
			min	= c_min;
			max	= c_max;
		} // otherwise go through on the rest of the prms and calculate extremes
		else {
			GetDockingSiteExtents( (*strPrmFilesIter) );
			// change min extents to current values if required
			if(c_min.x < min.x)
				min.x	= c_min.x;
			if(c_min.y < min.y)
				min.y	= c_min.y;
			if(c_min.z < min.z)
				min.z	= c_min.z;
			// change max extents
			if(c_max.x > max.x)
				max.x	= c_max.x;
			if(c_max.y > max.y)
				max.y	= c_max.y;
			if(c_max.z > max.z)
				max.z	= c_max.z;
		}
	}
}

// writes the (binary format) grid 
// and returns the number of gridpoints written
// that should be the same as the number of datapoints (a_data.size() above)
// mode stands for ios_base::trunc or ios:app
long RbtMOEGrid::WriteGrid(std::_Ios_Openmode mode)
{
	// file to save (binary, overidden)
	ofstream	ofstr;
	ofstr.open(stream_name.c_str(),ios_base::out|mode|ios_base::binary);
	long 		written		= 0;									// data written
	long 		ndata		= (long)(myData.size());				// data vector size
	long 		dimension	= (long) myShape.GetOrigin().size();	// space dim
	// writing stuff
	// if we are starting a new grid, write signature
	if(ios_base::trunc == mode)
		ofstr.write(GRID_SIGNATURE.c_str(), GRID_SIGNATURE.size());	// signature
	svl::write_int4m(ofstr,VERSION);								// VERSION
	svl::write_int4m(ofstr,ndata);									// size of data vector
	svl::write_int4m(ofstr,dimension);								// size of shape vector/origin
	svl::write_int4m(ofstr,0);										// the number of skipped bytes (is zero)
	vector<int>	i_ext	= myShape.GetIExtents();
	for(vector<int>::iterator i_ext_iter=i_ext.begin();
		i_ext_iter != i_ext.end();
		++i_ext_iter) {
		svl::write_int4m(ofstr,(*i_ext_iter));						// write shape vector extents
	}
	vector<double>	grid_origin		= myShape.GetOrigin();			// cat shape :
	vector<double>	grid_extents	= myShape.GetExtents();			// list all the values in the shape vector
	vector<double>::iterator	o_iter,e_iter;						// origin and extents
	for(	o_iter=grid_origin.begin(), e_iter=grid_extents.begin();
			o_iter != grid_origin.end() && e_iter != grid_extents.end();
			++o_iter, ++e_iter) {
		for(	double grid_coord= (*o_iter); 
				grid_coord <= (*e_iter);
				grid_coord += myShape.GetStepSize()) {
			svl::write_IEEE8m(ofstr,(grid_coord));					// write shape vector data
		}
	}
	// write data vector
	for(	RbtMOEGridData::iterator d_iter =  myData.begin();
			d_iter != myData.end();
			++d_iter) {
			svl::write_IEEE8m(ofstr,(*d_iter).GetValue());
			++written;
	}
	ofstr.close();
	return written;
}

long RbtMOEGrid::ReadGrid()
{
	std::ifstream   ifstr(stream_name.c_str());          // file to read
	return 0;
}
