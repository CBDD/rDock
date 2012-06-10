/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtMOEGrid.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

// MOE grid description at ${MOE}/html/svl/manual/gridman.htm
// here you can read:
// [...]
// Creating a grid consists of creating both a shape vector and a data vector.
// The shape vector stores gridpoint positions along each dimension. The number of
// dimensions must be between 2 and 24 with each dimension containing at least 2
// points. Each element of the shape vector is a flat vector of strictly
// increasing grid positions along the corresponding dimension. The values of the
// field at the gridpoints are stored in row-major order, i.e., with the index of
// the last dimension changing most rapidly. 
// [...]
// Generating a shape vector consists of creating coordinate vectors along each of
// the n dimensions of an n-dimensional lattice. The intervals between the
// coordinates need not be regular. To generate the shape vector, you need the
// origin of the lattice and the length of each side of the lattice, along with
// the coordinate spacings in each dimension. 
// [...]
// The values of the field at the gridpoints are stored in a data vector. There
// must be one data value for every gridpoint. The length of the data vector
// must therefore equal the product of the number of gridpoints along each
// dimension. 

#include <vector>
#include <fstream>
#include <cstdlib>
#include "Rbt.h"
#include "RbtCoord.h"

#include "RbtBiMolWorkSpace.h"
// single point in a multi-dimensional grid
// the coords contains N values where N is the dimension of the grid
class RbtMOEGridPoint {
		vector<double>	coords;
		double			value;

	public:
		// default ctor with a 1D gridpoint at 0.0 value 1.0
		// see also default ctor for the shape and data vector below
		RbtMOEGridPoint();
		// parametrized ctor when supplying data as a vector + value pair
		RbtMOEGridPoint(vector<double>,double);
		// parametrized ctor getting coord values from RbtCoord (default 3D
		// grid for convenience)
		RbtMOEGridPoint(RbtCoord&, double);
		// data access methods
		double					GetValue()	{ return value; 	};
		const vector<double>&	GetCoords()	{ return coords;	};
		void					SetCoords(double a_val)	{ coords.push_back(a_val); };
		void					SetValue(double a_val)	{ value = a_val; };
};

// MOE grid shape vector
class RbtMOEGridShape
{
		// data that is provided usually in the SVL code
		// a shape vector is defined by	
		//		i)		starting coordinates
		//		ii)		X,Y,Z extension (in absolute coordinates)
		//		iii)	grid stepsize 
		// data vector defined by 
		//		i)		a 1D vector containing the same number 
		//				of elements than the total number of the
		//				shape vector

		vector<double>	origin;		// the starting point
		vector<double>	extents;	// the final points
		vector<int>		i_extents;	// extents in integers: number of 
									// points at the edges (ie 20 X 30 X 40 grid)
		double			step_size;	// grid step size
		long			data_size;	// total number of points

	public:
	
		// default ctor allocates a 1D grid starting at 0.0 to 1.0 
		// with 1.0 stepsize (2 elements at least needed for MOE)
		RbtMOEGridShape();
		
		// parametrized ctor receives grid as the origin vector, the extent
		// vector and stepsize
		RbtMOEGridShape(vector<double>, vector<double>, double);

		const vector<double>&	GetOrigin()		{ return origin;	};
		const vector<double>&	GetExtents()	{ return extents;	};
		const vector<int>&		GetIExtents()	{ return i_extents;	};
		double					GetStepSize()	{ return step_size;	};
		long					GetDataSize()	{ return data_size;	};
};

typedef	vector<RbtMOEGridPoint> RbtMOEGridData;

// MOE grid data class for IO
//
// the MOE grid file format (in the present stage at 06/03/2003)
class RbtMOEGrid
{
		RbtMOEGridShape		myShape;
		RbtMOEGridData		myData;
		RbtString			stream_name;	// filename
		RbtCoord			min,max;		// common extent min and max coords
		RbtCoord			c_min,c_max;		// current extent min and max coords
	public:
		// default ctor
		RbtMOEGrid();
		// parameterised ctor for a single grid 
		// pass shape, data and outfilename
		RbtMOEGrid(RbtMOEGridShape&, RbtMOEGridData&, RbtString);
		// or pass parameters using these functions (when it is about 
		// to save multiple grids into a single file)
		void	SetShape(RbtMOEGridShape& a_shape)			{myShape		= a_shape;};
		void	SetData(RbtMOEGridData& a_data)				{myData			= a_data;};
		void	SetOutputFileName(RbtString& a_filename)	{stream_name	= a_filename;};
		// methods to calculate grid extents
		void	CalculateCommonExtents(vector<RbtString> strPrmFiles);	// for all
		void	GetDockingSiteExtents(RbtString& a_strPrmFile);			// for a single
		// get methods
		RbtCoord	GetMinExtents()	{return min;};
		RbtCoord	GetMaxExtents()	{return max;};
		
		// iofstream methods
		long		WriteGrid(std::_Ios_Openmode mode);	// returns the number of points written
		long		ReadGrid();				// ditto (fields read)
};

// the relevant SVL code in gridshow.svl
//
// basic double should be the IEEE 8 byte double MSB first in Linux
// on big-endian try swab with dd and/or cpio
// MOE grid format:
//
// Note: the internal grid format (subject to change without notice!)
// Header for the grid-file only:
//  GRID_SIGNATURE          8 bytes
// Each grid:
//  version number          	int4m
//  ndata (= length data)       int4m
//  ndim (= length size)        int4m
//  skip                int4m
//  ... (reserved for future use)   byte * skip
//  size (= app length shape)   int4m * ndim
//  shape               ieee8m * (add size)
//  data                ieee8m * ndata
//
// On reading, we check that ndata, ndim, and size are 24-bit non-negatives,
// but we don't check "ndata == (mul size)"
//
// snippets from '$MOE/lib/svl/run/gridshow.svl'
// const GRID_SIGNATURE = "SVLgrid\n";
// const VERSION = 100;
// const SKIP = 0;
//              
// Check file signature and return non-zero if the file of file name fn
// is a "grid file", i.e. a special file that contains exactly one grid
// and nothing else
//
//local function read_gridsign fn =
//    freadb [fn, 'char', length GRID_SIGNATURE] === GRID_SIGNATURE;

// Write grid to a file of given filename or filehandle
//
//local function writeb_fnum_grid [fn, [data, shape]]
//    local ndata = length data;
//    local ndim = length shape;
//    local size = app length shape;
//
//    fwriteb [fn, 'int4m', [VERSION, ndata, ndim, SKIP]];
//    fwriteb [fn, 'char', rep ["\0", SKIP]];
//    fwriteb [fn, 'int4m', size];
//    fwriteb [fn, 'IEEE8m', cat shape];
//    fwriteb [fn, 'IEEE8m', data];
//endfunction
//
//function grid_fwriteb [fn, grid]
//    if type fn === 'num' then writeb_fnum_grid [fn, grid]; return; endif
//    fn = fopenw fn;
//    writeb_fnum_grid [fn, grid];
//    fclose fn;
//endfunction
//
// Create/rewrite "grid-file" of given name that willcontain a given grid
//
//function grid_writefile [fn, grid]
//    fn = fopenw fn;
//    fwriteb [fn, 'char', GRID_SIGNATURE];
//    writeb_fnum_grid [fn, grid];
//    fclose fn;
//endfunction


//    fwriteb [fn, 'char', GRID_SIGNATURE];
//    fwriteb [fn, 'int4m', [VERSION, ndata, ndim, SKIP]];
//    fwriteb [fn, 'char', rep ["\0", SKIP]];
//    fwriteb [fn, 'int4m', size];
//    fwriteb [fn, 'IEEE8m', cat shape];
//    fwriteb [fn, 'IEEE8m', data];
