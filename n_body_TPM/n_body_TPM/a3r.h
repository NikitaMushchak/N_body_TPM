//-----------------------------------------------------------------------------------
//
//	Saving and loading particle coordinates in A3R format
//
//	Anton M. Krivtsov
//
//	11.04.2001
//
//-----------------------------------------------------------------------------------

#ifndef ___a3r_h___
#define ___a3r_h___

//#define __USE_STD_IOSTREAM

#include "vect3d.h"
#include "util.h"
#include "a3r.h"

//-----------------------------------------------------------------------------------

struct A3R_HEADER
{
	A3R_HEADER();		// Initialization of the structure members

	char file_type[4];	// "a3r"
	int count;			// number of particles
	int data_start;	// address of the start of the particles data
	char version[10];	// version of a3r format
	double r;			// particle radius
	int count_1;		// reserved for the future use;
};

//-------------------------------------------------------------------------------------------

int Save_A3R(const char* file_name, Vect3D* start, float* colors, int n, double r);

Vect3D* Load_A3R(const char* file_name, int& n, double& r);	

//-------------------------------------------------------------------------------------------

#endif //___a3r_h___

