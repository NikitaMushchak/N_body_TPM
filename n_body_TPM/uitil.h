#pragma once
//-----------------------------------------------------------------------------------
//
//	Mathematical utilities
//
//	Anton M. Krivtsov
//
//	11.04.2001
//
//-----------------------------------------------------------------------------------

#ifndef ___UTIL_H___
#define ___UTIL_H___
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <sstream>

using namespace std;
//-----------------------------------------------------------------------------------

static const double PI = 3.1415926535897932;
static const double DEG = PI / 180;

//-----------------------------------------------------------------------------------

int is_inverse_byte_order();


//-----------------------------------------------------------------------------------
//	RAND_MAX = 32767 = 2^15 - 1
//-----------------------------------------------------------------------------------

inline double rand(const double max)			//	Random double value from 0 to max
{
	return max * (double)rand() / RAND_MAX;
}

//-----------------------------------------------------------------------------------

inline double rand(const double min, const double max)	//	Random double value from min to max
{
	return min + (max - min) * (double)rand() / RAND_MAX;
}

//-----------------------------------------------------------------------------------

inline int rand(int max)		//	Random integer value from 0 to (max - 1)
{
	return max * rand() >> 15;	//	x >> 15 = x / ((int)RAND_MAX + 1)
}

//-----------------------------------------------------------------------------------

//template <class T> inline void swap(T& x1, T& x2) { T x = x1; x1 = x2; x2 = x; }

//-----------------------------------------------------------------------------------

template <class T> inline T sqr(const T x) { return x * x; }

//-----------------------------------------------------------------------------------

template <class T> string str(const T& x)
{
	stringstream buf;
	buf << x;
	return buf.str();
}

//-----------------------------------------------------------------------------------

template <class T> string strs(const T& x)
{
	stringstream buf;
	buf << x << " ";
	return buf.str();
}

//-----------------------------------------------------------------------------------

template <class type> type max_num(const type a, const type b) { return a>b ? a : b; }
template <class type> type min_num(const type a, const type b) { return a<b ? a : b; }

template <class type> type max_num(const type a, const type b, const type c) { return max_num(a, max_num(b, c)); }
template <class type> type min_num(const type a, const type b, const type c) { return min_num(a, min_num(b, c)); }

template <class type> type min_max_num(const type x, const type a, const type b) { return max_num(min_num(x, b), a); }
//-----------------------------------------------------------------------------------

inline void swap_double(double &x)
{
	char *p = reinterpret_cast<char*>(&x);
	char q[8] = { p[7], p[6], p[5], p[4], p[3], p[2], p[1], p[0] };
	memcpy((void*)(&x), (void*)q, 8);
}

//------------------------------------------------------------------------------------//

void seed_uni_rand(const long zidnum);

float uni_rand();

inline float uni_rand(const float a, const float b)
{
	return a + (b - a) * uni_rand();
}

//------------------------------------------------------------------------------------//
/*

REAL FUNCTION GSU2R(ISEED)
INTEGER ISEED,D2P32M
DOUBLE PRECISION Z,D2P31M,D2PN31,DMOD,DFLOAT
DATA  D2PN31/4.656612873077393D-10/,D2P31M/2147483647.D0/,D2P32M/16807/
Z=DFLOAT(ISEED)
Z=DMOD(D2P32M*Z,D2P31M)
GSU2R=Z*D2PN31
ISEED=Z
RETURN
END

*/
//------------------------------------------------------------------------------------//

#ifdef _WIN32
inline int ut_isnan(double x)
{
	return _isnan(x);
}
#else
#ifdef unix
inline int ut_isnan(double x)
{
	return isnan(x);
}
#else
#error "Unknown operation system. Function ut_isnan() defined only for \"unix\" and \"WIN32\""
#endif
#endif
//------------------------------------------------------------------------------------//
#ifdef _WIN32
inline int ut_isinf(double x)
{
	return !_finite(x);
}
#else
#ifdef unix
inline int ut_isinf(double x)
{
	return isinf(x);
}
#else
#error "Unknown operation system. Function ut_isinf() defined only for \"unix\" and \"WIN32\""
#endif
#endif

//------------------------------------------------------------------------------------//
#endif //___UTIL_H___
