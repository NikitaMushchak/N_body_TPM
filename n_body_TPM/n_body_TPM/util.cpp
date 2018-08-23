#include "uitil.h"

//------------------------------------------------------------------------------------//

int is_inverse_byte_order()
{
	short x = 1;
	char *p = reinterpret_cast<char*> (&x);
	return 0 != (int)(*p);
}


//------------------------------------------------------------------------------------//
/*
“Minimal” random number generator of Park and Miller.
Returns a uniform random deviate between 0.0 and 1.0.
Set or reset idum to any integer value (except the unlikely value MASK)
to initialize the sequence;
idum must not be altered between calls for successive deviates in
a sequence.
*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

static long idum = MASK;

//------------------------------------------------------------------------------------//

void seed_uni_rand(const long zidnum)
{
	idum = zidnum ^ MASK;
}

//------------------------------------------------------------------------------------//

float uni_rand()
{
	long k = idum / IQ;
	// Compute idum=(IA*idum) % IM without overflows
	// by Schrage’s method
	idum = IA * (idum - k * IQ) - IR * k;
	if (idum < 0) idum += IM;
	// Convert idum to a floating result.
	float ans = (float)AM * idum;
	return ans;
}

//------------------------------------------------------------------------------------//

#undef IA 
#undef IM 
#undef AM 
#undef IQ 
#undef IR 
#undef MASK 

//------------------------------------------------------------------------------------//
