#pragma once
#include "ai.hh"



//cojugate vector with std::vector<double> data
inline void conjugate(std::vector< std::vector<double> > &x, size_t N)
{
	size_t i = 0;

	for (; i <= N - 4; i += 4)
	{
		x[i][0] = x[i][0];
		x[i][1] = -x[i][1];

		x[i+1][0] =  x[i+1][0];
		x[i+1][1] = -x[i+1][1];

		x[i + 2][0] = x[i + 2][0];
		x[i + 2][1] = -x[i + 2][1];

		x[i + 3][0] = x[i + 3][0];
		x[i + 3][1] = -x[i + 3][1];
	}
	for (; i < N; ++i)
	{
		x[i][0] = x[i][0];
		x[i][1] = -x[i][1];
	}
}
// Cooley-Tukey FFT

/**/
inline void fft(std::vector< std::vector<double> > &x, size_t N)
{
	// DFT
	size_t k = N, n;
	double thetaT = 3.14159265358979 / ((double) N);

	std::vector<double> phiT = { cos(thetaT), -sin(thetaT) };
	std::vector<double> T;
	while (k > 1)
	{
		n = k;
		k >>= 1;

		phiT = { phiT[0] * phiT[0] - phiT[1] * phiT[1], 2 * phiT[0] * phiT[1] };
		std::vector<double> T { 1. , 0. };

		for (size_t l = 0; l < k; l++)
		{
			for (size_t a = l; a < N; a += n)
			{
				size_t b = a + k;
				std::vector<double> t = { x[a][0] - x[b][0], x[a][1] - x[b][1] };
				x[a][0] += x[b][0];
				x[a][1] += x[b][1];
				x[b][0] = t[0] * T[0] - t[1] * T[1];
				x[b][1] = t[0] * T[1] + t[1] * T[0];

			}
			double tempT = T[0];
			T[0] = T[0]* phiT[0] - T[1] * phiT[1];
			T[1] = tempT* phiT[1] + T[1] * phiT[0];

		}
	}
	// Decimate
	size_t m = (size_t)log2(N);
	for (size_t a = 0; a < N; ++a)
	{
		size_t b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			std::vector<double> t { x[a][0] , x[a][1] };
			x[a][0] = x[b][0];
			x[a][1] = x[b][1];
			x[b][0] = t[0];
			x[b][1] = t[1];

		}
	}
}
// inverse fft
/**/
/**/
inline void ifft(std::vector< std::vector<double> >& x, size_t N)
{
	// conjugate the std::vector<double> numbers
	conjugate(x, N);
	// fft
	fft(x, N);
	// conjugate the std::vector<double> numbers again
	conjugate(x, N);
	// scale the numbers
	for (size_t i = 0; i < N; i++)
	{
		x[i][0] /= N;
		x[i][1] /= N;

	}
}
/**/
inline void conjugate2(std::vector< std::vector<double> > &vector){
    std::size_t length = vector.size();

    for(std::size_t i = 0; i < length; ++i){
        vector[i][1] = -vector[i][1];
    }
}

inline void fft2(std::vector< std::vector<double> > &vector){
    std::size_t k = vector.size();
    std::size_t j = 0;
    std::size_t n = 0;

    const double length = (double) k;

    double thetaT = 3.14159265358979 / length;
    double swap0 = 0.;
    double swap1 = 0.;
    double T0 = 1.;
    double T1 = 0.;
    double phiT0 = cos(thetaT);
    double phiT1 = -sin(thetaT);

    while(k > 1){
        n = k;

        k >>= 1;

        swap0 = phiT0;
        swap1 = phiT1;

        phiT0 = swap0 * swap0 - swap1 * swap1;
        phiT1 = 2. * swap0 * swap1;

        T0 = 1.;
        T1 = 0.;

        for(std::size_t l = 0; l < k; ++l)
        {
            for(std::size_t i = l; i < length; i += n)
            {
                j = i + k;

                swap0 = vector[i][0] - vector[j][0];
                swap1 = vector[i][1] - vector[j][1];

                vector[i][0] += vector[j][0];
                vector[i][1] += vector[j][1];

                vector[j][0] = swap0 * T0 - swap1 * T1;
                vector[j][1] = swap0 * T1 + swap1 * T0;
            }

            swap0 = T0;

            T0 = swap0 * phiT0 - T1 * phiT1;
            T1 = swap0 * phiT1 + T1 * phiT0;
        }
    }

    std::size_t m = (std::size_t) log2(length);

    for(std::size_t i = 0; i < length; ++i){
        j = i;

        j = (((j & 0xaaaaaaaa) >> 1) | ((j & 0x55555555) << 1));
        j = (((j & 0xcccccccc) >> 2) | ((j & 0x33333333) << 2));
        j = (((j & 0xf0f0f0f0) >> 4) | ((j & 0x0f0f0f0f) << 4));
        j = (((j & 0xff00ff00) >> 8) | ((j & 0x00ff00ff) << 8));
        j = ((j >> 16) | (j << 16)) >> (32 - m);

        if(j > i){
            swap0 = vector[i][0];
            swap1 = vector[i][1];

            vector[i][0] = vector[j][0];
            vector[i][1] = vector[j][1];

            vector[j][0] = swap0;
            vector[j][1] = swap1;
        }
    }
}

inline void ifft2(std::vector< std::vector<double> > &vector){
    const double length = (double) vector.size();

    conjugate2(vector);

    fft2(vector);
		conjugate2(vector);
		
    for(std::size_t i = 0; i < length; ++i){
        vector[i][0] /= length;
				vector[i][1] /= length;
    }
}
