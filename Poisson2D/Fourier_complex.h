#pragma once
#pragma once
#include "ai.hh"
#include <complex>

//generate circulant
inline void createCirculantMatrix(
	std::vector< std::vector<double> > &matrix,
	std::vector< double> &vector
) {
	size_t length = vector.size();

	size_t displacement = 0;

	matrix.resize(length);

	for (size_t i = 0; i < length; ++i) {
		matrix[i].resize(length);

		size_t k = length - 1;
		for (size_t j = displacement; j < length; ++j) {
			matrix[i][k] = vector[j];

			--k;
		}
		for (size_t j = 0; j < displacement; ++j) {
			matrix[i][k] = vector[j];

			--k;
		}
		++displacement;
	}
}
///wise element miltiplication vector-vector
inline void getWiseElement(std::vector< std::vector<double> > &a, std::vector< std::vector<double> > &b, size_t m)
//inline void getWiseElement(std::complex<double> &a, std::vector< std::vector<double> > &b , size_t m)
{
	for (size_t i = 0; i < m; ++i)
	{
		double f = b[i][0];
		b[i][0] = a[i][0] * b[i][0] - a[i][1] * b[i][1];
		b[i][1] = a[i][0] * b[i][1] + a[i][1] * f;
	}
}

//cojugate vector with std::vector<double> data
//inline void conjugate(std::vector< std::vector<double> > &x, size_t N)
inline void conjugate(std::vector<std::complex<double>> &x, size_t N)
{
	/*size_t i = 0;

	for (; i <= N - 4; i += 4)
	{
		x[i][0] = x[i][0];
		x[i][1] = -x[i][1];

		x[i + 1][0] = x[i + 1][0];
		x[i + 1][1] = -x[i + 1][1];

		x[i + 2][0] = x[i + 2][0];
		x[i + 2][1] = -x[i + 2][1];

		x[i + 3][0] = x[i + 3][0];
		x[i + 3][1] = -x[i + 3][1];
	}
	for (; i < N; ++i)
	{
		x[i][0] = x[i][0];
		x[i][1] = -x[i][1];
	}/**/
	for (size_t i = 0; i < N; ++i)
	{
		//x[i] = 0;
		//std::conj(a);
		 std::conj(x[i]);

	}/**/
}
// Cooley-Tukey FFT

/**/
inline void fft(std::vector< std::complex<double> > &x, size_t N)
{
	// DFT
	size_t k = N, n;
	double thetaT = 3.14159265358979323846264338328L / N;
	std::complex<double> phiT = std::complex<double>(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (size_t l = 0; l < k; l++)
		{
			for (size_t a = l; a < N; a += n)
			{
				size_t b = a + k;
				std::complex<double> t = x[a] - x[b];
				x[a] += x[b];
				x[b] = t * T;
			}
			T *= phiT;
		}
	}
	// Decimate
	size_t m = (size_t)log2(N);
	for (size_t a = 0; a < N; a++)
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
			std::complex<double> t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
}
// inverse fft
/**/
/**/
inline void ifft(std::vector< std::complex<double> >& x, size_t N)
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
		x[i] /= N;
		

	}
}
/**/
//inline void SaveMartix_XYZ_Complex(std::vector<std::complex<double>> & x, size_t N)
//{
//	const char thr = Bond3Dv2::GapM();
//	int number = 0;
//	char str[100];
//	sprintf(str, "output/crack%07u.xyz", n_step);
//	int count = 0;
//	std::vector<int> numPartGap;
//
//	for (int i = 0; i < n; i++)
//	{
//		if (as[i].N_Gap() < thr) { continue; }
//		number++;
//		numPartGap.push_back(i);
//	}
//	/*Bond3Dv2* pB = bonds.GetData();
//	Bond3Dv2* pE = pB + bonds.GetCount();*/
//	//double d_V = 0;
//	//for (; pB<pE; pB++)           // цикл по всем связям
//	//{
//	//	if (pB->Shifted()) { continue; }
//	//	//if (pB->Gap() < Bond3Dv2::GapM()) { continue; }
//	//	if (pB->Gap() == thr)
//	//	{
//	//		int i = pB->First();
//	//		int j = pB->Second();
//	//		double L = pB->dR(as, area).Abs();         //  значение dR сохраняется в переменной (*b_i)
//	//		d_V += abs(L - pB->Length())*d_s;
//	//	}
//	//}
//	//	
//	FILE* f = fopen(str, "w");
//	fprintf(f, "%i\n", number);
//	//fprintf(f, "Volume  \t%f\n", dV);
//	int k = 1;
//	for (int it = 0; it<numPartGap.size(); it++)
//	{
//		int i = numPartGap.at(it);
//		fprintf(f, "%i %f %f %f %f %\n", k, as[i].R().x, as[i].R().y, as[i].R().z, 1.0);
//	}
//	fclose(f);
//}
