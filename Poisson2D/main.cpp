#include <vector>
#include <complex>
#include "particle_struct.h"
//#include "a3r.h"
#include "ai.hh"
#include "particles_ut.h"
#include "box_struct.h"
#include "Fourier_complex.h"

int main()
{
	std::vector <Particle> partilces;
	std::vector <BOX> boxes;

	double H = 1.;

	const size_t dim = 64;
	size_t number_particles = 1;
	double mass = 10;
	double h = 1 / double(dim - 1);


	std::vector<
		std::vector<
			std::complex <double>>> rho;
	rho.resize(dim);
	for (size_t i = 0; i < dim; ++i) {
		rho[i].resize(dim);
		
	}

	//std::complex <double> rho[dim][dim];
	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = 0; j < dim; ++j)
		{
			if (j == dim/2 && i == dim/2)
				rho[i][j] = mass / (h*h);
				//rho[i][j] = mass / (h*h);
			else
				//rho[i][j] = 0.;
				rho[i][j] = 0.;
		}
	}
	std::cout << "density field" << std::endl;
	for (size_t j = 0; j < dim; ++j)
		for (size_t k = 0; k < dim; ++k)
		{
			std::cout << "(" << j << " , " << k << " )" << " dens = " << " ( " << rho[j][k] << " , "
				 << " )" << std::endl;
		}

	/*std::cout << "density field" << std::endl;
	for (size_t j = 0; j < dim; ++j)
		for (size_t k = 0; k < dim; ++k)
		{
			std::cout << "(" << j << " , " << k << " )" << " dens = " << " ( " << rho[j][k].real << " , "
				<< rho[j][k].imag << " )" << std::endl;
		}*/
	//std::vector <
	//	std::vector<double>> f(dim); //for rows and colomns
	//for (size_t k = 0; k < dim; ++k) { f[k].resize(2); } //resize f as complex

	std::vector<std::complex<double>> f;
	f.resize(dim);
	//FFT for rows
	for (size_t j = 0; j < dim; j++) {
		for (size_t k = 0; k < dim; k++) {
			f[k] = rho[j][k];
			//f[k][1] = rho[j][k][1];
			//f[k] = rho[j][k];
		}
		fft(f, dim);
		for (size_t k = 0; k < dim; k++) {
			rho[j][k] = f[k];
			//rho[j][k][1] = f[k][1];
		}
	}
	// FFT columns of rho 
	for (size_t k = 0; k < dim; k++) {
		for (size_t j = 0; j < dim; j++)
		{
			f[j] = rho[j][k];
			//f[j][1] = rho[j][k][1];
		}
				fft(f , dim); 
				for (size_t j = 0; j < dim; j++) {
					rho[j][k] = f[j];
					//rho[j][k][1] = f[j][1];
				}
	}


	/*double pi = 3.14159265358979;
	std::vector<double> W = { 0. , std::exp(2.0 * pi / double(dim)) };
	std::vector<double> Wn = { 1. , 0. };
	std::vector<double> Wm = { 1. , 0. };
	for (size_t m = 0; m < dim; m++) {
		for (size_t n = 0; n < dim; n++) {

			std::vector<double> d = { 4.0, 0.0 };
			d[0] -= Wm[0] + (Wm[0] / (Wm[0] * Wm[0] + Wm[1] * Wm[1])) + Wn[0] +  (Wn[0]/(Wn[0] * Wn[0] + Wn[1] * Wn[1]));
					
			d[1] -= Wm[1] -  (Wm[1]/(Wm[0] * Wm[0] + Wm[1] * Wm[1])) + Wn[1] - (Wn[1]/(Wn[0] * Wn[0] + Wn[1] * Wn[1]));
					
			if (d[0] != 0.0 && d[1] != 0.0)
					{
						double q = rho[m][n][0];
						
						rho[m][n][0] = rho[m][n][0] * (h * h) * (d[0]/(d[0]*d[0]+d[1]*d[1])) +
							rho[m][n][1]*(h*h)* (d[1]/(d[0] * d[0] + d[1] * d[1]));
						
						
						rho[m][n][1] = rho[m][n][1] * (h * h) * (d[0] / (d[0] * d[0] + d[1] * d[1])) -
							q * (h*h)* (d[1] / (d[0] * d[0] + d[1] * d[1]));
					}
						
					double y = Wn[0];
					
					Wn[0] = Wn[0]*W[0] - Wn[1]*W[1]; 
					Wn[1] = Wn[1]*W[0] + y*W[1]; 
			} 
				double r = Wm[0];
			Wm[0] = Wm[0] * W[0] - Wm[1] * W[1];
			Wm[1] = Wm[1] * W[0] + r * W[1];
	}*/
	////Solve  equation in Fourier space
	//double pi = 3.14159265358979;
	//std::vector<double> W = { 0. , std::exp(2.0 * pi / double(dim)) };
	////std::vector<double> Wl = { 1. , 0. };
	//std::vector<double> Wn = { 1. , 0. };
	//std::vector<double> Wm = { 1. , 0. };
	////for (size_t l = 0; l < 0; ++l) { // for 3D
	//for (size_t m = 0; m < dim; ++m) {
	//	for (size_t n = 0; n < dim; ++n) {
	//		//std::vector<double> d = { 6., 0. }; // for 3D
	//		std::vector<double> d = { 4., 0. }; // for slices

	//											//d[0] -= Wl[0] + 1. / Wl[0] + Wn[0] + 1. / Wn[0] + Wm[0] + 1. / Wm[0];

	//											/*d[0] =d[0] -( Wl[0] + (Wl[0]/(Wl[0]* Wl[0] + Wl[1]* Wl[1])) + Wn[0] + (Wn[0]/(Wn[0] * Wn[0] + Wn[1] * Wn[1])) +
	//											Wm[0] + (Wm[0] / (Wm[0] * Wm[0] + Wm[1] * Wm[1])));  for 3D /**/
	//		d[0] -= Wn[0] + (Wn[0] / (Wn[0] * Wn[0] + Wn[1] * Wn[1])) +
	//			Wm[0] + (Wm[0] / (Wm[0] * Wm[0] + Wm[1] * Wm[1]));
	//		/*
	//		d[1] =d[0] - (Wl[1] - (Wl[1] / (Wl[0] * Wl[0] + Wl[1] * Wl[1])) + Wn[1] -(Wn[1] / (Wn[0] * Wn[0] + Wn[1] * Wn[1])) +
	//		Wm[1] - (Wm[1] / (Wm[0] * Wm[0] + Wm[1] * Wm[1]))); for 3D/**/
	//		d[1] -= Wn[1] - (Wn[1] / (Wn[0] * Wn[0] + Wn[1] * Wn[1])) +
	//			Wm[1] - (Wm[1] / (Wm[0] * Wm[0] + Wm[1] * Wm[1]));

	//		if (d[0] != 0. && d[1] != 0.) {
	//			double q = -rho[m][n][0];
	//			//density[l][n][m][0] = density[l][n][m][0]*((h*h*h) / d[0]) - density[l][n][m][1] * ((h*h*h) / d[1]);
	//			rho[m][n][0] = rho[m][n][0] * ((h*h)*d[0] / (d[0] * d[0] + d[1] * d[1])) -
	//				rho[m][n][1] * ((h*h)* d[1] / (d[0] * d[0] + d[1] * d[1]));

	//			rho[m][n][1] = q * ((h*h) * d[1] / (d[0] * d[0] + d[1] * d[1])) +
	//				rho[m][n][1] * ((h*h) * d[0] / (d[0] * d[0] + d[1] * d[1]));
	//		}

	//		double s = Wn[0];
	//		Wn[0] = Wn[0] * W[0] - Wn[1] * W[1];
	//		Wn[1] = s * W[1] + Wn[1] * W[0];
	//	}
	//	double b = Wm[0];
	//	Wm[0] = Wm[0] * W[0] - Wm[1] * W[1];
	//	Wm[1] = b * W[1] + Wm[1] * W[0];
	//}

	// solve equation in Fourier space 
	std::complex<double> i(0.0, 1.0); 
		double pi = 4 * std::atan(1.0); 
		std::complex<double> W = std::exp(2.0 * pi * i / double(dim)); 
		std::complex<double> Wm = 1.0, Wn = 1.0; 
		for (int m = 0; m < dim; m++) {
			
				for (int n = 0; n < dim; n++) {
					
						std::complex<double> denom = 4.0; 
						denom -= Wm + 1.0 / Wm + Wn + 1.0 / Wn; 
						if (denom != 0.0) 
							rho[m][n] *= h * h / denom; 
							Wn *= W; 
				} 
					Wm *= W; 
		}
	// inverse FFT rows of rho 
	for (size_t j = 0; j < dim; j++) {
		for (size_t k = 0; k < dim; k++) {
			f[k] = rho[j][k];
			//f[k][1] = rho[j][k][1];
		}
				ifft(f , dim); 
				for (size_t k = 0; k < dim; k++)
				{
					rho[j][k] = f[k];
					//rho[j][k][1] = f[k][1];
				}
	}

	// inverse FFT columns of rho 
	for (size_t k = 0; k < dim; k++) {
		for (size_t j = 0; j < dim; j++)
		{
			f[j]= rho[j][k];
			//f[j][1] = rho[j][k][1];
		}
				ifft(f, dim); 
				for (size_t j = 0; j < dim; j++)
				{
					rho[j][k] = f[j];
					//rho[j][k][1] = f[j][1];
				}
	}
	//ai::saveMatrix("./rho", rho);
	//output
	std::cout << "Potential field  =" << std::endl;
	
		for (size_t j = 0; j < dim; ++j)
			for (size_t k = 0; k < dim; ++k)
			{
				std::cout << "(" << j << " , " << k << " )" << " phi = " << " ( " << rho[j][k] << " , "
					<< " )" << std::endl;
			}
	return 0;

}