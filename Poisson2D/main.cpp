#include <vector>
#include "particle_struct.h"
//#include "a3r.h"
#include "ai.hh"
#include "particles_ut.h"
#include "box_struct.h"
#include "Fourier_tools.hh"

int main()
{
	std::vector <Particle> partilces;
	std::vector <BOX> boxes;

	double H = 1.;

	const size_t dim = 4;
	size_t number_particles = 1;
	double mass = 10;
	double h = 1 / double(dim - 1);


	std::vector<
		std::vector<
			std::vector< 
				std::vector<double>>>> rho;
	rho.resize(dim);
	for (size_t i = 0; i < dim; ++i) {
		rho[i].resize(dim);
		for (size_t j = 0; j < dim; j++)
		{
			rho[i][j].resize(dim);
			for (size_t  k = 0; k <dim; k++)
			{
				rho[i][j][k].resize(2);
			}
		}
	}

	//std::complex <double> rho[dim][dim];
	
	
	
	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = 0; j < dim; ++j)
			for (size_t k = 0; k < dim; k++)
		{
			if (j == dim / 2 && i == dim / 2 && k==dim/2)
				rho[i][j][k][0] = mass / (h*h);
			//rho[i][j] = mass / (h*h);
			
		}
	}
	std::cout << "density field" << std::endl;
	for (size_t i = 0; i < dim; i++)
		for (size_t j = 0; j < dim; ++j)
			for (size_t k = 0; k < dim; ++k)
		{
			std::cout << "(" << i <<" , "<< j << " , " << k << " )" << " dens = " << " ( " << rho[i][j][k][0] << " , "<< rho[i][j][k][1]
				<< " )" << std::endl;
		}

	

	//calculate FFT for each layer
	for (size_t i = 0; i < dim; ++i)
	{
		std::vector <
			std::vector<double>> f; //for rows and colomns
		f.resize(dim);
		for (size_t l = 0; l < dim; l++)
		{
			f[l].resize(2);
		}

		//FFT rows of density
		for (size_t j = 0; j < dim; ++j) {
			for (size_t k = 0; k < dim; ++k) {
				f[k][0] = rho[i][j][k][0];
				f[k][1] = rho[i][j][k][1];
				
			}
			fft(f, dim);
			for (size_t k = 0; k < dim; ++k) {
				rho[i][j][k][0] = f[k][0];
				rho[i][j][k][1] = f[k][1];
				
			}
		}
		//FFT coloumns of density
		for (size_t k = 0; k < dim; ++k) {
			for (size_t j = 0; j < dim; ++j) {
				f[j][0] = rho[i][j][k][0];
				f[j][1] = rho[i][j][k][1];
				
			}
			fft(f, dim);
			for (size_t j = 0; j < dim; ++j) {
				rho[i][j][k][0] = f[j][0];
				rho[i][j][k][1] = f[j][1];
			}
		}

	}
	std::cout << "FFT for layers ....DONE" << std::endl;
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
		//std::complex<double> z(0.0, 1.0);
		double pi = 3.14159265358979;
		std::vector<double> NUL = { 0., 0. };
		//std::complex<double> W = std::exp(2.0 * pi * z / double(dim));
		std::vector<double> W = { cos(2.0 * pi / double(dim)) , sin(2.0 * pi / double(dim)) };
		//std::complex<double> Wm = 1.0, Wn = 1.0 ,  Wk = 1.0;
		std::vector<double> Wm = { 1.0, 0. }, Wn = { 1.0, 0. }, Wk = { 1.0, 0. };
		for (size_t k = 0; k < dim; k++) {
			for (size_t m = 0; m < dim; m++) {
				for (size_t n = 0; n < dim; n++) {

					//std::complex<double> denom = 6.0;
					std::vector<double> denom = { 6.0 ,0.};
					
					
					std::vector < double> i_Wk = { Wk[0] / (Wk[0] * Wk[0] + Wk[1]*Wk[1]), -(Wk[1] / (Wk[0] * Wk[0] + Wk[1] * Wk[1])) };
					std::vector < double> i_Wm = { Wm[0] / (Wm[0] * Wm[0] + Wm[1] * Wm[1]), -(Wm[1] / (Wm[0] * Wm[0] + Wm[1] * Wm[1])) };
					std::vector < double> i_Wn = { Wn[0] / (Wn[0] * Wn[0] + Wn[1] * Wn[1]), -(Wn[1] / (Wn[0] * Wn[0] + Wn[1] * Wn[1])) };
					
					//denom -= Wk + 1.0 / Wk + Wm + 1.0 / Wm + Wn + 1.0 / Wn;
					denom[0] -=  ( Wk[0] + i_Wk[0] + Wm[0] + i_Wm[0] + Wn[0] + i_Wn[0]);
					denom[1] -=  ( Wk[1] + i_Wk[1] + Wm[1] + i_Wm[1] + Wn[1] + i_Wn[1]);
					if (denom != NUL) {

						//std::vector<double> hh_denom = { (h*h*denom[0]) / (denom[0] * denom[0] + denom[1] * denom[1]),-(h * h*denom[1]) / (denom[0] * denom[0] + denom[1] * denom[1]) };

						double q = rho[k][m][n][0];
						rho[k][m][n][0] = h * h*(rho[k][m][n][0] * denom[0] + rho[k][m][n][1] * denom[1]) / (denom[0] * denom[0]+denom[1]*denom[1]);
						rho[k][m][n][1] = h*h*(rho[k][m][n][1] * denom[0] - q * denom[1])/ (denom[0] * denom[0] + denom[1] * denom[1]);
						//rho[k][m][n] = rho[k][m][n] * h * h / denom;
					}
					//Wn = Wn * W;
					double b = Wn[0];
					Wn[0] = Wn[0] * W[0] - Wn[1] * W[1];
					Wn[1] = Wn[1] * W[0] + b * W[1];
				}
				//Wm = Wm * W;
				double r = Wm[0];
				Wm[0] = Wm[0] * W[0]- Wm[1]*W[1];
				Wm[1] = Wm[1] * W[0] + r * W[1];
			}
			//Wk *= W;
			double a = Wk[0];
			Wk[0] =Wk[0] * W[0] - Wk[1]*W[1];
			Wk[1] =Wk[1] * W[0] + a*W[1];
		}
		std::cout << "Solve Poisson eq. in Fourier space" << std::endl;
		for (size_t i = 0; i < dim; i++)
		{
			std::vector<std::vector<double>> f;
			f.resize(dim);
			for (size_t y = 0; y < dim; y++)
			{
				f[y].resize(2);
			}
			// inverse FFT rows of rho 
			for (size_t j = 0; j < dim; j++) {
				for (size_t k = 0; k < dim; k++) {
					f[k][0] = rho[i][j][k][0];
					f[k][1] = rho[i][j][k][1];
				}
				ifft(f, dim);
				for (size_t k = 0; k < dim; k++)
				{
					rho[i][j][k][0] = f[k][0];
					rho[i][j][k][1] = f[k][1];
				}
			}

			// inverse FFT columns of rho 
			for (size_t k = 0; k < dim; k++) {
				for (size_t j = 0; j < dim; j++)
				{
					f[j][0] = rho[i][j][k][0];
					f[j][1] = rho[i][j][k][1];
				}
				ifft(f, dim);
				for (size_t j = 0; j < dim; j++)
				{
					rho[i][j][k][0] = f[j][0];
					rho[i][j][k][1] = f[j][1];
				}
			}
		}
		std::cout << "Solving ... DONE" << std::endl;
	//ai::saveMatrix("./rho", rho);
	//output
	std::cout << "Potential field  =" << std::endl;
	for (size_t i = 0; i < dim; i++)
		for (size_t j = 0; j < dim; ++j)
			for (size_t k = 0; k < dim; ++k)
		{
			std::cout << "(" <<i<<" , "<< j << " , " << k << " )" << " phi = " << " ( " << rho[i][j][k][0] << " , "<<rho[i][j][k][1]
				<< " )" << std::endl;
		}
	return 0;

}