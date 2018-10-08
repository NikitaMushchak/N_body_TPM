#include "CaclModule.h"
#include <vector>
#include "Fourier_tools.hh"
#include "particle_struct.h"
#include "box_struct.h"



 void CaclDensity(std::vector<
					std::vector<double>>& Particles , 
				 std::vector<
					std::vector <
						std::vector<
							std::vector<double>>>>& density,
				 float mass,
  			  	 double H, 
				 size_t dim)
{
	//std::fill(density.begin(), density.end(), 0.);
	//for ()
	size_t x, y, z = 0;
	double dx, dy, dz, tx, ty, tz;
	
	
	for (size_t i = 0; i < Particles.size(); ++i)
	{
		//calculating indexes
		dx = Particles[i][0] - std::floor(Particles[i][0] / H);
		dy = Particles[i][1] - std::floor(Particles[i][1] / H);
		dz = Particles[i][2] - std::floor(Particles[i][2] / H);
		tx = 1. - dx;
		ty = 1. - dy;
		tz = 1. - dz;
		if (std::floor(Particles[i][0] / H) < dim &&
			std::floor(Particles[i][1] / H)  < dim &&
			std::floor(Particles[i][2] / H)  < dim) {
			x = std::floor(Particles[i][0] / H);
			y = std::floor(Particles[i][1] / H);
			z = std::floor(Particles[i][2] / H);
			
			density[x][y][z][0] += mass * tx *ty*tz;

			if (y + 1 >= dim) { y = y % (dim - 1); density[x][y][z][0] += mass * tx *dy*tz; }
			else { density[x][y + 1][z][0] += mass * tx *dy*tz; }
			if (z + 1 >= dim) { z = z % (dim - 1); density[x][y][z][0] += mass * tx *ty*dz; }
			else { density[x][y][z + 1][0] += mass * tx *ty*dz; }
			if (y + 1 >= dim && z + 1 >= dim) { z = z % (dim - 1); y = y % (dim - 1); density[x][y][z][0] += mass * tx *dy*dz; }
			else { density[x][y + 1][z + 1][0] += mass * tx *dy*dz; }
			if (x + 1 >= dim) { x = x % (dim - 1); density[x][y][z][0] += mass * dx *ty*tz; }
			else { density[x + 1][y][z][0] += mass * dx *ty*tz; }
			if (x + 1 >= dim && y + 1 >= dim) { x = x % (dim - 1); y = y % (dim - 1); density[x][y][z][0] += mass * dx *dy*tz; }
			else { density[x + 1][y + 1][z][0] += mass * dx *dy*tz; }
			if (x + 1 >= dim && z + 1 >= dim) { x = x % (dim - 1); z = z % (dim - 1); density[x][y][z][0] += mass * dx *ty*dz; }
			else { density[x + 1][y][z + 1][0] += mass * dx *ty*dz; }
			if (x + 1 >= dim && y + 1 >= dim && z + 1 >= dim) { x = x % (dim - 1); y = y % (dim - 1); z = z % (dim - 1); density[x][y][z][0] += mass * dx * dy * dz; }
			else { density[x + 1][y + 1][z + 1][0] += mass * dx * dy * dz; }

		}

	}
}


 void CalcPotential(std::vector<
					std::vector<
						std::vector<
							std::vector<double>>>>& rho,
							size_t dim )
{
	double h = 1 / double(dim - 1);
	
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

	//}
	//std::cout << "FFT for layers ....DONE" << std::endl;

	// solve equation in Fourier space 
	//std::complex<double> z(0.0, 1.0);
	double pi = 3.14159265358979;
	std::vector<double> NUL = { 0., 0. };
	std::vector<double> W = { cos(2.0 * pi / double(dim)) , sin(2.0 * pi / double(dim)) };
	std::vector<double> Wm = { 1.0, 0. }, Wn = { 1.0, 0. };//, Wk = { 1.0, 0. };
	//for (size_t k = 0; k < dim; k++) {
		for (size_t m = 0; m < dim; m++) {
			for (size_t n = 0; n < dim; n++) {
				std::vector<double> denom = { 4.0 ,0. };

				//std::vector < double> i_Wk = { Wk[0] / (Wk[0] * Wk[0] + Wk[1] * Wk[1]), -(Wk[1] / (Wk[0] * Wk[0] + Wk[1] * Wk[1])) };
				std::vector < double> i_Wm = { Wm[0] / (Wm[0] * Wm[0] + Wm[1] * Wm[1]), -(Wm[1] / (Wm[0] * Wm[0] + Wm[1] * Wm[1])) };
				std::vector < double> i_Wn = { Wn[0] / (Wn[0] * Wn[0] + Wn[1] * Wn[1]), -(Wn[1] / (Wn[0] * Wn[0] + Wn[1] * Wn[1])) };

				denom[0] -= (/*Wk[0] + i_Wk[0] +*/ Wm[0] + i_Wm[0] + Wn[0] + i_Wn[0]);
				denom[1] -= (/*Wk[1] + i_Wk[1] +*/ Wm[1] + i_Wm[1] + Wn[1] + i_Wn[1]);
				if (denom != NUL) {

					

					double q = rho[i][m][n][0];
					rho[i][m][n][0] = 
						//h * h*
						(rho[i][m][n][0] * denom[0] + rho[i][m][n][1] * denom[1]) / (denom[0] * denom[0] + denom[1] * denom[1]);
					rho[i][m][n][1] =
						//h * h*
						(rho[i][m][n][1] * denom[0] - q * denom[1]) / (denom[0] * denom[0] + denom[1] * denom[1]);
				
				}
				double b = Wn[0];
				Wn[0] = Wn[0] * W[0] - Wn[1] * W[1];
				Wn[1] = Wn[1] * W[0] + b * W[1];
			}
			double r = Wm[0];
			Wm[0] = Wm[0] * W[0] - Wm[1] * W[1];
			Wm[1] = Wm[1] * W[0] + r * W[1];
		}
	/*	double a = Wk[0];
		Wk[0] = Wk[0] * W[0] - Wk[1] * W[1];
		Wk[1] = Wk[1] * W[0] + a * W[1];
	}*/
	std::cout << "Solve Poisson eq. in Fourier space" << std::endl;
	/*for (size_t i = 0; i < dim; i++)
	{*/
		//std::vector<std::vector<double>> f;
		//f.resize(dim);
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
}

 void GetAccel(std::vector<std::vector<double>>& Particles, 
						std::vector<std::vector<std::vector<std::vector<double >>>>& density,
							std::vector<std::vector<double>>& a, double H)
{
	std::vector <
		std::vector <
		std::vector <
		std::vector<double>>>> g;
	size_t dim = density.size();
	g.resize(dim);
	for (size_t i = 0; i < dim; ++i)
	{
		g[i].resize(dim);
		for (size_t j = 0; j < dim; ++j)
		{
			g[i][j].resize(dim);
			for (size_t k = 0; k < dim; ++k)
			{
				g[i][j][k].resize(3);
			}
		}
	}
	std::cout << "making acceleretion fields ....";
	//acceleration field
	for (size_t i = 1; i < dim - 1; ++i)
	{
		for (size_t j = 1; j < dim - 1; ++j)
		{
			for (size_t k = 1; k < dim - 1; ++k)
			{
				g[i][j][k][0] = -0.5*(density[i + 1][j][k][0] - density[i - 1][j][k][0]);
				g[i][j][k][1] = -0.5*(density[i][j + 1][k][0] - density[i][j - 1][k][0]);
				g[i][j][k][2] = -0.5*(density[i][j][k + 1][0] - density[i][j][k - 1][0]);
			}
		}
	}
	std::cout << "Done." << std::endl;
	//accleration of particles
	/*std::vector<
		std::vector<double>> a;*/
	a.resize(Particles.size());
	for (size_t i = 0; i < Particles.size(); ++i)
	{
		a[i].resize(3);
	}
	//particles cicle
	//double dx, dy, dz, tx, ty, tz;
	std::cout << "Interpolating forces into particles ";
	size_t x, y, z = 0;
	for (size_t i = 0; i < Particles.size(); ++i)
	{
		//calculating indexes of parent cell
		double dx = Particles[i][0] - std::floor((Particles[i][0]) / H);
		double dy = Particles[i][1] - std::floor((Particles[i][1]) / H);
		double dz = Particles[i][2] - std::floor((Particles[i][2]) / H);
		double tx = 1. - dx;
		double ty = 1. - dy;
		double tz = 1. - dz;
		if (std::floor(Particles[i][0] / H) < dim &&
			std::floor(Particles[i][1] / H) < dim &&
			std::floor(Particles[i][2] / H) < dim) {
			x = std::floor(Particles[i][0] / H);
			y = std::floor(Particles[i][1] / H);
			z = std::floor(Particles[i][2] / H);
		}
		std::cout << ".";
		a[i][0] = g[x][y][z][0] * tx*ty*tz + g[x + 1][y][z][0] * dx*ty*tz +
			g[x][y + 1][z][0] * tx*dy*tz + g[x + 1][y + 1][z][0] * dx*dy*tz +
			g[x][y][z + 1][0] * tx*ty*dz + g[x + 1][y][z + 1][0] * dx*ty*dz +
			g[x][y + 1][z + 1][0] * tx*dy*dz + g[x + 1][y + 1][z + 1][0] * dx*dy*dz;
		std::cout << ".";
		a[i][1] = g[x][y][z][1] * tx*ty*tz + g[x + 1][y][z][1] * dx*ty*tz +
			g[x][y + 1][z][1] * tx*dy*tz + g[x + 1][y + 1][z][1] * dx*dy*tz +
			g[x][y][z + 1][1] * tx*ty*dz + g[x + 1][y][z + 1][1] * dx*ty*dz +
			g[x][y + 1][z + 1][1] * tx*dy*dz + g[x + 1][y + 1][z + 1][1] * dx*dy*dz;
		std::cout << ".";
		a[i][2] = g[x][y][z][2] * tx*ty*tz + g[x + 1][y][z][2] * dx*ty*tz +
			g[x][y + 1][z][2] * tx*dy*tz + g[x + 1][y + 1][z][2] * dx*dy*tz +
			g[x][y][z + 1][2] * tx*ty*dz + g[x + 1][y][z + 1][2] * dx*ty*dz +
			g[x][y + 1][z + 1][2] * tx*dy*dz + g[x + 1][y + 1][z + 1][2] * dx*dy*dz;

	}
	g.clear();
	std::cout << " DONE." << std::endl;
	
}