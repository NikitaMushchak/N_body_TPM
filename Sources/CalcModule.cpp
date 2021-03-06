#include "CaclModule.h"
#include <vector>
//#include <algorithm>
#include "Fourier_tools.hh"
#include "particle_struct.h"
#include "box_struct.h"

void ScalePos(std::vector<std::vector<double> >& Particles, double scale)
{
    std::size_t size = Particles.size();
    for (size_t i = 0; i < size; i++) {
        Particles[i][0]/=scale;
        Particles[i][1]/=scale;
        Particles[i][2]/=scale;
    }
}

 void CaclDensity(std::vector<
					std::vector<double>>& Particles ,
				 std::vector<
					std::vector <
						std::vector<
							std::vector<double> > > >& density,
				 float mass,
  			  	 double H,
				 size_t dim)
{
	//std::fill(density.begin(), density.end(), 0.);
	//for ()
	size_t x, y, z = 0;
	double dx, dy, dz, tx, ty, tz;
    double h = 1./(double)(dim-1);
    mass/=(float)(h*h*h);
    // std::cout <<"mass scaled"<<mass<<std::endl;

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
	double h = 1. / double(dim - 1);

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
			fft2(f);
			for (size_t k = 0; k < dim; ++k) {
				rho[i][j][k][0] = f[k][0];
				rho[i][j][k][1] = f[k][1];

			}
		}
	}

	for (size_t i = 0; i < dim; ++i)
	{
		std::vector <
			std::vector<double>> f; //for rows and colomns
		f.resize(dim);
		for (size_t l = 0; l < dim; l++)
		{
			f[l].resize(2);
		}
		//FFT coloumns of density
		for (size_t k = 0; k < dim; ++k) {
			for (size_t j = 0; j < dim; ++j) {
				f[j][0] = rho[i][j][k][0];
				f[j][1] = rho[i][j][k][1];

			}
			fft2(f);
			for (size_t j = 0; j < dim; ++j) {
				rho[i][j][k][0] = f[j][0];
				rho[i][j][k][1] = f[j][1];
			}
		}

	}

	//FFT in 3rd dimention
	for (size_t k = 0; k < dim; ++k) {
		std::vector <
			std::vector<double>> f;
		f.resize(dim);
		for (size_t l = 0; l < dim; l++)
		{
			f[l].resize(2);
		}

		for (size_t j = 0; j < dim; j++) {
			for (size_t i = 0; i< dim; i++)
			{
				f[i][0] = rho[i][j][k][0];
				f[i][1] = rho[i][j][k][1];
				//f[j][1] = rho[j][k][1];
			}
			fft2(f);
			for (size_t i = 0; i < dim; i++) {
				rho[i][j][k][0] = f[i][0];
				rho[i][j][k][1] = f[i][1];
				//rho[j][k][1] = f[j][1];
			}
		}
	}
	/*std::cout << "Density in Fourier space" << std::endl;
	for (size_t i = 0; i < dim; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			for (size_t k = 0; k < dim; k++)
			{
				std::cout << "(" << i << " , " << j << " , " << k << " )" << "dens = " << "(" << rho[i][j][k][0] << " , " << rho[i][j][k][1] << ")" << std::endl;
			}
		}

	}*/
		//std::cout << "FFT for layers ....DONE" << std::endl;

		// solve equation in Fourier space

		double pi = 3.14159265358979;

		//std::vector<double> W = { cos(2.0 * pi / double(dim)) , sin(2.0 * pi / double(dim)) };

		double W0 = cos(2.0 * pi / double(dim));
		double W1 = sin(2.0 * pi / double(dim));
		//std::vector<double> Wm = { 1.0, 0. }, Wn = { 1.0, 0. }, Wk = { 1.0, 0. };
		double Wm0 = 1.;
		double Wm1 = 0.;
		double Wn0 = 1.;
		double Wn1 = 0.;
		double Wk0 = 1.;
		double Wk1 = 0.;

		for (size_t k = 0; k < dim; k++) {
			for (size_t m = 0; m < dim; m++) {
				for (size_t n = 0; n < dim; n++) {

				double denom0 = 6.;
				double denom1 = 0.;
				//std::vector < double> i_Wk = { Wk0 / (Wk0 * Wk0 + Wk1 * Wk1), -(Wk1 / (Wk0 * Wk0 + Wk1 * Wk1)) };
				double i_Wk0 = Wk0 / (Wk0 * Wk0 + Wk1 * Wk1);
				double i_Wk1 = -(Wk1 / (Wk0 * Wk0 + Wk1 * Wk1));
				//std::vector < double> i_Wm = { Wm0 / (Wm0 * Wm0 + Wm1 * Wm1), -(Wm1 / (Wm0 * Wm0 + Wm1 * Wm1)) };
				double i_Wm0 = Wm0 / (Wm0 * Wm0 + Wm1 * Wm1);
				double i_Wm1 = -(Wm1 / (Wm0 * Wm0 + Wm1 * Wm1));
				//std::vector < double> i_Wn = { Wn0 / (Wn0 * Wn0 + Wn1 * Wn1), -(Wn1 / (Wn0 * Wn0 + Wn1 * Wn1)) };
				double i_Wn0 = Wn0 / (Wn0 * Wn0 + Wn1 * Wn1);
				double i_Wn1 = -(Wn1 / (Wn0 * Wn0 + Wn1 * Wn1));


				denom0 -= (Wk0 + i_Wk0 + Wm0 + i_Wm0 + Wn0 + i_Wn0);
				denom1 -= (Wk1 + i_Wk1 + Wm1 + i_Wm1 + Wn1 + i_Wn1);
				if (denom0 !=0. )
				{

					double q =  -rho[k][m][n][0];
          double b =  -rho[k][m][n][1];
					rho[k][m][n][0] =  4.*pi*h*h*h*(q * denom0 + b * denom1) / (denom0 * denom0 + denom1 * denom1);
					rho[k][m][n][1] =  4.*pi*h*h*h*(b * denom0 - q * denom1) / (denom0 * denom0 + denom1 * denom1);

				}
				double b = Wn0;
				Wn0 = Wn0 * W0 - Wn1 * W1;
				Wn1 = Wn1 * W0 + b * W1;
			}
			double r = Wm0;
			Wm0 = Wm0 * W0 - Wm1 * W1;
			Wm1 = Wm1 * W0 + r * W1;
		}

		double a = Wk0;
		Wk0 = Wk0 * W0 - Wk1 * W1;
		Wk1 = Wk1 * W0 + a * W1;
	}

		/*std::cout << "Potential in Fourier space" << std::endl;
		for (size_t i = 0; i < dim; i++)
		{
			for (size_t j = 0; j < dim; j++)
			{
				for (size_t k = 0; k < dim; k++)
				{
					std::cout << "(" << i << " , " << j << " , " << k << " )" << "dens = " << "(" << rho[i][j][k][0] << " , " << rho[i][j][k][1] << ")" << std::endl;
				}
			}
		}*/
	//std::cout << "Solve Poisson eq. in Fourier space" << std::endl;
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
				ifft2(f);
				for (size_t k = 0; k < dim; k++)
				{
					rho[i][j][k][0] = f[k][0];
					rho[i][j][k][1] = f[k][1];
				}
			}
		}
		for (size_t i = 0; i < dim; i++)
		{
			std::vector<std::vector<double>> f;
			f.resize(dim);
			for (size_t y = 0; y < dim; y++)
			{
				f[y].resize(2);
			}
		// inverse FFT columns of rho
		for (size_t k = 0; k < dim; k++) {
			for (size_t j = 0; j < dim; j++)
			{
				f[j][0] = rho[i][j][k][0];
				f[j][1] = rho[i][j][k][1];
			}
			ifft2(f);
			for (size_t j = 0; j < dim; j++)
			{
				rho[i][j][k][0] = f[j][0];
				rho[i][j][k][1] = f[j][1];
			}
		}
	}

		//IFFT in 3rd dimention

		for (size_t k = 0; k < dim; ++k) {
			std::vector<std::vector<double>> f;
			f.resize(dim);
			for (size_t y = 0; y < dim; y++)
			{
				f[y].resize(2);
			}

			for (size_t j = 0; j < dim; j++) {
				for (size_t i = 0; i< dim; i++)
				{
					f[i][0] = rho[i][j][k][0];
					f[i][1] = rho[i][j][k][1];
					//f[j][1] = rho[j][k][1];
				}
				ifft2(f);
				for (size_t i = 0; i < dim; i++) {
					rho[i][j][k][0] = f[i][0];
					rho[i][j][k][1] = f[i][1];
					//rho[j][k][1] = f[j][1];
				}
			}
		}
	//std::cout << "Solving ... DONE" << std::endl;
}

  void GetAccel(std::vector<std::vector<double>>& Particles,
						std::vector<std::vector<std::vector<std::vector<double >>>>& density,
							std::vector<std::vector<double>>& a, double H)
{

	size_t dim = density.size();

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

            // x componet
            	a[i][0] = -0.5*(density[x + 1][y][z][0] - density[x - 1][y][z][0])*tx*ty*tz- 0.5*(density[x + 2][y][z][0] - density[x][y][z][0])*dx*ty*tz
            		- 0.5*(density[x + 1][y + 1][z][0] - density[x - 1][y + 1][z][0])*tx*dy*tz - 0.5*(density[x + 2][y + 1][z][0] - density[x][y + 1][z][0])*dx*dy*tz
            		- 0.5*(density[x + 1][y][z + 1][0] - density[x - 1][y][z + 1][0])*tx*ty*dz - 0.5 *(density[x + 2][y][z + 1][0] - density[x][y][z + 1][0])*dx*ty*dz
            		- 0.5*(density[x + 1][y + 1][z + 1][0] - density[x - 1][y + 1][z + 1][0])*tx*dy*dz - 0.5*(density[x + 2][y + 1][z + 1][0] - density[x][y + 1][z + 1][0])*dx*dy*dz;
            	// std::cout << "g[i].x = " << g[i][0] << std::endl;
            	//y component
            	a[i][1] = -0.5*(density[x][y + 1][z][0] - density[x][y - 1][z][0])*tx*ty*tz - 0.5*(density[x + 1][y+1][z][0] - density[x+1][y-1][z][0])*dx*ty*tz
            		- 0.5*(density[x][y + 2][z][0] - density[x][y][z][0])*tx*dy*tz - 0.5*(density[x + 1][y + 2][z][0] - density[x+1][y][z][0])*dx*dy*tz
            		- 0.5*(density[x][y+1][z + 1][0] - density[x][y-1][z + 1][0])*tx*ty*dz - 0.5 *(density[x + 1][y+1][z + 1][0] - density[x+1][y-1][z + 1][0])*dx*ty*dz
            		- 0.5*(density[x][y + 2][z + 1][0] - density[x][y][z + 1][0])*tx*dy*dz - 0.5*(density[x + 1][y + 2][z + 1][0] - density[x+1][y][z + 1][0])*dx*dy*dz;
            	// std::cout << "g[i].y = " << g[i][1] << std::endl;
            	//z component
            	a[i][2] = -0.5*(density[x][y][z + 1][0] - density[x][y][z - 1][0])*tx*ty*tz -0.5*(density[x + 1][y][z + 1][0] - density[x + 1][y][z - 1][0])*dx*ty*tz
            		- 0.5*(density[x][y + 1][z+1][0] - density[x][y+1][z-1][0])*tx*dy*tz - 0.5*(density[x + 1][y + 1][z+1][0] - density[x + 1][y+1][z-1][0])*dx*dy*tz
            		- 0.5*(density[x][y][z + 2][0] - density[x][y][z][0])*tx*ty*dz - 0.5 *(density[x + 1][y][z + 2][0] - density[x + 1][y][z][0])*dx*ty*dz
            		- 0.5*(density[x][y + 1][z + 2][0] - density[x][y+1][z][0])*tx*dy*dz - 0.5*(density[x + 1][y + 1][z + 2][0] - density[x + 1][y+1][z][0])*dx*dy*dz;
            	// std::cout << "g[i].z = " << g[i][2] << std::endl;
            //}
	}

	//std::cout << " DONE." << std::endl;

}

double Signum(double  x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}
void Direct(std::vector<std::vector<double>> & Particles, std::vector<std::vector<double>> & dir, double mass)
{
  // for (size_t i =0 ; i<Particles.size(); ++i )
  // {
  //   for (size_t j =0; j< Particles.size(); ++j)
  //   {
  //     if (i != j)
  //     {
  //       dir[i][0]+= 1. /((Particles[i][0]-Particles[j][0])*(Particles[i][0]-Particles[j][0])*(Particles[i][0]-Particles[j][0]));
  //       dir[i][1]+= 1. /((Particles[i][1]-Particles[j][1])*(Particles[i][1]-Particles[j][1])*(Particles[i][1]-Particles[j][1]));
  //       dir[i][2]+= 1. /((Particles[i][2]-Particles[j][2])*(Particles[i][2]-Particles[j][2])*(Particles[i][2]-Particles[j][2]));
  //       // std::cout<<"dir"<< " i = "<<i<<" j = "<<j <<"   dir = "<<dir[i][0]<<std::endl;
  //     }
  //   }
  // }

  for(size_t i =0 ; i<Particles.size(); ++i)
  {
    for(size_t j = i+1 ; j<Particles.size(); ++j)
    {
        double dx = Particles[i][0] - Particles[j][0];
        double dy = Particles[i][1] - Particles[j][1];
        double dz = Particles[i][2] - Particles[j][2];
        double distij = sqrt(dx*dx + dy*dy +dz*dz);
        double magi = (1.0*mass) /(distij*distij*distij);
        dir[i][0] -= magi*dx;
        dir[i][1] -= magi*dy;
        dir[i][2] -= magi*dz;
        double magj = (1.0*mass)/(distij*distij*distij);
        dir[j][0] += magi*dx;
        dir[j][1] += magi*dy;
        dir[j][2] += magi*dz;
    }

  }
}
