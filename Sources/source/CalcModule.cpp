#include "CaclModule.h"
#include <vector>
//#include <algorithm>
#include "Fourier_tools.hh"
#include "particle_struct.h"
#include "box_struct.h"

void ScalePos(std::vector<std::vector<double> >& Particles, double scale)
{
    std::size_t size = Particles.size();
    double max = ai::max(Particles);
	// double min = ai::min(Particles);
	std::cout << "diff = "<<max <<std::endl;
	for (size_t i = 0; i < size; i++) {
        Particles[i][0]/=scale;
        Particles[i][1]/=scale;
        Particles[i][2]/=scale;
    }
}


void SetSun(std::vector<
 std::vector <
   std::vector<
     std::vector<double> > > >& density,
   double mass,
 size_t dim)
 {
   size_t center = 0.5 * dim;
   density[center][center][center][0] = 1000. * mass;
 }

void GenRing(std::vector<std::vector<double> >& Particles, std::vector<std::vector<double> >& vel, double number_particles, double L)
{
    srand (time (NULL) );
    double R1 = L / 8.;
    double R2 =  L / 4.;
    double Di = R2-R1;
    double center =0.5*L;
    double x;

    for (size_t i = 0; i < number_particles ; ++i)
    {

            x = R1 + rand() % (int)Di;
            // x = R1+2.;
                //std::cout<<"x = "<<x<<std::endl;
                    //
                    Particles[i][0] =x * cos(i) +center;
                    Particles[i][1] = x * sin(i)+center;
                    Particles[i][2] = center;

                    vel[i][0] =  1.0*sqrt(1000./x) * cos(0.5*3.14159265359 - (double)i);
                    vel[i][1] = -1.0*sqrt(1000./x) * sin(0.5* 3.14159265359 - (double)i);
                    //vel[i][2] = 0.;//0.3 *cos(3.14159265359 - (double)i);
    }


}
void PuttoBox(std::vector<
					std::vector<double> >& Particles ,
					std::vector<std::vector<size_t> >&box,
					double H,
				 size_t dim)
{
		H*=1; //cell size is  H
		size_t x, y, z=0;

		for (size_t i =0 ; i< Particles.size(); ++i)
		{
			if (std::floor(Particles[i][0] / H) < dim &&
			std::floor(Particles[i][1] / H)  < dim &&
			std::floor(Particles[i][2] / H)  < dim) {
			x = std::floor(Particles[i][0] / H);
			y = std::floor(Particles[i][1] / H);
			z = std::floor(Particles[i][2] / H);
			box[i][0] = x;
			box[i][1] = y;
			box[i][2] = z;
			// std::cout<<" in ";
			}
		}
		// std::cout<<"puttobox"<<std::endl;
		// ai::printMatrix(box);
	size_t it =0;
	for (size_t i = 0 ; i<Particles.size(); ++i)
	{
		for (size_t j = 1 ; j <Particles.size(); ++j )
		{
			if(i!=j)
			if (box[i][0] == box[j][0] && box[i][1] == box[j][1] && box[i][2] == box[j][2])
                  if( box[j][3] == 0 && box[i][3]==0)
                  {
                      it++;
                    box[i][3] = box[j][3] = it;
                  }
                  else
                  {
                    box[i][3] = box[j][3] = it;
                  }

        }

    }
}

 void CaclDensity(std::vector<
					std::vector<double> >& Particles ,
				 std::vector<
					std::vector <
						std::vector<
							std::vector<double> > > >& density,
				 double mass,
  			  	 double H,
				 size_t dim)
{
	//std::fill(density.begin(), density.end(), 0.);
	//for ()
	size_t x, y, z = 0;
	double dx, dy, dz, tx, ty, tz;
    double h = 1./(double)(dim-1);
    mass/=(h*h*h);
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
			//std::cout<<"x = "<<x <<"   y = "<<y << "  z = "<<z<<std::endl;

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






void GetAccelPM(std::vector<std::vector<double>>& Particles,
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
				//std::cout<<"x = "<<x <<"  y = "<<y <<"  z = "<<z <<std::endl;
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


void GetAccel(std::vector<std::vector<double>>& Particles,
           std::vector<std::vector<std::vector<std::vector<double > > > >& density,
                      std::vector<std::vector<size_t> > &box,
            std::vector<std::vector<double>>& a, double H)
{
  size_t dim = density.size();
  // double h = 1./H;
std::vector<std::vector<double > > as;
//as-матрица с посчитанными значениями ускорений
as.resize(Particles.size());
for (size_t i =0 ; i< Particles.size() ; ++i)
{
    as[i].resize(3);
}
  as=a;
  double mass = 1.0;
  for (size_t i = 0; i < Particles.size(); ++i)
  {
        //calculating indexes of parent cell

      //std::cout<<"x = "<<x <<"  y = "<<y <<"  z = "<<z <<std::endl;

      // find number of boxes with multiple particles


            std::vector <std::vector <double > > par;
            std::vector < std::vector <double > > ac;
            //Записываем частицы, которые попадают в радиус в 2 ячейки
            par.push_back(std::vector<double>{Particles[i][0] , Particles[i][1] , Particles[i][2]});
            ac.push_back(std::vector<double>{a[i][0] , a[i][1] , a[i][2]});
            // ac.push_back(std::vector<double>{0. , 0. , 0.});
            for (size_t j = 0 ; j<box.size(); ++j)
            {
                if (i!=j){
                if (box[i][3]!= 0 && box[i][3] == box[j][3]) // Записываем частицы в в i-ой ячейке
                {
                  par.push_back(std::vector<double>{ Particles[j][0] , Particles[j][1] , Particles[j][2]});
                  ac.push_back(std::vector<double>{ a[j][0] , a[j][1] , a[j][2]});
                  // ac.push_back(std::vector<double>{0. , 0. , 0.});

                }
                if (box[i][0]!=box[j][0] &&
                  box[i][1]!=box[j][1] &&
                  box[i][2]!=box[j][2] &&
                  abs((int)box[i][0] - (int)box[j][0]) <= 3 &&
                  abs((int)box[i][1] - (int)box[j][1]) <= 3 &&
                  abs((int)box[i][2] - (int)box[j][2]) <= 3 )
                  {
                    //std::cout<<"in";
                    par.push_back(std::vector<double>{ Particles[j][0] , Particles[j][1] , Particles[j][2]} );
                    ac.push_back(std::vector<double>{ a[j][0] , a[j][1] , a[j][2]});
                    // ac.push_back(std::vector<double>{0. , 0. , 0.});

                  }

                }

              }
            //  std::cout <<"Nearest particles coords"<<std::endl;
            // ai::printMatrix(par);

            // std::cout<<"their accel"<<std::endl;
            // ai::printMatrix(ac);
            as[i] = DirectPM(par, ac , mass);

            // as.push_back(std::vector<double>{ac[i][0], ac[i][1], ac[i][2]});
          //  as[i] = ac[i];
            // std::cout<<"accel step"<<std::endl;
            // ai::printMatrix(as);
  }
      // std::cout<<"a matrix"<<std::endl;
      // ai::printMatrix(a);
      a = as;
      // for (size_t i = 0 ; i<a.size(); ++i)
      // {
      //   a[i][0]-=as[i][0];
      //   a[i][1]-=as[i][1];
      //   a[i][2]-=as[i][1];
      // }
}

 std::vector<double> DirectPM(std::vector<std::vector<double> > & Particles, std::vector<std::vector<double> > & dir, double mass)
	{
	 double rsr = 3.0;
	 double ksix ;
	 double ksiy ;
	 double ksiz ;
	 double dx;
	 double dy;
	 double dz;
	 double distij;
	 double magi;
	 double nagi;
	 double Rx , Ry, Rz;
	 double gx, gy, gz;
	 double a = 0.3; //равновесное расстояние
	 // std::cout<<"input accel"<<std::endl;
	 // ai::printMatrix(dir);
	  for(size_t i = 0 ; i<1; ++i)
	  {
		for(size_t j = i+1 ; j<Particles.size(); ++j)
		{

		   dx = Particles[i][0] - Particles[j][0];
		   dy = Particles[i][1] - Particles[j][1];
		   dz = Particles[i][2] - Particles[j][2];

		   // std::cout<<"dx = "<<dx<<std::endl;
		   distij = sqrt(dx*dx + dy*dy +dz*dz);
		   magi = (a*a) /(distij*distij*distij);
		   nagi = (1. *std::pow(a,6))/(std::pow(distij, 7));
		   std::cout<<"distij = "<<distij<<std::endl;

				 ksix = 2.*std::abs(dx)/rsr;
				 ksiy = 2.*std::abs(dy)/rsr;
				 ksiz = 2.*std::abs(dz)/rsr;
				 // std::cout<<"ksix = "<<ksix<<std::endl;
				if ( distij >0 && distij <= a)
				{
					Rx = (1./(35.* rsr*rsr))*(224.* ksix -
					224.*ksix*ksix*ksix +
					70.*ksix*ksix*ksix*ksix +
					48.*ksix*ksix*ksix*ksix*ksix -
					21.*ksix*ksix*ksix*ksix*ksix*ksix);
					Ry = (1./(35.* rsr*rsr))*(224.* ksiy -
					224.*ksiy*ksiy*ksiy +
					70.*ksiy*ksiy*ksiy*ksiy +
					48.*ksiy*ksiy*ksiy*ksiy*ksiy -
					21.*ksiy*ksiy*ksiy*ksiy*ksiy*ksiy);
					Rz = (1./(35.* rsr*rsr))*(224.* ksiz -
					224.*ksiz*ksiz*ksiz +
					70.*ksiz*ksiz*ksiz*ksiz +
					48.*ksiz*ksiz*ksiz*ksiz*ksiz -
					21.*ksiz*ksiz*ksiz*ksiz*ksiz*ksiz);

					// std::cout<<"PM accel dir["<<i<<"][0] = "<<dir[i][0]<<std::endl;
					// std::cout<<"nagi*dx = "<<nagi*dx<<"  magi*dx = "<<magi*dx<<std::endl;
					std::cout<<"Rx = "<<Rx<<std::endl;
					//std::cout<<"gx = "<<gx<<std::endl;
					dir[i][0] += (nagi*dx - dx*magi)/(a * a) + Signum(dx)*Rx;//direct force
					dir[i][1] += (nagi*dy - magi*dy)/(a * a) +  Signum(dy)*Ry;
					dir[i][2] += (nagi*dz - magi*dz)/(a * a) +  Signum(dz)*Rz;

					std::cout<<"dir["<<i<<"][0] = "<<dir[i][0]<<std::endl;


				}

				if ( distij > a && distij <= rsr/2.)
				{

					Rx = (1./(35.* rsr*rsr))*(224.* ksix -
					224.*ksix*ksix*ksix +
					70.*ksix*ksix*ksix*ksix +
					48.*ksix*ksix*ksix*ksix*ksix -
					21.*ksix*ksix*ksix*ksix*ksix*ksix);
					Ry = (1./(35.* rsr*rsr))*(224.* ksiy -
					224.*ksiy*ksiy*ksiy +
					70.*ksiy*ksiy*ksiy*ksiy +
					48.*ksiy*ksiy*ksiy*ksiy*ksiy -
					21.*ksiy*ksiy*ksiy*ksiy*ksiy*ksiy);
					Rz = (1./(35.* rsr*rsr))*(224.* ksiz -
					224.*ksiz*ksiz*ksiz +
					70.*ksiz*ksiz*ksiz*ksiz +
					48.*ksiz*ksiz*ksiz*ksiz*ksiz -
					21.*ksiz*ksiz*ksiz*ksiz*ksiz*ksiz);
					// std::cout<<"Rx = "<<Rx<<std::endl;

					//std::cout<<"gx = "<<gx<<std::endl;
					dir[i][0] += dx*(nagi - magi)/(a * a) + Signum(dx)*Rx;//direct force
					dir[i][1] += dy*(nagi - magi)/(a * a) + Signum(dy)*Ry;
					dir[i][2] += dz*(nagi - magi)/(a * a) + Signum(dz)*Rz;
					std::cout<<"Rx = "<<Rx<<std::endl;
					std::cout<<"dir[i][0] = "<<dir[i][0]<<std::endl;


				}
				if(distij > 0.5*rsr && distij <= rsr)
				{
					Rx = (1./(35.* rsr*rsr))*(12./(ksix*ksix) - 224. + 896.*ksix -
					840.*ksix*ksix +
					224. *ksix*ksix*ksix +
					70.*ksix*ksix*ksix*ksix -
					48.*ksix*ksix*ksix*ksix*ksix +
					7.*ksix*ksix*ksix*ksix*ksix*ksix);
					Ry = (1./(35.* rsr*rsr))*(12./(ksiy*ksiy) - 224. + 896.*ksiy -
					840.*ksiy*ksiy +
					224. *ksiy*ksiy*ksiy +
					70.*ksiy*ksiy*ksiy*ksiy -
					48.*ksiy*ksiy*ksiy*ksiy*ksiy +
					7.*ksiy*ksiy*ksiy*ksiy*ksiy*ksiy);

					Rz = (1./(35.* rsr*rsr))*(12./(ksiz*ksiz) - 224. + 896.*ksiz -
					840.*ksiz*ksiz +
					224. *ksiz*ksiz*ksiz +
					70.*ksiz*ksiz*ksiz*ksiz -
					48.*ksiz*ksiz*ksiz*ksiz*ksiz +
					7.*ksiz*ksiz*ksiz*ksiz*ksiz*ksiz);

					//  std::cout<<"gx = "<<gx<<std::endl;
					dir[i][0] += dx*(nagi - magi)/(a * a) + Signum(dx)*Rx; //direct force
					dir[i][1] += dy*(nagi - magi)/(a * a) + Signum(dy)*Ry;
					dir[i][2] += dz*(nagi - magi)/(a * a) + Signum(dz)*Rz;
					std::cout<<"Rx = "<<Rx<<std::endl;
					std::cout<<"dir[i][0] = "<<dir[i][0]<<std::endl;

				}
		}

	  }
	  // std::cout<<"dir matrix"<<std::endl;
	  // ai::printMatrix(dir);
	  return std::vector<double>{dir[0][0], dir[0][1], dir[0][2]};
	}


double Signum(double  x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}


void Direct(std::vector<std::vector<double>> & Particles, std::vector<std::vector<double>> & dir, double mass)
{
			double dx ;
			double dy ;
			double dz ;
			double distij;
			double magj, magi;
			double a = 0.3;//равновесное расстояние
			double close; // слагаемое, отвечающее за близкодействие
            //double s = 1500;

	  for(size_t i =0 ; i<Particles.size(); ++i)
		for(size_t j = i+1 ; j<Particles.size(); ++j)
		{
			dx = Particles[i][0] - Particles[j][0];
			dy = Particles[i][1] - Particles[j][1];
			dz = Particles[i][2] - Particles[j][2];
			distij = sqrt(dx*dx + dy*dy +dz*dz);
			magi = (a*a) /(distij*distij*distij);
			close = (1. *std::pow(a,6))/(std::pow(distij, 7));
            // close =std::exp(-2.0*a*s*(distij - a))/distij;;
			// close =0.;
			dir[i][0] -= -mass*mass*(close*dx - magi*dx)/(a*a);
			dir[i][1] -= -mass*mass*(close*dy - magi*dy)/(a*a);
			dir[i][2] -= -mass*mass*(close*dz - magi*dz)/(a*a);
            // dir[i][0] -=  (std::pow(a,13)/std::pow(dx,13) - Signum(dx)*(a*a)/(dx*dx))/(a*a);
            // dir[i][1] -=  (std::pow(a,13)/std::pow(dy,13) - Signum(dy)*(a*a)/(dy*dy))/(a*a);
            // dir[i][2] -=  (std::pow(a,13)/std::pow(dz,13) - Signum(dz)*(a*a)/(dz*dz))/(a*a);
			magj = (1.0*mass*a*a)/(distij*distij*distij);
			dir[j][0] += -mass*mass*(close*dx - magj*dx)/(a*a);
			dir[j][1] += -mass*mass*(close*dy - magj*dy)/(a*a);
			dir[j][2] += -mass*mass*(close* dz - magj*dz)/(a*a);
            // dir[j][0] +=  (std::pow(a,13)/std::pow(dx,13) - Signum(dx)*(a*a)/(dx*dx))/(a*a);
            // dir[j][1] +=  (std::pow(a,13)/std::pow(dy,13) - Signum(dy)*(a*a)/(dy*dy))/(a*a);
            // dir[j][2] +=  (std::pow(a,13)/std::pow(dz,13) - Signum(dz)*(a*a)/(dz*dz))/(a*a);

            std::cout<<"DIST ij = "<<distij<<std::endl;
            std::cout<<"dirij x= "<<dir[i][0] <<"    dirji x=  "<<dir[j][0]<<std::endl;
			std::cout<<"dirij y= "<<dir[i][1] <<"    dirji y=  "<<dir[j][1]<<std::endl;
			std::cout<<"dirij z= "<<dir[i][2] <<"    dirji z=  "<<dir[j][2]<<std::endl;

		}
}
