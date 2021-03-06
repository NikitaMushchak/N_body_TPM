#include "CaclModule.h"
#include <vector>
#include <random>
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
   density[center][center][center][0] = 195000000. * mass;
 }

void GenRing(std::vector<std::vector<double> >& Particles, std::vector<std::vector<double> >& vel, double number_particles, double L)
{
    std::random_device random_device;
    std::mt19937 generator(random_device());
    double R1 = L / 16.;
    double R2 =  L / 6.;
    double Di = R2-R1;
    double z, rh ;
    double h = 0.1* R2;
    double center =0.5*L;
    double x;

    for (size_t i = 0; i < number_particles ; ++i)
    {
            std::uniform_real_distribution<> dis(R1,R2);
            x = dis(generator);
			
            rh = h * std::sqrt(1 - x*x/(R2*R2));
            std::uniform_real_distribution<> dis1(-rh,rh);
            //std::cout<<"rh "<<rh<<std::endl;
            z = dis1(generator);
            //x = R1+2.;
                //std::cout<<"x = "<<x<<std::endl;
                    //
                    Particles[i][0] = x * cos(i) +center;
                    Particles[i][1] = x * sin(i) +center;
                    Particles[i][2] = z + center;

                    vel[i][0] =  0.17*sqrt((10.*number_particles)/std::sqrt(x*x + z*z)) * cos(0.5*3.14159265359 - (double)i);
                    vel[i][1] = -0.17*sqrt((10.*number_particles)/std::sqrt(x*x + z*z)) * sin(0.5* 3.14159265359 - (double)i);


                    //std::cout<<"vel[i][0] "<<vel[i][0]<<  "  "<<"vel[i][1] "<<vel[i][1]<<std::endl;
                    //vel[i][2] = 0.;//0.3 *cos(3.14159265359 - (double)i);
    }


}


void GenEllipsoid(std::vector<std::vector<double> >& Particles,
					std::vector<std::vector<double> >& vel, double number_particles, double L)
{
    srand (time (NULL) );
    double R1 = L / 4.;
    double R2 =  L / 2.;
	double pi = 3.14159265359;

	double ws ;

    double Di = R2-R1;
    double center =0.5*L;
    double x; // генерированный радиус

    for (size_t i = 0; i < number_particles ; ++i)
    {

            x = R1 + rand() % (int)Di;

			ws =
			std::sqrt((3.*pi*number_particles)/(4.* std::pow(x , 3))) +
				 std::sqrt( (10.*number_particles)/(std::pow(x , 3)));
				//ws= 0.6*ws;
            //x = R1+2.;
                 //std::cout<<"x = "<<x<<std::endl;
                // std::cout<<"ws = "<<ws<<std::endl;
                    //
                    Particles[i][0] =x * cos(i) +center;
                    Particles[i][1] = x * sin(i)+center;
                    Particles[i][2] = center;

                    // vel[i][0] =  ws *x* cos(0.5* pi - (double)i);
                    // vel[i][1] = -ws *x* sin(0.5* pi - (double)i);
                    //vel[i][2] = 0.;//0.3 *cos(3.14159265359 - (double)i);

					vel[i][0] = - ws*x*sin(i);
					vel[i][1] = ws*x*cos(i);
					// std::cout<<"vel x= "<<vel[i][0]<<std::endl;
					// std::cout<<"vel y= "<<vel[i][1]<<std::endl;
					//std::cout<<"ws = "<<std::sqrt(vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1])/x<<std::endl;
    }


}


void Geninit(std::vector<std::vector<double> >& Particles,
					std::vector<std::vector<double> >& vel, double number_particles, double L)
{
    std::random_device random_device;
    std::mt19937 generator(random_device());
    double R1 = L / 4.;
    double R2 =  L / 2.;
	double pi = 3.14159265359;
    double h = 0.1 * R2;
	double ws;
    double rh;
    double Di = R2-R1;
    double center = 0.5*L;
    double x; // генерированное расстояние по плоскости
    double z; // генергированаая рсстояние по глубине


    for (size_t i = 0; i < number_particles ; ++i)
    {

            std::uniform_real_distribution<> dis(R1,R2);
            x = dis(generator);
            rh = h * std::sqrt(1 - x*x/(R2*R2));
            std::uniform_real_distribution<> dis1(-rh,rh);
            //std::cout<<"rh "<<rh<<std::endl;
            z = dis1(generator);
            //std::cout<<"z = "<<z<<std::endl;
			ws =
			//std::sqrt((3.*pi*number_particles)/(4.* std::pow(x , 3))) +
				 std::sqrt( (10.*number_particles) / (std::pow(std::sqrt(x*x + z*z) , 3) ) );
				ws= 0.6*ws;
            //x = R1+2.;
                 //std::cout<<"x = "<<x<<std::endl;
                // std::cout<<"ws = "<<ws<<std::endl;
                    //
                    Particles[i][0] = x * cos(i) +center;
                    Particles[i][1] = x * sin(i)+center;
                    Particles[i][2] = z + center;

                    // vel[i][0] =  ws *x* cos(0.5* pi - (double)i);
                    // vel[i][1] = -ws *x* sin(0.5* pi - (double)i);
                    //vel[i][2] = 0.;//0.3 *cos(3.14159265359 - (double)i);

					vel[i][0] = -ws*x*sin(i);
					vel[i][1] = ws*x*cos(i);
					// std::cout<<"vel x= "<<vel[i][0]<<std::endl;
					// std::cout<<"vel y= "<<vel[i][1]<<std::endl;
					//std::cout<<"ws = "<<std::sqrt(vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1])/x<<std::endl;
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

void CaclDensitySun(
                   std::vector<double>  Sun ,
                std::vector<
                   std::vector <
                       std::vector<
                           std::vector<double> > > >& density,
                size_t number_particles,
                double mass,
                double H,
                size_t dim)
{
   //std::fill(density.begin(), density.end(), 0.);
   //for ()
   size_t x, y, z = 0;
   double dx, dy, dz, tx, ty, tz;
   double h = 1./(double)(dim-1);
   mass*=10.*number_particles;
   mass/=(h*h*h);
   // std::cout <<"mass scaled"<<mass<<std::endl;

   // for (size_t i = 0; i < Particles.size(); ++i)
   // {
       //calculating indexes
       dx = Sun[0] - std::floor(Sun[0] / H);
       dy = Sun[1] - std::floor(Sun[1] / H);
       dz = Sun[2] - std::floor(Sun[2] / H);
       tx = 1. - dx;
       ty = 1. - dy;
       tz = 1. - dz;
       if (std::floor(Sun[0] / H) < dim &&
           std::floor(Sun[1] / H)  < dim &&
           std::floor(Sun[2] / H)  < dim) {
           x = std::floor(Sun[0] / H);
           y = std::floor(Sun[1] / H);
           z = std::floor(Sun[2] / H);
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

   // }
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

			density[x][y][z][0] = mass * tx *ty*tz;

			if (y + 1 >= dim) { y = y % (dim - 1); density[x][y][z][0] = mass * tx *dy*tz; }
			else { density[x][y + 1][z][0] = mass * tx *dy*tz; }
			if (z + 1 >= dim) { z = z % (dim - 1); density[x][y][z][0] = mass * tx *ty*dz; }
			else { density[x][y][z + 1][0] = mass * tx *ty*dz; }
			if (y + 1 >= dim && z + 1 >= dim) { z = z % (dim - 1); y = y % (dim - 1); density[x][y][z][0] = mass * tx *dy*dz; }
			else { density[x][y + 1][z + 1][0] = mass * tx *dy*dz; }
			if (x + 1 >= dim) { x = x % (dim - 1); density[x][y][z][0] = mass * dx *ty*tz; }
			else { density[x + 1][y][z][0] = mass * dx *ty*tz; }
			if (x + 1 >= dim && y + 1 >= dim) { x = x % (dim - 1); y = y % (dim - 1); density[x][y][z][0] = mass * dx *dy*tz; }
			else { density[x + 1][y + 1][z][0] = mass * dx *dy*tz; }
			if (x + 1 >= dim && z + 1 >= dim) { x = x % (dim - 1); z = z % (dim - 1); density[x][y][z][0] = mass * dx *ty*dz; }
			else { density[x + 1][y][z + 1][0] = mass * dx *ty*dz; }
			if (x + 1 >= dim && y + 1 >= dim && z + 1 >= dim) { x = x % (dim - 1); y = y % (dim - 1); z = z % (dim - 1); density[x][y][z][0] = mass * dx * dy * dz; }
			else { density[x + 1][y + 1][z + 1][0] = mass * dx * dy * dz; }

		}
	}
}


  void CalcPotential(std::vector<
					std::vector<
						std::vector<
							std::vector<double>>>>& rho,
							size_t dim )
{
	double h = 1. / double(dim-1);

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


		double pi = 3.14159265358979;

		std::vector<double> W = { cos(2.0 * pi / double(dim)) , sin(2.0 * pi / double(dim)) };

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
				if (denom0 != 0. && denom1 != 0.)
				{

					double q =  -rho[k][m][n][0];
					double b =  -rho[k][m][n][1];
					rho[k][m][n][0] =  h*h*(q * denom0 + b * denom1) / ((denom0 * denom0 + denom1 * denom1));
					rho[k][m][n][1] =  h*h*(b * denom0 - q * denom1) / ((denom0 * denom0 + denom1 * denom1));

                    rho[k][m][n][0] = (rho[k][m][n][0] ) / 4.;
                    rho[k][m][n][1] = (rho[k][m][n][1] ) / 4.;

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

void CalcPotentialB(std::vector<
					std::vector<
						std::vector<
							std::vector<double> > > >& rho,
							size_t dim )
{
	double pi = 3.14159265358979;
	double h = 1. / double(dim-1);
	double sinx , siny, sinz;
	double fx ;
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

		for (size_t j = 0; j < dim; j++){
			for (size_t i = 0; i< dim; i++){
				f[i][0] = rho[i][j][k][0];
				f[i][1] = rho[i][j][k][1];
				//f[j][1] = rho[j][k][1];
			}
			fft2(f);

			sinz = 4.*sin(k * pi/dim)*sin(k*pi/dim);
			siny = 4.*sin(j * pi/dim)*sin(j*pi/dim);

			for (size_t i = 0; i < dim; i++) {

				sinx = 4.*sin(i * pi/dim)*sin(i*pi/dim);
				rho[i][j][k][0] = f[i][0];
				rho[i][j][k][1] = f[i][1];
				//rho[j][k][1] = f[j][1];
				fx = std::fabs(sinx + siny + sinz) < std::pow(10 , -8) ? std::pow(10 , -8) : sinx + siny + sinz;

				rho[i][j][k][0] /= -fx/(4.*pi*h*h*h) ;//* (4188.818841*pi);
				rho[i][j][k][1] /= -fx/(4.*pi*h*h*h) ;//* (4188.818841*pi);
			}
		}
	}
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
}




void GetAccelPM(std::vector<std::vector<double> >& Particles,
						std::vector<std::vector<std::vector<std::vector<double > > > >& density,
							std::vector<std::vector<double> >& a, double H)
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


void GetAccel(std::vector<std::vector<double> >& Particles,
				std::vector<std::vector<double> >& a,
				double H,
				double r,
				size_t dim)
{

  // double h = 1./H;
  // r - радиус применения расчета близкодействия
std::vector<std::vector<double > > as;
//as-матрица с посчитанными значениями ускорений по близкодействию
as.resize(Particles.size());
for (size_t i =0 ; i< Particles.size() ; ++i)
{
    as[i].resize(3);
}

std::vector<std::vector<double > > ab;
//ab-матрица- буфер значений
ab.resize(Particles.size());
for (size_t i =0 ; i< Particles.size() ; ++i)
{
    ab[i].resize(3);
}


  double eps = std::pow(1 , -10);
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
  double Rx ;//, Ry, Rz;
  double Sh;//поправочный коэффициент
  //double gx, gy, gz;
  double acut = 0.3; //равновесное расстояние
  double mass = 1.0;
  for (size_t i = 0; i < Particles.size(); ++i)
  {
			  for(std::size_t j = 0 ; j < Particles.size(); ++j)
			  {
				  if(i!=j)
				  {
					  if(r >= std::sqrt(
					  std::pow(Particles[i][0]-Particles[j][0] , 2) +
					  std::pow(Particles[i][1]-Particles[j][1] , 2) +
					  std::pow(Particles[i][2]-Particles[j][2] , 2) )	)
					  {
						  // par.push_back(std::vector<double>{ Particles[j][0] , Particles[j][1] , Particles[j][2]} );
						  // ac.push_back(std::vector<double>{ a[j][0] , a[j][1] , a[j][2] });
                          dx = Particles[i][0] - Particles[j][0];
               		      dy = Particles[i][1] - Particles[j][1];
               		      dz = Particles[i][2] - Particles[j][2];

               		   // std::cout << "Particles[i][0] = " << Particles[i][0] << "  Particles[j][0] = " << Particles[j][0] << std::endl;
               		   // std::cout << "Particles[i][1] = " << Particles[i][1] << "  Particles[j][1] = " << Particles[j][1] << std::endl;
               		   // std::cout << "Particles[i][2] = " << Particles[i][2] << "  Particles[j][2] = " << Particles[j][2] << std::endl;
               		   // std::cout<<"dx = "<<dx<<std::endl;
               		   // std::cout<<"dy = "<<dy<<std::endl;
               		   // std::cout<<"dz = "<<dz<<std::endl;
               		   distij = sqrt(dx*dx + dy*dy +dz*dz);
               		   magi = (acut*acut) /(distij*distij*distij);
               		   nagi = (1. *std::pow(acut,6))/(std::pow(distij, 7));
                       //nagi=0.;
               		   // std::cout<<"distij = "<<distij<<std::endl;
               				 ksix = 2.*distij/rsr;
               				if ( distij <= rsr*0.5)
               				{
               					// std::cout<<" ! "<<std::endl;
               					Rx = (1./(35.* rsr*rsr))*(224.* ksix -
               					224.*ksix*ksix*ksix +
               					70.*ksix*ksix*ksix*ksix +
               					48.*ksix*ksix*ksix*ksix*ksix -
               					21.*ksix*ksix*ksix*ksix*ksix*ksix)/ ksix;
               					// std::cout<<"Rx = "<<Rx<<std::endl;
                                   // std::cout<<" Pm dir[i][0] = "<<dir[i][0]<<std::endl;
               					//std::cout<<"gx = "<<gx<<std::endl;
               				    as[i][0] += dx*(nagi - magi)/(acut * acut) + dx*Rx ;//direct force
               					as[i][1] += dy*(nagi - magi)/(acut * acut) + dy*Rx;
               					as[i][2] += dz*(nagi - magi)/(acut * acut) + dz*Rx;
               					 // std::cout<<"Rx = "<<Rx<<std::endl;
               					 // std::cout<<"dir[i][0] = "<<dir[i][0]<<std::endl;
               					 // std::cout<<"dir[i][1] = "<<dir[i][1]<<std::endl;
               					 // std::cout<<"dir[i][2] = "<<dir[i][2]<<std::endl;
               				}

               				if(distij > 0.5*rsr && distij <= rsr)
               				{
               					// std::cout<<" ! 111111"<<std::endl;

               					Rx =ksix < std::pow(10,-8) ? 0. : (1./(35.* rsr*rsr))*(12./(ksix*ksix) - 224. + 896.*ksix -
               					840.*ksix*ksix +
               					224. *ksix*ksix*ksix +
               					70.*ksix*ksix*ksix*ksix -
               					48.*ksix*ksix*ksix*ksix*ksix +
               					7.*ksix*ksix*ksix*ksix*ksix*ksix) / ksix;
                                   // std::cout<<" PM dir[i][0] = "<<dir[i][0]<<std::endl;

               					//  std::cout<<"gx = "<<gx<<std::endl;
               					as[i][0] += dx*(nagi - magi)/(acut * acut) + dx*Rx; //direct force
               					as[i][1] += dy*(nagi - magi)/(acut * acut) + dy*Rx;
               					as[i][2] += dz*(nagi - magi)/(acut * acut) + dz*Rx;
               					 // std::cout<<"Rx = "<<Rx<<std::endl;
               					 // std::cout<<"dir[i][0] = "<<dir[i][0]<<std::endl;
               					 // std::cout<<"dir[i][1] = "<<dir[i][1]<<std::endl;
               					 // std::cout<<"dir[i][2] = "<<dir[i][2]<<std::endl;

               				}
                            if(distij < 2.)
                               {
                                   Sh = (0.0122242 + 0.609459*distij - 0.432782 * distij*distij + 0.189694 * distij*distij*distij
                                     - 0.0610997 * distij*distij * distij*distij) / distij;
                                     as[i][0] +=  Sh*dx;
                                     as[i][1] +=  Sh*dy;
                                     as[i][2] +=  Sh*dz;
                               }
                            if(distij <=3. && distij >=2)
                               {
                                   Sh = (-0.137179 + 0.158876*distij - 0.0348711 *distij*distij) / distij;
                                   as[i][0] +=  Sh*dx;
                                   as[i][1] +=  Sh*dy;
                                   as[i][2] +=  Sh*dz;
                               }
					  }
				  }
            }
              ab[i][0]=as[i][0];
              ab[i][1]=as[i][1];
              ab[i][2]=as[i][2];
  }
  for(size_t i =0 ; i< Particles.size(); ++i)
  {
     a[i][0]+=ab[i][0];
     a[i][1]+=ab[i][1];
     a[i][2]+=ab[i][2];
  }
 }

 std::vector<double> DirectPM(std::vector<std::vector<double> > & Particles, std::vector<std::vector<double> > & dir, double mass)
	{
	 double eps = std::pow(1 , -10);
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
	 double Rx ;//, Ry, Rz;
     double Sh;//поправочный коэффициент
	 //double gx, gy, gz;
	 double a = 0.3; //равновесное расстояние
	 // std::cout<<"DIRECT INNNNN !!!!!!"<<std::endl;
	 // ai::printMatrix(dir);
	  for(size_t i = 0 ; i<1; ++i)
	  {
		for(size_t j = i+1 ; j<Particles.size(); ++j)
		{
			// std::cout<<"i = "<<i<<"   j = "<<j<<std::endl;
		   dx = Particles[i][0] - Particles[j][0];
		   dy = Particles[i][1] - Particles[j][1];
		   dz = Particles[i][2] - Particles[j][2];

		   // std::cout << "Particles[i][0] = " << Particles[i][0] << "  Particles[j][0] = " << Particles[j][0] << std::endl;
		   // std::cout << "Particles[i][1] = " << Particles[i][1] << "  Particles[j][1] = " << Particles[j][1] << std::endl;
		   // std::cout << "Particles[i][2] = " << Particles[i][2] << "  Particles[j][2] = " << Particles[j][2] << std::endl;
           //
           //
		   // std::cout<<"dx = "<<dx<<std::endl;
		   // std::cout<<"dy = "<<dy<<std::endl;
		   // std::cout<<"dz = "<<dz<<std::endl;
		   distij = sqrt(dx*dx + dy*dy +dz*dz);

		   magi = (a*a) /(distij*distij*distij);
		   nagi = (1. *std::pow(a,6))/(std::pow(distij, 7));
		   // std::cout<<"distij = "<<distij<<std::endl;

				 ksix = 2.*distij/rsr;
				 // ksiy = 2.*distij/rsr;
				 // ksiz = 2.*distij/rsr;

				if ( distij <= rsr*0.5)
				{

					// std::cout<<" ! "<<std::endl;
					Rx = (1./(35.* rsr*rsr))*(224.* ksix -
					224.*ksix*ksix*ksix +
					70.*ksix*ksix*ksix*ksix +
					48.*ksix*ksix*ksix*ksix*ksix -
					21.*ksix*ksix*ksix*ksix*ksix*ksix)/ ksix;

					// Ry = (1./(35.* rsr*rsr))*(224.* ksiy -
					// 224.*ksiy*ksiy*ksiy +
					// 70.*ksiy*ksiy*ksiy*ksiy +
					// 48.*ksiy*ksiy*ksiy*ksiy*ksiy -
					// 21.*ksiy*ksiy*ksiy*ksiy*ksiy*ksiy) / ksiy;
                    //
					// Rz = (1./(35.* rsr*rsr))*(224.* ksiz -
					// 224.*ksiz*ksiz*ksiz +
					// 70.*ksiz*ksiz*ksiz*ksiz +
					// 48.*ksiz*ksiz*ksiz*ksiz*ksiz -
					// 21.*ksiz*ksiz*ksiz*ksiz*ksiz*ksiz) / ksiz;

					// std::cout<<"Rx = "<<Rx<<std::endl;
                    // std::cout<<" Pm dir[i][0] = "<<dir[i][0]<<std::endl;
					//std::cout<<"gx = "<<gx<<std::endl;
					dir[i][0] += dx*(nagi - magi)/(a * a) + dx*Rx;//direct force
					dir[i][1] += dy*(nagi - magi)/(a * a) + dy*Rx;
					dir[i][2] += dz*(nagi - magi)/(a * a) + dz*Rx;
					 // std::cout<<"Rx = "<<Rx<<std::endl;
					 // std::cout<<"dir[i][0] = "<<dir[i][0]<<std::endl;
					 // std::cout<<"dir[i][1] = "<<dir[i][1]<<std::endl;
					 // std::cout<<"dir[i][2] = "<<dir[i][2]<<std::endl;


				}
				if(distij > 0.5*rsr && distij <= rsr)
				{
					// std::cout<<" ! 111111"<<std::endl;

					Rx =ksix < std::pow(10,-8) ? 0. : (1./(35.* rsr*rsr))*(12./(ksix*ksix) - 224. + 896.*ksix -
					840.*ksix*ksix +
					224. *ksix*ksix*ksix +
					70.*ksix*ksix*ksix*ksix -
					48.*ksix*ksix*ksix*ksix*ksix +
					7.*ksix*ksix*ksix*ksix*ksix*ksix) / ksix;

					// Ry = ksiy < std::pow(10,-8)? 0. : (1./(35.* rsr*rsr))*(  12./(ksiy*ksiy) - 224. + 896.*ksiy -
					// 840.*ksiy*ksiy +
					// 224. *ksiy*ksiy*ksiy +
					// 70.*ksiy*ksiy*ksiy*ksiy -
					// 48.*ksiy*ksiy*ksiy*ksiy*ksiy +
					// 7.*ksiy*ksiy*ksiy*ksiy*ksiy*ksiy) / ksiy;
                    //
					// Rz = ksiz < std::pow(10,-8) ? 0. : (1./(35.* rsr*rsr))*(12./(ksiz*ksiz) - 224. + 896.*ksiz -
					// 840.*ksiz*ksiz +
					// 224. *ksiz*ksiz*ksiz +
					// 70.*ksiz*ksiz*ksiz*ksiz -
					// 48.*ksiz*ksiz*ksiz*ksiz*ksiz +
					// 7.*ksiz*ksiz*ksiz*ksiz*ksiz*ksiz) / ksiz;
                    // std::cout<<" PM dir[i][0] = "<<dir[i][0]<<std::endl;

					//  std::cout<<"gx = "<<gx<<std::endl;
					dir[i][0] += dx*(nagi - magi)/(a * a) + dx*Rx; //direct force
					dir[i][1] += dy*(nagi - magi)/(a * a) + dy*Rx;
					dir[i][2] += dz*(nagi - magi)/(a * a) + dz*Rx;
					 // std::cout<<"Rx = "<<Rx<<std::endl;
					 // std::cout<<"dir[i][0] = "<<dir[i][0]<<std::endl;
					 // std::cout<<"dir[i][1] = "<<dir[i][1]<<std::endl;
					 // std::cout<<"dir[i][2] = "<<dir[i][2]<<std::endl;

				}
                if(distij < 2.)
                {
                    Sh = (0.0122242 + 0.609459*distij - 0.432782 * distij*distij + 0.189694 * distij*distij*distij
                      - 0.0610997 * distij*distij * distij*distij) / distij;
                      dir[i][0]+=Sh*dx;
                      dir[i][1]+=Sh*dy;
                      dir[i][2]+=Sh*dz;
                }
                if(distij <=3. && distij >=2)
                {
                    Sh = (-0.137179 + 0.158876*distij - 0.0348711 *distij*distij) / distij;
                      dir[i][0]+=Sh*dx;
                      dir[i][1]+=Sh*dy;
                      dir[i][2]+=Sh*dz;
                }
		}

	  }
	  // std::cout<<"dir matrix"<<std::endl;
	  // ai::printMatrix(dir);
	  return std::vector<double>{ dir[0][0], dir[0][1], dir[0][2] };
	}


double Signum(double  x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}


void Direct(std::vector<std::vector<double> > & Particles, std::vector<std::vector<double> > & dir, double mass)
{
			double dx ;
			double dy ;
			double dz ;
			double distij;
			double magj, magi;
			double a = 0.3;//равновесное расстояние
			double close; // слагаемое, отвечающее за близкодействие
            //double s = 1500;
			std::vector<std::vector<double> > as;
			as.resize(dir.size());
			for(size_t i =0 ; i< as.size() ; ++i) as[i].resize(3);

	  for(size_t i = 0 ; i<Particles.size(); ++i)
		for(size_t j = i+1 ; j<Particles.size(); ++j)
		{
			dx = Particles[i][0] - Particles[j][0];
			dy = Particles[i][1] - Particles[j][1];
			dz = Particles[i][2] - Particles[j][2];
			distij = sqrt(dx*dx + dy*dy +dz*dz);
			magi = (mass*a*a)/(distij*distij*distij);
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
			dir[j][2] += -mass*mass*(close*dz - magj*dz)/(a*a);
            // dir[j][0] +=  (std::pow(a,13)/std::pow(dx,13) - Signum(dx)*(a*a)/(dx*dx))/(a*a);
            // dir[j][1] +=  (std::pow(a,13)/std::pow(dy,13) - Signum(dy)*(a*a)/(dy*dy))/(a*a);
            // dir[j][2] +=  (std::pow(a,13)/std::pow(dz,13) - Signum(dz)*(a*a)/(dz*dz))/(a*a);

            /*std::cout<<" direct DIST ij = "<<distij<<std::endl;
            std::cout<<" direct dirij x= "<<dir[i][0] <<"    dirji x=  "<<dir[j][0]<<std::endl;
			std::cout<<" direct dirij y= "<<dir[i][1] <<"    dirji y=  "<<dir[j][1]<<std::endl;
			std::cout<<" direct dirij z= "<<dir[i][2] <<"    dirji z=  "<<dir[j][2]<<std::endl;*/

		}
}

// расчет ускорений для 1 частицы вокруг солнца
void DirectSun(std::vector<std::vector<double> > & Particles,
				std::vector<double>& Sun,
					std::vector<std::vector<double> > & dir,
						double mass, double L)
{
	double dx;
	double dy;
	double dz;
	double distij;
	double magj, magi;
	double a = 0.3;//равновесное расстояние
	double close; // слагаемое, отвечающее за близкодействие
				  //double s = 1500;
	double num = Particles.size();

	for (size_t i = 0; i<Particles.size(); ++i)
		//for (size_t j = i + 1; j<Particles.size(); ++j)
		{
			dx = Particles[i][0] - Sun[0];
			dy = Particles[i][1] - Sun[1];
			dz = Particles[i][2] - Sun[2];
			distij = sqrt(dx*dx + dy * dy + dz * dz);
			magi = (a*a) / (distij*distij*distij);
			close = (1. *std::pow(a, 6)) / (std::pow(distij, 7));
			// close =std::exp(-2.0*a*s*(distij - a))/distij;;
			// close =0.;
			dir[i][0] -= - 10.*num*mass * mass*(close*dx - magi * dx) / (a*a);
			dir[i][1] -= - 10.*num*mass * mass*(close*dy - magi * dy) / (a*a);
			dir[i][2] -= - 10.*num*mass * mass*(close*dz - magi * dz) / (a*a);
			// dir[i][0] -=  (std::pow(a,13)/std::pow(dx,13) - Signum(dx)*(a*a)/(dx*dx))/(a*a);
			// dir[i][1] -=  (std::pow(a,13)/std::pow(dy,13) - Signum(dy)*(a*a)/(dy*dy))/(a*a);
			// dir[i][2] -=  (std::pow(a,13)/std::pow(dz,13) - Signum(dz)*(a*a)/(dz*dz))/(a*a);
			/*magj = (1.0*mass*a*a) / (distij*distij*distij);
			dir[j][0] += -mass * mass*(close*dx - magj * dx) / (a*a);
			dir[j][1] += -mass * mass*(close*dy - magj * dy) / (a*a);
			dir[j][2] += -mass * mass*(close* dz - magj * dz) / (a*a);*/
			// dir[j][0] +=  (std::pow(a,13)/std::pow(dx,13) - Signum(dx)*(a*a)/(dx*dx))/(a*a);
			// dir[j][1] +=  (std::pow(a,13)/std::pow(dy,13) - Signum(dy)*(a*a)/(dy*dy))/(a*a);
			// dir[j][2] +=  (std::pow(a,13)/std::pow(dz,13) - Signum(dz)*(a*a)/(dz*dz))/(a*a);

			// std::cout << "direct DIST ij = " << distij << std::endl;
			// std::cout << "direct dirij x= " << dir[i][0] << std::endl;//"    dirji x=  " << dir[j][0] << std::endl;
			// std::cout << "direct dirij y= " << dir[i][1] << std::endl;// "    dirji y=  " << dir[j][1] << std::endl;
			// std::cout << "direct dirij z= " << dir[i][2] << std::endl;// "    dirji z=  " << dir[j][2] << std::endl;

		}
}
