#include <vector>
#include "particle_struct.h"
#include "Integrator.h"
#include "CaclModule.h"
#include "ai.hh"
//input forces

inline double Max_num(double f, double m)
{
	return (f < m) ? m :f;
}
inline double Min_num(double f, double m)
{
	return (f < m) ? f : m ;
}
inline double Get_Step(std::vector<std::vector<double>> &a, double mass)
{
	double  u_max_a0_ratio = 100.;// maximum possible dv of particles
	double f_max = 0.;
	double a0 = 1.0;
		double gamma0 = 1.;
	double PI = 3.14159265358979;
	double sqrt11 = 3.316624790355399849;
	double dt_max_ratio = 0.01;
	double t0 = 2 * PI / sqrt(gamma0 * mass / a0 / a0 / a0) /sqrt11;
	double dt_max = dt_max_ratio * t0;
	 for (size_t i = 0; i < a.size(); i++)
	 {
		 f_max = Max_num(f_max, (a[i][0] * a[i][0] + a[i][1] * a[i][1] + a[i][2] * a[i][2]));
	 }
	 return  Min_num(dt_max, sqrt((a0  * u_max_a0_ratio) / f_max));
}
 void Step_PM(std::vector<std::vector<double>>& Particles , std::vector<std::vector<double>>& vel , std::vector<std::vector<double>>& a, double mass, double dt)
{
	//std::cout << "Par = " << Particles.size() << "  vel= " << vel.size() << "   a  = " << a.size() << std::endl;
	//double dt = 0;
	
	for (size_t i = 0; i < Particles.size(); ++i)
	{
		vel[i][0] += a[i][0] * dt;
		vel[i][1] += a[i][1] * dt;
		vel[i][2] += a[i][2] * dt;
	}

	for (size_t i = 0; i < Particles.size(); ++i)
	{
		Particles[i][0] += vel[i][0] * dt;
		Particles[i][1] += vel[i][1] * dt;
		Particles[i][2] += vel[i][2] * dt;
	}
	ai::printMatrix(Particles);
	ai::printMatrix(a);

	//return dt;
}
 double StepP(std::vector<std::vector<double>>& Particles,
	std::vector<std::vector<double>>& vel,
	std::vector<std::vector<double>>& a,
	double mass, size_t dim, double H )
{
	std::vector<
		std::vector<
		std::vector<
		std::vector<double>>>> density;
	//resize 3D matrix
	density.resize(dim);
	for (size_t i = 0; i < dim; ++i) {
		density[i].resize(dim);
		for (size_t j = 0; j < dim; ++j) {
			density[i][j].resize(dim);
			for (size_t k = 0; k < dim; ++k) {
				density[i][j][k].resize(2);
			}
		}
	}
	//CIC assigment
	CaclDensity(Particles, density, mass, H, dim);
	//potetial field
	CalcPotential(density, dim);
	//Acceleration
	/*std::vector<std::vector<double >>a;
	a.resize(Particles.size());
	for (size_t i = 0; i < Particles.size(); i++)
	{
		a[i].resize(3);
		
	}*/
	 GetAccel(Particles, density,a, H);
	density.clear();

	/*std::vector<std::vector<double >> vel;
	vel.resize(Particles.size());
	for (size_t i = 0; i < Particles.size(); i++)
	{
		vel[i].resize(3);
	}*/
	double dt;
	double time;
	dt = Get_Step(a, mass);
	for (size_t i = 0; i < Particles.size(); ++i)
	{
		vel[i][0] += a[i][0] * dt;
		vel[i][1] += a[i][1] * dt;
		vel[i][2] += a[i][2] * dt;
	}

	for (size_t i = 0; i < Particles.size(); ++i)
	{
		Particles[i][0] += vel[i][0] * dt;
		Particles[i][1] += vel[i][1] * dt;
		Particles[i][2] += vel[i][2] * dt;
	}
	std::cout << "Get_Step" << dt << std::endl;
	//Step_PM(Particles, vel, a, mass, dt);
	  time = dt;
	//double t += dt;
	 std::cout << "Step_PM = " << dt << std::endl;
	 char cbuf[256];
	std::string suff = "p.a3r";
	sprintf(cbuf, "./Results/%09d_", dt);
	std::string filename = cbuf + suff;
	
	//print results
	ai::saveA3R(filename , Particles);
	ai::saveMatrix(filename, Particles);
	return dt;
}