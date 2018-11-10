#pragma once
#include <vector>
//#include <algorithm>

void ScalePos(std::vector<std::vector<double> >& Particles, double scale);

void  Genconfig(std::vector<std::vector<double> >& Particles, double number_particles,double L);

 void CaclDensity(std::vector<std::vector<double> >&Particles,std::vector<std::vector <std::vector<std::vector<double> > > >&density,
   // std::vector<std::vector<size_t> >& box,
	double mass,
	double H,
	size_t dim,
	double L);
	
	void PuttoBox(std::vector<
					std::vector<double> >& Particles ,
					std::vector<std::vector<size_t> >&box,
					double H,
				 size_t dim
					);


  void CalcPotential(std::vector<
					std::vector<
						std::vector<
							std::vector<double> > > >& rho,
				   size_t dim,
				   double H,
				   double L);
 void GetAccel(std::vector<std::vector<double> >& Particles,
							   std::vector<std::vector<std::vector<std::vector<double > > > >& density,
                  std::vector<std::vector<size_t> >& box,
								         std::vector<std::vector<double> >& a,
									double H);

void Direct(std::vector<std::vector<double> > & Particles, std::vector<std::vector<double> > & dir, double mass);
double Signum(double  x);
