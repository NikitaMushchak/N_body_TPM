#pragma once
#include <vector>
//#include <algorithm>

void ScalePos(std::vector<std::vector<double>>& Particles, double scale);

 void CaclDensity(std::vector<std::vector<double>>&Particles,std::vector<std::vector <std::vector<std::vector<double>>>>&density,
	double mass,
	double H,
	size_t dim);


  void CalcPotential(std::vector<
					std::vector<
						std::vector<
							std::vector<double>>>>& rho,
				   size_t dim);
 void GetAccelPM(std::vector<std::vector<double>>& Particles,
							std::vector<std::vector<std::vector<std::vector<double >>>>& density,
								std::vector<std::vector<double>>& a,
									double H);
void GetAccel(std::vector<std::vector<double>>& Particles,
           std::vector<std::vector<std::vector<std::vector<double > > > >& density,
                      std::vector<std::vector<size_t> > &box,
            std::vector<std::vector<double> >& a, double H);
			
			
std::vector<double> DirectPM(std::vector<std::vector<double> > & Particles, std::vector<std::vector<double> > & dir, double mass);
 
void Direct(std::vector<std::vector<double>> & Particles, std::vector<std::vector<double>> & dir, double mass);
double Signum(double  x);


void SetSun(std::vector<
 std::vector <
   std::vector<
     std::vector<double> > > >& density,
   double mass,
	size_t dim);
 
 
 void GenRing(std::vector<std::vector<double> >& Particles,
				 std::vector<std::vector<double> >& vel,
				 double number_particles, 
				 double L);
				 
				 
void PuttoBox(std::vector<
					std::vector<double> >& Particles ,
					std::vector<std::vector<size_t> >&box,
					double H,
				 size_t dim);
void DirectSun(std::vector<std::vector<double>> & Particles,
	std::vector<std::vector<double>> & dir,
	double mass, double L);