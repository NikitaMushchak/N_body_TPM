#include <vector>
#include <iostream>
#include "tpm.h"
#include "CaclModule.h"
#include "Integrator.h"
#include "ai.hh"

int TPM(const double H, const double dim,const double number_particles,const double mass,const double T1){

std::vector<std::vector<double>> Particles;
for (size_t i = 0; i < number_particles ; ++i) {
    std::vector<double> a = { (double) 1.+3.1*i,  1.+3.1*i, 1.+3.1*i  };
    Particles.push_back(a);
}


// const double T1 = 3.1;
  double dt=0;
  double time=0;
 size_t it = 0 ;
//integrator here
while (time <T1)
    {
    it++;
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
     std::vector<std::vector<double >>a;
     a.resize(Particles.size());
     for (size_t i = 0; i < Particles.size(); i++)
     {
     a[i].resize(3);
     }
     GetAccel(Particles, density, a, H);

     std::vector<std::vector<double >>dir;
     dir.resize(Particles.size());
     for (size_t i = 0; i < Particles.size(); i++)
     {
     dir[i].resize(3);
     }
     Direct(Particles, dir, mass);
     std::cout<<"Dir accel "<<std::endl;
     ai::printMatrix(dir);
     density.clear();
     std::cout<<"Acceleration"<<std::endl;
     ai::printMatrix(a);

     std::vector<std::vector<double >> vel;
     vel.resize(Particles.size());
     for (size_t i = 0; i < Particles.size(); i++)
     {
     vel[i].resize(3);
     }/**/

     dt =  10*Get_Step(dir, mass);
     //std::cout << "dt = " << dt << std::endl;
     // time += dt;
     for (size_t i = 0; i < Particles.size(); ++i)
     {
         vel[i][0] += dir[i][0] * dt;
         vel[i][1] += dir[i][1] * dt;
         vel[i][2] += dir[i][2] * dt;
     }
     // std::cout << "Vel\n";
     // ai::printMatrix(vel);
     for (size_t i = 0; i < Particles.size(); ++i)
     {
         Particles[i][0] += vel[i][0] * dt;
         Particles[i][1] += vel[i][1] * dt;
         Particles[i][2] += vel[i][2] * dt;
     }
     // std::cout << "Get_Step" << dt << std::endl;
     //Step_PM(Particles, vel, a, mass, dt);
     time += dt;

     std::cout << "Step_PM # = " << it <<"      dt = "<<dt<<"   time = "<<time  <<std::endl;
     // std::cout<<"Particles pos = "<<std::endl;
     // ai::printMatrix(Particles);
     char cbuf[512];
     std::string suff = "p.a3r";
     sprintf(cbuf, "./Results/%09d_", it);
     std::string filename = cbuf + suff;

     //print results
     ai::saveA3R(filename, Particles);
     //ai::saveMatrix(filename, Particles);
     vel.clear();
     a.clear();
    }
return 0;
}
