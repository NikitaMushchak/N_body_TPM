#include <vector>
#include <iostream>
#include "tpm.h"
#include "CaclModule.h"
#include "Integrator.h"
#include "ai.hh"


int TPM(const double H, const double L , const double dim,const double number_particles,const double mass,const double T1){

std::cout<<" _ _       _           _      " <<std::endl;
std::cout         <<"| \\ | ___ | |_  ___  _| | _ _ "<<std::endl;
std::cout         <<"|   ||___|| . \\/ . \\/ . || | |"<<std::endl;
std::cout         <<"|_\\_|     |___/\\___/\\___|`_. |"<<std::endl;
std::cout         <<"                         <___'"<<std::endl;
std::cout<<" _______ _____  __  __ "<<std::endl;
std::cout        <<"|__   __|  __ \\|  \\/  |"<<std::endl;
std::cout        <<"   | |  | |__) | \\  / |"<<std::endl;
std::cout        <<"   | |  |  ___/| |\\/| |"<<std::endl;
std::cout        <<"   | |  | |    | |  | |"<<std::endl;
std::cout        <<"   |_|  |_|    |_|  |_|"<<std::endl;
std::cout<<"--Developed by Nikita Mushchak, 2018--"<<std::endl;
std::cout<<"Department of Theoretical Mechanics, SPbSTU, Russia"<<std::endl;


std::cout << "Parameters: "<<std::endl
          <<" H = "<<H<<std::endl
          <<" L = "<<L<<std::endl
          <<" dim = "<<dim <<std::endl
          <<"number_particles = "<<number_particles<<std::endl
          <<"mass = "<<mass<<std::endl
          <<"T1 = "<<T1<<std::endl;

std::vector<std::vector<double>> Particles;
for (size_t i = 0; i < number_particles ; ++i) {
    std::vector<double> a = { (double) 16.+1.9*i,  16.+1.9*i, 16.+1.9*i };
    Particles.push_back(a);
}
std::cout << "non-Scaled pos"<<std::endl;
ai::printMatrix(Particles);

double scale = L/dim;

ScalePos(Particles,scale);

std::cout << "Scaled pos"<<std::endl;
ai::printMatrix(Particles);
// const double T1 = 3.1;
  double dt=0;
  double time=0;
 size_t it = 0 ;


 std::vector<
     std::vector<
     std::vector<
     std::vector<double> > > > density;
 //resize 3D matrix

 std::vector<std::vector<double> >a;


 std::vector<std::vector<double> >dir;


 std::vector<std::vector<double> > vel;
 /**/

//integrator here
auto start = ai::time();
while (time <T1)
    {
    it++;
    auto st = ai::time();
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

    a.resize(Particles.size());
    for (size_t i = 0; i < Particles.size(); i++)
    {
    a[i].resize(3);
    }
    dir.resize(Particles.size());
    for (size_t i = 0; i < Particles.size(); i++)
    {
    dir[i].resize(3);
    }

    vel.resize(Particles.size());
    for (size_t i = 0; i < Particles.size(); i++)
    {
    vel[i].resize(3);
    }
    auto en = ai::time();
    std::cout<<"Matrix time = "<<ai::duration(st, en ,"ms")<<" ms"<<std::endl;
     // std::cout<<"Particles coords"<<std::endl;
     // ai::printMatrix(Particles);

     //CIC assigment
     auto t1 =ai::time();
     CaclDensity(Particles, density, mass, H, dim);
     auto t2=ai::time();
     std::cout<<"Calcdensity time = "<<ai::duration(t1, t2 , "ms")<<" ms"<<std::endl;
     //potetial field
     auto t3 = ai::time();
     CalcPotential(density, dim);
     auto t4 =ai::time();
     std::cout<<"CalcPotential time = "<<ai::duration(t3,t4, "ms")<<" ms"<<std::endl;
     //Acceleration

     auto t5 = ai::time();
     GetAccel(Particles, density, a, H);
     auto t6 = ai::time();
     std::cout <<"Acclel time = "<<ai::duration(t5,t6,"ms")<<" ms"<<std::endl;

     auto t7 = ai::time();

     Direct(Particles, dir, mass);
     auto t8 = ai::time();
     std::cout <<"Direct time= "<<ai::duration(t7,t8,"ms")<<" ms"<<std::endl;
     std::cout<<"Dir accel "<<std::endl;
     ai::printMatrix(dir);

     std::cout<<"Acceleration"<<std::endl;
     ai::printMatrix(a);
     std::cout<<"  \nDir/ FFT"<<dir[0][0]/a[0][0]<<std::endl;

     auto t10 = ai::time();
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
     auto t11 = ai::time();
     std::cout <<"update part = "<<ai::duration(t10,t11, "ms")<<" ms"<<std::endl;
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
     // vel.clear();
     // a.clear();

     // std::fill(Particles.begin(), Particles.end(), 0.);
     // std::fill(density.begin(), density.end(), 0.);
     // std::fill(a.begin(), a.end(), 0.);
     // std::fill(dir.begin(), dir.end(), 0.);
     // std::fill(vel.begin(), vel.end(), 0);

    }

    auto finish = ai::time();
    std::cout <<"Time caclulation = "<<ai::duration (start , finish, "ms")<<" ms"<<std::endl;
return 0;
}
