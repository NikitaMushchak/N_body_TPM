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

	std::vector<std::vector<double> > Particles;
	Particles.resize(number_particles);

	for(size_t i = 0 ; i< number_particles; ++i )
	{
		Particles[i].resize(3);
	}



	 for (size_t i = 0; i < number_particles ; ++i) {
	 	std::vector<double> a = { (double) 32.- 13./2. + 13.*i, 32. , 32. };
	 	Particles[i] = a;
	 }
// 2 - 35% 3 - 38% 3.1 - 37%
	std::vector<std::vector<double> > PM;

	std::vector<std::vector<double> >  APM;

	std::vector<std::vector<double> > Dir;

	std::vector<std::vector<double> > sh;

	std::vector<double> Sun = {L/2. , L/2. , L/2.};

	// std::cout << "non-Scaled pos"<<std::endl;
	// ai::printMatrix(Particles);

	double scale = L/dim;
	//std::cout<<"scale = "<<scale<<std::endl;

	ScalePos(Particles,scale);

	 //std::cout << "Scaled pos"<<std::endl;
	 //ai::printMatrix(Particles);
	  double dt=0;
	  double time=0;
	 size_t it = 0 ;


	 std::vector<
		 std::vector<
		 std::vector<
		 std::vector<double> > > > density;
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


		std::vector<
		 std::vector<
		 std::vector<
		 std::vector<double> > > > null4;
		 null4.resize(dim);
		for (size_t i = 0; i < dim; ++i) {
			null4[i].resize(dim);
			for (size_t j = 0; j < dim; ++j) {
				null4[i][j].resize(dim);
				for (size_t k = 0; k < dim; ++k) {
					null4[i][j][k].resize(2);
				}
			}
		}

	// std::vector<std::vector<size_t> > nuls; //for particles in same box
	//
	//
	// 	  nuls.resize(Particles.size());
	//
	// 	  for (size_t i =0 ; i<Particles.size(); ++i)
	//  {
	//
	// 		 nuls[i].resize(4);
	//
	//  }



	 std::vector<std::vector<double> >a;
	a.resize(Particles.size());
		for (size_t i = 0; i < Particles.size(); i++)
		{
			a[i].resize(3);
		}


		std::vector<std::vector<double> >apm;
		apm.resize(Particles.size());
		   for (size_t i = 0; i < Particles.size(); i++)
		   {
		   		apm[i].resize(3);
		   }
	  std::vector<std::vector<double> >dir;
	 dir.resize(Particles.size());
		 for (size_t i = 0; i < Particles.size(); i++)
		 {
		 	dir[i].resize(3);
		 }

	 std::vector<std::vector<double> > vel;
	  vel.resize(Particles.size());
		for (size_t i = 0; i < Particles.size(); i++)
		{
			vel[i].resize(3);
		}

		std::vector< std::vector<double > > null2 ;
		null2.resize( Particles.size());
		for (size_t i = 0 ; i<  Particles.size() ; ++i)
		{
			null2[i].resize(3);
		}
	 /**/

	 //Функция генерируюшая кольцо частиц
	 // GenRing(Particles,  vel,  number_particles,  L);

	double r = 3.0; //радиус близкодействия
	double energy;
	// TODO вынести константу!!!!
	//integrator here
	auto start = ai::time();
	while (time < T1)
		{
		it++;

		 energy=0.;

		//Добавление в центр сетки гравитирующенго тела

		// CaclDensitySun(Sun ,
		//                  density,
		//                  mass,
		//                 H,
		//                  dim);
		 //CIC assigment
		 auto t1 =ai::time();
		 CaclDensity(Particles, density, mass, scale, dim );
		 auto t2=ai::time();
		 std::cout<<"Calcdensity time = "<<ai::duration(t1, t2 , "ms")<<" ms"<<std::endl;
		 //potetial field


		 auto t3 = ai::time();

		 CalcPotential(density, dim);

		 auto t4 =ai::time();
		 std::cout<<"CalcPotential time = "<<ai::duration(t3,t4, "ms")<<" ms"<<std::endl;
		 //Acceleration

		 auto t5 = ai::time();
		 //Расчет ускорений по сеточному методу
		 GetAccelPM(Particles, density, a, scale);
		 apm=a;
		  // Direct(Particles , a , mass);
		 GetAccel( Particles, a, H, r, dim);//рачет ускорений по методу PPPM



		 //расчет ускорений по прямому алгоритму
		// DirectSun(Particles, a, mass, L);
		 Direct(Particles, dir, mass);

		  PM.push_back( std::vector<double>{std::abs(Particles[0][0]-Particles[1][0]) , a[0][0]});
		  APM.push_back( std::vector<double>{std::abs(Particles[0][0]-Particles[1][0]) , apm[0][0]});
		  Dir.push_back(std::vector<double>{std::abs(Particles[0][0]-Particles[1][0]) , dir[0][0]});
		  //Dir.push_back(std::vector<double>{Particles[0][0], Particles[0][1], Particles[0][2] , dir[0][0], dir[0][1], dir[0][2]});


		 //std::cout<<"PM accel"<<std::endl;
		 //ai::printMatrix(a);


		 std::cout<<"Dir accel"<<std::endl;
		 //ai::printMatrix(dir);
		 auto t6 = ai::time();
		 std::cout <<"Acclel time = "<<ai::duration(t5,t6,"ms")<<" ms"<<std::endl;

		 auto t7 = ai::time();

		 //Direct(Particles, dir, mass);
		 auto t8 = ai::time();
		 std::cout <<"Direct time= "<<ai::duration(t7,t8,"ms")<<" ms"<<std::endl;
		 //std::cout<<"Dir accel "<<std::endl;
		 //ai::printMatrix(dir);

		 //std::cout<<"Acceleration"<<std::endl;
		 //ai::printMatrix(a);
		 std::cout<<"  \nDir / FFT  = "<< std::abs(1-(dir[0][0]/a[0][0]))*100 <<" %"<<std::endl;

		 auto t10 = ai::time();
		 dt =  1*Get_Step(a, mass);
		//dt=0.0001;
		 //std::cout << "dt = " << dt << std::endl;
		 // time += dt;
		 for (size_t i = 0; i < Particles.size(); ++i)
		 {
			 vel[i][0] += a[i][0] * dt;
			 vel[i][1] += a[i][1] * dt;
			 vel[i][2] += a[i][2] * dt;

			 //energy += std::sqrt( vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2] );
		 }
		 sh.push_back(std::vector<double>{it , std::abs(Particles[0][0]-Particles[1][0]) ,dir[0][0]/a[0][0] , energy});
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
		 if (it % 1 == 0 || it==1)
		 {
			 char cbuf[512];
			 std::string suff = "p.a3r";
			 sprintf(cbuf, "./Results/%09d_", it);
			 std::string filename = cbuf + suff;

			 //print results
			 ai::saveA3R(filename, Particles , L);

		 }
		  auto st = ai::time();

		  a=null2;
		  dir=null2;
		  //vel=null2;
		  density=null4;
		 //box=nuls;

		  auto en = ai::time();
		std::cout<<"Matrix time = "<<ai::duration(st, en ,"ms")<<" ms"<<std::endl;

		}

		auto finish = ai::time();
		std::cout <<"Time caclulation = "<<ai::duration (start , finish, "ms")<<" ms"<<std::endl;

		std::cout<<"Saving results ..."<<std::endl;
		ai::saveMatrix("./dir", Dir);
		ai::saveMatrix("./pm", PM);
		ai::saveMatrix("./apm", APM);
		ai::saveMatrix("./sh",sh);

		std::cout<<"Done"<<std::endl;

		/*int iy;
		std::cin >> iy;*/
	return 0;
	}
