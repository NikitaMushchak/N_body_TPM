#include <vector>
#include <algorithm>
// #include "particle_struct.h"
//#include "a3r.h"
#include "ai.hh"
#include "particles_ut.h"
// #include "box_struct.h"
#include "Fourier_tools.hh"
#include "CaclModule.h"
#include "Integrator.h"
#include "tpm.h"


int main(const int argc, const char *argv[]) {

	// std::vector< Particle > particles;
	// std::vector<BOX> boxes;

	 double H = 1.0; //dimention of the cubic BOX
   //number of boxes and particle must be N^3, where N=2^n
	 double L = 64; // dimention of calculation area
	 double dim = 64; //number of BOXes power of 2
	 double number_particles = 2;
	 double mass = 1.0;
	 double T1 = 0.5;


	 for(int i = 1; i < argc; ++i){
           // if("-v" == std::string(argv[i]) || "--version" == std::string(argv[i])){
           //     std::cout << "  Build: "  << buildVersion << "." << std::endl;
           //     std::cout << "  Compiler: "  << compiler << "." << std::endl;
           //     std::cout << "  Target: "  << target << "." << std::endl;
           //     std::cout << "  AiLibrary: " << ai::getVersion() << "."
           //         << std::endl;
           //     std::cout << "  Compilation timestamp: "  << timestamp << "."
           //         << std::endl;
		   //
           //     return 0;
           // }

           if("-h" == std::string(argv[i]) || "--help" == std::string(argv[i])){
               std::cout << "usage: tpm [options]"
                   << std::endl
                   << "    -h  --help            print this usage and exit"
                   << std::endl
                   << "    -v  --version         print build info and exit"
                   << std::endl
                   << "    --list-errors         print possible errors ans exit"
                   << std::endl << std::endl

                   << "  Program parameters" << std::endl
                   << "    --number_particles=<value>           number of particles"
                   << std::endl
                   << "    --H=<value>          dimention of the cubic BOX"
                   << std::endl
									 << "    --L=<value>          dimention of the calculation area"
                   << std::endl
                   << "    --dim=<value>           number of BOXes power of 2 "
                   << "[double, n/d]"
                   << std::endl << std::endl

                   << "  particles parameters" << std::endl
                   << "    --mass=<value>           mass of partilces "

                   << std::endl << std::endl

                   << "  Time parameters" << std::endl
                   << "    --T1=<value>        modeling time [double]"

                   << std::endl << std::endl;
									 	return 0;
								}
				   if(
				              ai::assignAbsDoubleParameter(argv[i], "--number_particles=", number_particles)
				              || ai::assignAbsDoubleParameter(argv[i], "--H=", H)
											|| ai::assignAbsDoubleParameter(argv[i], "--L=", L)
				              || ai::assignAbsDoubleParameter(argv[i], "--dim=", dim)
				              || ai::assignAbsDoubleParameter(argv[i], "--T1=", T1)
				              || ai::assignAbsDoubleParameter(argv[i], "--mass=", mass)

				              )



				          {
				              continue;
				          }

}
							// auto start = ai::time();
               return TPM(H,L,dim,number_particles,mass, T1);

							 // auto
           }
