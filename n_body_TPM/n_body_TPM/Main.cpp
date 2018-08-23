#include <vector>
#include "particle_struct.h"
//#include "a3r.h"
#include "ai.hh"
#include "particles_ut.h"



int main() {
	
	std::vector< Particle > particles;
	
	//size_t iterator = 1;
	for (size_t i = 0; i < 5; ++i) {
		Particle part; //create Particle c
		std::vector<double> a = { (double)i, 5. - i, 1.7 * i };
		part.r = a;
		particles.push_back(part);
		//++iterator;
		
	}
	
	
		ai::saveA3R("./test.a3r", particle_ut::Particles_to_coords(particles));
	return 0;
}
