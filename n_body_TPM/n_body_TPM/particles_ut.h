

#include <vector>
//#include "particle_ut.h"


struct Particle;
using std::vector;


namespace particle_ut
{

	vector <std::vector<double>> Particles_to_coords(vector <Particle>& particles);
	vector <float> Particles_to_mass(vector <Particle>& particles);
	vector <float> Particles_to_radius(vector <Particle> & particles);
	vector<char> Particles_to_color(vector <Particle>& particles);


}
