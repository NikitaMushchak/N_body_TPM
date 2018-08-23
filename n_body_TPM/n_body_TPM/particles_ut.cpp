#include "particles_ut.h"
#include "particle_struct.h"

#include "ai.hh"
#include "uitil.h"


namespace particle_ut
{
	vector <std::vector<double>> Particles_to_coords(vector <Particle>& particles)
	{
		vector< vector<double> > coords;
		vector< Particle > ::iterator p;
		for (p = particles.begin(); p < particles.end(); ++p)
		{
			coords.push_back(p->r);
		}
		return coords;
	}
}
