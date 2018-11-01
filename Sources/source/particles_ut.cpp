
#include "particles_ut.h"
// #include "particle_struct.h"
// #include "box_struct.h"

// #include "ai.hh"
//#include "uitil.h"


// namespace particle_ut
// {
// 	vector <std::vector<double>> Particles_to_coords(vector <Particle>& particles)
// 	{
// 		vector< vector<double> > coords;
// 		vector< Particle > ::iterator p;
// 		for (p = particles.begin(); p < particles.end(); ++p)
// 		{
// 			coords.push_back(p->r);
// 		}
// 		return coords;
// 	}
// 	vector <std::vector<double>> Mesh_to_coords(vector <BOX>& boxes)
// 	{
// 		vector< vector<double> > coords;
// 		vector< BOX > ::iterator p;
// 		for (p = boxes.begin(); p < boxes.end(); ++p)
// 		{
// 			coords.push_back(p->pos);
// 		}
// 		return coords;
// 	}
// }

;
inline void ScalePos(std::vector<std::vector<double>> &Coordinates, double scale)
{
	std::size_t size= Coordinates.size();//number of particles
	for (size_t i = 0 ; i<size ;++i)
	{

		Coordinates[i][0]/=scale;
		Coordinates[i][1]/=scale;
		Coordinates[i][2]/=scale;
	}
}
