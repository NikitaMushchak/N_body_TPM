//#include <vector>
//#include "ai.hh" //input-output library
//#include "particle_struct.h"
//#include "ai.hh"
//#include "particles_ut.h"
//#include "box_struct.h"
//struct Particle
//{
//	std::vector<double> r, v, a;
//	double mass;
//	float radius;
//	char color;
//	Particle() : mass(1.0), radius(1.0), color(0) {
//		std::fill(r.begin(), r.end(), 0.),
//			std::fill(v.begin(), v.end(), 0.),
//			std::fill(a.begin(), a.end(), 0.);
//	}
//};
//struct BOX
//{
//	std::vector<double> pos; // position of center
//	double mass;
//	double H; // dimention of the box
//	double density = mass / (H*H*H);
//	BOX() : mass(1.0), H(1.0), density(0) {
//		std::fill(pos.begin(), pos.end(), 0.);
//	}
//};
//namespace particle_ut
//{
//	vector <std::vector<double>> Particles_to_coords(vector <Particle>& particles)
//	{
//		vector< vector<double> > coords;
//		vector< Particle > ::iterator p;
//		for (p = particles.begin(); p < particles.end(); ++p)
//		{
//			coords.push_back(p->r);
//		}
//		return coords;
//	}
//	vector <std::vector<double>> Mesh_to_coords(vector <BOX>& boxes)
//	{
//		vector< vector<double> > coords;
//		vector< BOX > ::iterator p;
//		for (p = boxes.begin(); p < boxes.end(); ++p)
//		{
//			coords.push_back(p->pos);
//		}
//		return coords;
//	}
//}
//int main() {
//	std::vector< Particle > particles;
//	std::vector<BOX> boxes;
//	double H = 1.0; //dimention of the cubic BOX
//	size_t dim = 5; //number of BOXes  
//	size_t number_particles = 5;
//	double mass = 1;
//	std::vector<double> r_box;
//	for (size_t i = 0; i < dim; ++i)
//		for (size_t j = 0; j < dim; ++j)
//			for (size_t k = 0; k < dim; ++k){
//				BOX box;
//				std::vector<double> r = { (double)i*H , j*H , k*H, 0 };
//				box.pos = r;
//				boxes.push_back(box);
//			}
//
//	for (size_t i = 0; i < number_particles; ++i) {
//		Particle part; //create Particle c
//		std::vector<double> a = { (double)i, 5. - i, 1.7 * i };
//		part.r = a;
//		particles.push_back(part);
//	}
//	std::vector<std::vector<double>> Particles = particle_ut::Particles_to_coords(particles); 
//	std::vector<std::vector<double>> Boxes = particle_ut::Mesh_to_coords(boxes); 
//	std::cout << "Particles" << std::endl;
//	ai::printMatrix(Particles);
//	std::cout << "Mesh" << std::endl;
//	ai::printMatrix(Boxes);
//	ai::saveA3R("./particles.a3r", Particles);
//	ai::saveA3R("./boxes.a3r", Boxes);
//	std::cout << "pa size " << Particles.size() << std::endl;
//	std::cout << "me size " << Boxes.size() << std::endl;
//
//	for (size_t i = 0; i < Particles.size(); ++i)
//	{
//		for (size_t j = 0; j < Boxes.size(); ++j)
//		{
//			if (std::abs(Boxes[j][0] - Particles[i][0]) <= 0.5*H &&
//				std::abs(Boxes[j][1] - Particles[i][1]) <= 0.5*H &&
//				std::abs(Boxes[j][2] - Particles[i][2]) <= 0.5*H){
//				Boxes[j][3] += mass / (H * H * H);
//			}
//		}
//	}
//	ai::saveMatrix("./Boxes.txt", Boxes); //output
//	return 0;
//}