#include <vector>
#include "particle_struct.h"
//#include "a3r.h"
#include "ai.hh"
#include "particles_ut.h"
#include "box_struct.h"



int main() {
	
	std::vector< Particle > particles;
	std::vector<BOX> boxes;
	double H = 1.0; //dimention of the cubic BOX
	size_t dim = 5; //number of BOXes  
	size_t number_particles = 5;
	double mass = 1;
	std::vector<double> r_box;
	for (size_t i = 0 ; i < dim; ++i)
		for (size_t j = 0 ; j < dim ;++j)
			for (size_t k = 0; k < dim ;++k)
			{
				BOX box;
				std::vector<double> r = { (double)i*H , j*H , k*H, 0 };
				box.pos = r;
				boxes.push_back(box);
				std::cout << "Box center pos = " << "(" << box.pos[0] << "," << box.pos[1] << "," << box.pos[2] << ") " 
					<<"  initial density = "<<box.pos[3]<< std::endl;
			}
	
	//std::cout << "Box center pos = " << "(" << boxes[0][0] << "," << box.pos[1] << "," << box.pos[2] << ") " << std::endl;
	
	//size_t iterator = 1;
	for (size_t i = 0; i < number_particles; ++i) {
		Particle part; //create Particle c
		std::vector<double> a = { (double)i, 5. - i, 1.7 * i };
		part.r = a;
		particles.push_back(part);
		//++iterator;
		std::cout << "particle pos = " << "(" << part.r[0] << "," << part.r[1] << "," << part.r[2] << ") " << std::endl;
	}
		std::vector<std::vector<double>> Particles = particle_ut::Particles_to_coords(particles); //delete particles
		std::vector<std::vector<double>> Boxes = particle_ut::Mesh_to_coords(boxes); //delete boxes
		std::cout << "Particles" << std::endl;
		ai::printMatrix(Particles);
		std::cout << "Mesh" << std::endl;
		ai::printMatrix(Boxes);
		ai::saveA3R("./particles.a3r", Particles);
		ai::saveA3R("./boxes.a3r", Boxes);
		std::cout << "pa size " << Particles.size() << std::endl;
		std::cout << "me size " << Boxes.size() << std::endl;
		
		for (size_t i = 0; i < Particles.size(); ++i)
		{
			for (size_t j = 0; j < Boxes.size(); ++j)
			{
				if (
					std::abs(Boxes[j][0] - Particles[i][0]) <= 0.5*H &&
					std::abs(Boxes[j][1] - Particles[i][1]) <= 0.5*H &&
					std::abs(Boxes[j][2] - Particles[i][2]) <= 0.5*H 
					)
				{
					Boxes[j][3] += mass / (H * H * H);
					std::cout << "Number particle" << i;
				}
			}
		
		}
		std::cout << "New Boxes" << std::endl;
		ai::printMatrix(Boxes);
		ai::saveMatrix("./Boxes.txt",Boxes);
		return 0;
}
