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
	//number of boxes and particle must be N^3, where N=2^n
	size_t dim = 5; //number of BOXes  
	size_t number_particles = 5;
	double mass = 1;
	for (size_t i = 0 ; i < dim; ++i)
		for (size_t j = 0 ; j < dim ;++j)
			for (size_t k = 0; k < dim ;++k){
				BOX box;
				std::vector<double> r = { (double)i*H , j*H , k*H, 0 };
				box.pos = r;
				boxes.push_back(box);
				}

	for (size_t i = 0; i < number_particles; ++i) {
		Particle part; //create Particle 
		std::vector<double> a = { (double)i, 5. - i, 1.7 * i };
		part.r = a;
		particles.push_back(part);
	}
		std::vector<std::vector<double>> Particles = particle_ut::Particles_to_coords(particles); //matrix of particles
		std::vector<std::vector<double>> Boxes = particle_ut::Mesh_to_coords(boxes); // matrix of grid points
		ai::saveA3R("./particles.a3r", Particles);
		ai::saveA3R("./boxes.a3r", Boxes);
		//auto Oldstart = ai::time();
		for (size_t i = 0; i < Particles.size(); ++i){
			for (size_t j = 0; j < Boxes.size(); ++j){
				if (std::abs(Boxes[j][0] - Particles[i][0]) <= 0.5*H &&
					std::abs(Boxes[j][1] - Particles[i][1]) <= 0.5*H &&
					std::abs(Boxes[j][2] - Particles[i][2]) <= 0.5*H ){
					Boxes[j][3] += mass / (H * H * H); //density
					}
			}
		}
		//auto Oldfinish = ai::time();
		//std::cout << "time duration Old = " << ai::duration(Oldstart, Oldfinish, "ms") << " ms" << std::endl;
		std::cout << "Particles= " << std::endl;
		ai::printMatrix(Particles);
		ai::printMatrix(Boxes); //output 
		//ai::saveMatrix("./output/Mesh.txt", Boxes);
		//ai::saveMatrix("./output/Particles.txt", Particles);
		
		std::vector <
			std::vector<
				std::vector<double>>> density;
		//resize 3D matrix
		density.resize(dim);
		for (size_t i = 0; i < dim; ++i){
			density[i].resize(dim);
			for (size_t j = 0; j < dim; ++j)
				density[i][j].resize(dim);
			}
		//now density dimXdimXdim
		//auto NewStart = ai::time();
		double shift = 0.49995; //
		for (size_t i = 0; i < Particles.size(); ++i)
		{
			//calculating indexes
			size_t x, y, z = 0;
			if (std::floor((Particles[i][0] / H) + shift) < dim &&
				std::floor((Particles[i][1] / H) + shift) < dim &&
				std::floor((Particles[i][2] / H) + shift) < dim) {
				x = std::floor((Particles[i][0] / H) + shift);
				y = std::floor((Particles[i][1] / H) + shift);
				z = std::floor((Particles[i][2] / H) + shift);

				density[x][y][z] += mass / (H * H * H);
			}

		}
		//auto NewFinish = ai::time();
		//std::cout << "time duration New = " << ai::duration(NewStart, NewFinish, "ms") << " ms" << std::endl;
		std::cout << "New density  =" << std::endl;
		for (size_t i = 0; i < dim; ++i)
			for (size_t j = 0; j < dim; ++j)
				for (size_t k = 0; k < dim; ++k)
				{
					std::cout <<"("<<i<<" , "<<j<<" , "<<k<<" )"<<" dens="<<density[i][j][k]<<std::endl;
				}
		
		
		return 0;
}
