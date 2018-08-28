#include <vector>
#include "particle_struct.h"
//#include "a3r.h"
#include "ai.hh"
#include "particles_ut.h"
#include "box_struct.h"
#include "Fourier_tools.hh"



int main() {
	
	std::vector< Particle > particles;
	std::vector<BOX> boxes;
	double H = 1.0; //dimention of the cubic BOX
	//number of boxes and particle must be N^3, where N=2^n
	size_t dim = 4; //number of BOXes power of 2 
	size_t number_particles = 1;
	double mass = 1;
	for (size_t i = 0 ; i < dim; ++i)
		for (size_t j = 0 ; j < dim ;++j)
			for (size_t k = 0; k < dim ;++k){
				BOX box;
				std::vector<double> r = { (double)i*H , j*H , k*H, 0 };
				box.pos = r;
				boxes.push_back(box);
				}

	for (size_t i = 1; i < number_particles+1; ++i) {
		Particle part; //create Particle 
		std::vector<double> a = { (double)i, 1.9 - i, 1.1*i };
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
		
		std::vector<
			std::vector <
				std::vector<
					std::vector<double>>>> density;
		//resize 3D matrix
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
		//now density dimXdimXdim
		//auto NewStart = ai::time();
		 
		size_t x, y, z = 0;
		double dx, dy, dz, tx, ty, tz;
		for (size_t i = 0; i < Particles.size(); ++i)
		{
			//calculating indexes
			dx = Particles[i][0] - std::floor(Particles[i][0] / H);
			dy = Particles[i][1] - std::floor(Particles[i][1] / H);
			dz = Particles[i][2] - std::floor(Particles[i][2] / H);
			tx = 1. - dx;
			ty = 1. - dy;
			tz = 1. - dz;
			if (std::floor(Particles[i][0] / H) < dim &&
				std::floor(Particles[i][1] / H)  < dim &&
				std::floor(Particles[i][2] / H)  < dim) {
				x = std::floor(Particles[i][0] / H);
				y = std::floor(Particles[i][1] / H);
				z = std::floor(Particles[i][2] / H);

				density[x][y][z][0] += mass * tx *ty*tz;
				
				if (y + 1 >= dim) { y = y % (dim-1); density[x][y][z][0] += mass * tx *dy*tz; }
				else { density[x][y+1][z][0] += mass * tx *dy*tz; }
				if (z + 1 >= dim) { z = z % (dim-1); density[x][y][z][0] += mass * tx *ty*dz; }
				else { density[x][y][z+1][0] += mass * tx *ty*dz; }
				if (y+1>=dim && z + 1 >= dim) { z = z % (dim-1); y = y % (dim-1); density[x][y][z][0] += mass * tx *dy*dz; }
				else{ density[x][y + 1][z + 1][0] += mass * tx *dy*dz; }
				if (x + 1 >= dim) { x = x % (dim - 1); density[x][y][z][0] += mass * dx *ty*tz; }
				else {density[x + 1][y][z][0] += mass * dx *ty*tz; }
				if (x+1>=dim && y+1>=dim) { x = x % (dim - 1); y = y % (dim - 1); density[x][y][z][0] += mass * dx *dy*tz; }
				else { density[x + 1][y + 1][z][0] += mass * dx *dy*tz; }
				if (x+1>=dim && z+1>=dim) { x = x % (dim - 1); z = z % (dim - 1); density[x][y][z][0] += mass * dx *ty*dz; }
				else{ density[x + 1][y][z + 1][0] += mass * dx *ty*dz; }
				if (x+1>=dim && y+1>=dim && z+1>=dim) { x = x % (dim - 1); y =y % (dim - 1); z = z % (dim - 1); density[x][y][z][0] += mass * dx * dy * dz; }
				else { density[x + 1][y + 1][z + 1][0] += mass * dx * dy * dz; }
			
			}

		}
		//auto NewFinish = ai::time();
		//std::cout << "time duration New = " << ai::duration(NewStart, NewFinish, "ms") << " ms" << std::endl;
		std::cout << "New density  =" << std::endl;
		for (size_t i = 0; i < dim; ++i)
			for (size_t j = 0; j < dim; ++j)
				for (size_t k = 0; k < dim; ++k)
				{
					std::cout <<"("<<i<<" , "<<j<<" , "<<k<<" )"<< " dens = " << " ( " << density[i][j][k][0] << " , "
						<< density[i][j][k][1] << " )" <<std::endl;
				}
		//module with Fourier
		

		std::cout << "dens dim" << density.size() << std::endl;

		//FFT for slices i - number of slices
		for (size_t i = 0; i < dim; ++i)
		{
			std::vector <
				std::vector<double>> f(dim); //for rows and colomns
			for (size_t k = 0; k < dim; ++k) { f[k].resize(2); } //resize f as complex
			//std::fill(f.begin(), f.end(), 0);
			
			//FFT rows of density
			for (size_t j = 0; j < dim; ++j) {
				for (size_t k = 0; k < dim; ++k) {
					f[k][0] = density[i][j][k][0];
					f[k][1] = density[i][j][k][1];
				}
				fft(f , dim);
				for (size_t k = 0; k < dim; ++k) {
					density[i][j][k][0] = f[k][0];
					density[i][j][k][1] = f[k][1];
				}
			}
			//FFT coloumns of density
			for (size_t k = 0; k < dim; ++k) {
				for (size_t j = 0; j < dim; ++j) {
					f[j][0] = density[i][j][k][0];
					f[j][1] = density[i][j][k][1];
				}
				fft(f, dim);
				for (size_t j = 0; j < dim; ++j) {
					density[i][j][k][0] = f[j][0];
					density[i][j][k][1] = f[j][1];
				}
			}
		//Solve  equation in Fourier space
		


		//IFFT rows
			for (size_t j = 0; j < dim; ++j) {
				for (size_t k = 0; k < dim; ++k) {
					f[k][0] = density[i][j][k][0];
					f[k][1] = density[i][j][k][1];
				}
				ifft(f, dim);
				for (size_t k = 0; k < dim; ++k)
				{
					density[i][j][k][0] = f[k][0];
					density[i][j][k][1] = f[k][1];					
				}
			}

			//IFFT coloums
			for (size_t k = 0; k < dim; ++k) {
				for (size_t j = 0; j < dim; ++j) {
					f[j][0] = density[i][j][k][0];
					f[j][1] = density[i][j][k][1];
				}
				ifft(f, dim);
				for (size_t j = 0; j < dim; ++j) {
					density[i][j][k][0] = f[j][0];
					density[i][j][k][1] = f[j][1];
				}
			}
		}// cicle for slices

		//output
		std::cout << "FFT and IFFT density  =" << std::endl;
		for (size_t i = 0; i < dim; ++i)
			for (size_t j = 0; j < dim; ++j)
				for (size_t k = 0; k < dim; ++k)
				{
					std::cout << "(" << i << " , " << j << " , " << k << " )" << " dens = " <<" ( " <<density[i][j][k][0]<<" , "
						<<density[i][j][k][1]<<" )"<< std::endl;
				}
		return 0;
}
