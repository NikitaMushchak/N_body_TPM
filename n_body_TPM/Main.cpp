#include <vector>
#include "particle_struct.h"
//#include "a3r.h"
#include "ai.hh"
#include "particles_ut.h"
#include "box_struct.h"
#include "Fourier_tools.hh"
#include "CaclModule.h"


int main() {
	
	std::vector< Particle > particles;
	std::vector<BOX> boxes;
	
	 const double H = 1.0; //dimention of the cubic BOX
	//number of boxes and particle must be N^3, where N=2^n
	size_t dim = 4; //number of BOXes power of 2 
	size_t number_particles = 1;
	double mass = 1.0;
	
	/*for (size_t i = 0 ; i < dim; ++i)
		for (size_t j = 0 ; j < dim ;++j)
			for (size_t k = 0; k < dim ;++k){
				BOX box;
				std::vector<double> r = { (double)i*H , j*H , k*H, 0 };
				box.pos = r;
				boxes.push_back(box);
				}*/

	for (size_t i = 1; i < number_particles+1; ++i) {
		Particle part; //create Particle 
		//std::vector<double> a = { (double) 1.*i,  1.* i, 1.*i };
		std::vector<double> a = { (double) 1.5*i,  1.5* i, 1.5*i };
		part.r = a;
		particles.push_back(part);
	}
		std::vector<std::vector<double>> Particles = particle_ut::Particles_to_coords(particles); //matrix of particles
		//std::vector<std::vector<double>> Boxes = particle_ut::Mesh_to_coords(boxes); // matrix of grid points
		//ai::saveA3R("./particles.a3r", Particles);
		//ai::saveA3R("./boxes.a3r", Boxes);
		//auto Oldstart = ai::time();
		//for (size_t i = 0; i < Particles.size(); ++i){
		//	for (size_t j = 0; j < Boxes.size(); ++j){
		//		if (std::abs(Boxes[j][0] - Particles[i][0]) <= 0.5*H &&
		//			std::abs(Boxes[j][1] - Particles[i][1]) <= 0.5*H &&
		//			std::abs(Boxes[j][2] - Particles[i][2]) <= 0.5*H ){
		//			Boxes[j][3] += mass / (H * H * H); //density
		//			}
		//	}
		//}
		//auto Oldfinish = ai::time();
		//std::cout << "time duration Old = " << ai::duration(Oldstart, Oldfinish, "ms") << " ms" << std::endl;
		//std::cout << "Particles= " << std::endl;
		ai::printMatrix(Particles);
		//ai::printMatrix(Boxes); //output 
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
		//now density dimXdimXdimX2
		
		CaclDensity(Particles, density, mass, H , dim);
		//auto NewFinish = ai::time();
		//std::cout << "time duration New = " << ai::duration(NewStart, NewFinish, "ms") << " ms" << std::endl;
		/*std::cout << "New density  =" << std::endl;
		for (size_t i = 0; i < dim; ++i)
			for (size_t j = 0; j < dim; ++j)
				for (size_t k = 0; k < dim; ++k)
				{
					std::cout <<"("<<i<<" , "<<j<<" , "<<k<<" )"<< " dens = " << " ( " << density[i][j][k][0] << " , "
						<< density[i][j][k][1] << " )" <<std::endl;
				}*/
		//module with Fourier
		

		std::cout << "dens dim = " << density.size() << std::endl;
		auto start = ai::time();
		CalcPotential(density, dim);
		auto finish = ai::time();

		std::cout << "TIME =  " << ai::duration(start, finish, "ms") << "ms" << std::endl;
		//output
		std::cout << "Potential field  =" << std::endl;
		for (size_t i = 0; i < dim; ++i)
			for (size_t j = 0; j < dim; ++j)
				for (size_t k = 0; k < dim; ++k)
				{
					//if (density[i][j][k][0] != 0.)
					std::cout << "(" << i << " , " << j << " , " << k << " )" << " Pot = " <<" ( " <<density[i][j][k][0]<<" , "
						<<density[i][j][k][1]<<" )"<< std::endl;
				}
		//updating particle s velocities and position
		// acceleraton
		
		std::vector <
			std::vector <
				std::vector <
					std::vector<double>>>> g;
		g.resize(dim);
		for (size_t i = 0; i < dim; ++i)
		{
			g[i].resize(dim);
			for (size_t j = 0; j < dim; ++j)
			{
				g[i][j].resize(dim);
				for (size_t k = 0; k < dim; ++k)
				{
					g[i][j][k].resize(3);
				}
			}
		}
		std::cout << "making acceleretion fields ....";
		//acceleration field
		for (size_t i = 1; i < dim-1; ++i)
		{
			for (size_t j = 1; j < dim-1; ++j)
			{
				for (size_t k = 1; k < dim-1; ++k)
				{
					g[i][j][k][0] = -0.5*(density[i + 1][j][k][0] - density[i - 1][j][k][0]);
					g[i][j][k][1] = -0.5*(density[i][j+1][k][0] - density[i][j-1][k][0]);
					g[i][j][k][2] = -0.5*(density[i][j][k+1][0] - density[i][j][k-1][0]);
				}
			}
		}
		std::cout << "Done." << std::endl;
		//accleration of particles
		std::vector<
			std::vector<double>> a;
		a.resize(Particles.size());
		for (size_t i = 0; i < Particles.size(); ++i)
		{
			a[i].resize(3);
		}
		//particles cicle
		//double dx, dy, dz, tx, ty, tz;
		std::cout << "Interpolating forces into particles ";
		size_t x, y, z = 0;
		for (size_t i = 0; i < Particles.size(); ++i)
		{
			//calculating indexes of parent cell
			double dx = Particles[i][0] - std::floor((Particles[i][0]) / H);
			double dy = Particles[i][1] - std::floor((Particles[i][1]) / H);
			double dz = Particles[i][2] - std::floor((Particles[i][2]) / H);
			double tx = 1. - dx;
			double ty = 1. - dy;
			double tz = 1. - dz;
			if (std::floor(Particles[i][0] / H) < dim &&
				std::floor(Particles[i][1] / H) < dim &&
				std::floor(Particles[i][2] / H) < dim) {
				x = std::floor(Particles[i][0] / H);
				y = std::floor(Particles[i][1] / H);
				z = std::floor(Particles[i][2] / H);
			}
			std::cout << ".";
			a[i][0] = g[x][y][z][0] * tx*ty*tz + g[x + 1][y][z][0] * dx*ty*tz +
				g[x][y + 1][z][0] * tx*dy*tz + g[x + 1][y + 1][z][0] * dx*dy*tz +
				g[x][y][z + 1][0] * tx*ty*dz + g[x + 1][y][z + 1][0] * dx*ty*dz +
				g[x][y + 1][z + 1][0] * tx*dy*dz + g[x + 1][y + 1][z + 1][0] * dx*dy*dz;
			std::cout << ".";
			a[i][1] = g[x][y][z][1] * tx*ty*tz + g[x + 1][y][z][1] * dx*ty*tz +
						g[x][y + 1][z][1] * tx*dy*tz + g[x + 1][y + 1][z][1] * dx*dy*tz +
							g[x][y][z + 1][1] * tx*ty*dz + g[x + 1][y][z + 1][1] * dx*ty*dz +
								g[x][y + 1][z + 1][1] * tx*dy*dz + g[x + 1][y + 1][z + 1][1] *dx*dy*dz;
			std::cout << ".";
			a[i][2] = g[x][y][z][2] * tx*ty*tz + g[x + 1][y][z][2] * dx*ty*tz +
						g[x][y + 1][z][2] * tx*dy*tz + g[x + 1][y + 1][z][2] * dx*dy*tz +
							g[x][y][z + 1][2] * tx*ty*dz + g[x + 1][y][z + 1][2] * dx*ty*dz +
								g[x][y + 1][z + 1][2] * tx*dy*dz + g[x + 1][y + 1][z + 1][2] * dx*dy*dz;

		}
		std::cout << " DONE."<<std::endl;
		//	std::vector<
		//	std::vector<double>> g;
		//g.resize(Particles.size());
		//for (size_t i = 0; i < Particles.size(); ++i) {
		//	g[i].resize(3); // x , y, z component of acceleraton
		//}
		//std::cout << "g size = " << g.size() << std::endl;
		//size_t x, y, z = 0;
		//double dx, dy, dz, tx, ty, tz;
		//std::cout << "Calculate accelerations .....";
		//for (size_t i = 0; i < Particles.size(); ++i)
		//{
		//	//std::fill(g.begin(), g.end(), 0.);
		//	for (size_t v = 0; v < g.size(); ++v)
		//	{
		//		g[v][0] = 0.;
		//		g[v][1] = 0.;
		//		g[v][2] = 0.;
		//	}
		//	
		//	//calculating indexes
		//	dx = Particles[i][0] - std::floor(Particles[i][0]);/// H);
		//	dy = Particles[i][1] - std::floor(Particles[i][1]);/// H);
		//	dz = Particles[i][2] - std::floor(Particles[i][2]);/// H);
		//	tx = 1. - dx;
		//	ty = 1. - dy;
		//	tz = 1. - dz;
		//	if (std::floor(Particles[i][0] / H) < dim &&
		//		std::floor(Particles[i][1] / H) < dim &&
		//		std::floor(Particles[i][2] / H) < dim) {
		//		x = std::floor(Particles[i][0] / H);
		//		y = std::floor(Particles[i][1] / H);
		//		z = std::floor(Particles[i][2] / H);

		//		std::cout << "x , y , z = " << x << " " << y << " " << z << std::endl;
		//	}
		//	
		//	// calculate g at particles
		//	// x componet
		//	g[i][0] = -0.5*(density[x + 1][y][z][0] - density[x - 1][y][z][0])*tx*ty*tz- 0.5*(density[x + 2][y][z][0] - density[x][y][z][0])*dx*ty*tz
		//		- 0.5*(density[x + 1][y + 1][z][0] - density[x - 1][y + 1][z][0])*tx*dy*tz - 0.5*(density[x + 2][y + 1][z][0] - density[x][y + 1][z][0])*dx*dy*tz
		//		- 0.5*(density[x + 1][y][z + 1][0] - density[x - 1][y][z + 1][0])*tx*ty*dz - 0.5 *(density[x + 2][y][z + 1][0] - density[x][y][z + 1][0])*dx*ty*dz
		//		- 0.5*(density[x + 1][y + 1][z + 1][0] - density[x - 1][y + 1][z + 1][0])*tx*dy*dz - 0.5*(density[x + 2][y + 1][z + 1][0] - density[x][y + 1][z + 1][0])*dx*dy*dz;
		//	std::cout << "g[i].x = " << g[i][0] << std::endl;
		//	//y component
		//	g[i][1] = -0.5*(density[x][y + 1][z][0] - density[x][y - 1][z][0])*tx*ty*tz - 0.5*(density[x + 1][y+1][z][0] - density[x+1][y-1][z][0])*dx*ty*tz
		//		- 0.5*(density[x][y + 2][z][0] - density[x][y][z][0])*tx*dy*tz - 0.5*(density[x + 1][y + 2][z][0] - density[x+1][y][z][0])*dx*dy*tz
		//		- 0.5*(density[x][y+1][z + 1][0] - density[x][y-1][z + 1][0])*tx*ty*dz - 0.5 *(density[x + 1][y+1][z + 1][0] - density[x+1][y-1][z + 1][0])*dx*ty*dz
		//		- 0.5*(density[x][y + 2][z + 1][0] - density[x][y][z + 1][0])*tx*dy*dz - 0.5*(density[x + 1][y + 2][z + 1][0] - density[x+1][y][z + 1][0])*dx*dy*dz;
		//	std::cout << "g[i].y = " << g[i][1] << std::endl;
		//	//z component
		//	g[i][2] = -0.5*(density[x][y][z + 1][0] - density[x][y][z - 1][0])*tx*ty*tz -0.5*(density[x + 1][y][z + 1][0] - density[x + 1][y][z - 1][0])*dx*ty*tz
		//		- 0.5*(density[x][y + 1][z+1][0] - density[x][y+1][z-1][0])*tx*dy*tz - 0.5*(density[x + 1][y + 1][z+1][0] - density[x + 1][y+1][z-1][0])*dx*dy*tz
		//		- 0.5*(density[x][y][z + 2][0] - density[x][y][z][0])*tx*ty*dz - 0.5 *(density[x + 1][y][z + 2][0] - density[x + 1][y][z][0])*dx*ty*dz
		//		- 0.5*(density[x][y + 1][z + 2][0] - density[x][y+1][z][0])*tx*dy*dz - 0.5*(density[x + 1][y + 1][z + 2][0] - density[x + 1][y+1][z][0])*dx*dy*dz;
		//	std::cout << "g[i].z = " << g[i][2] << std::endl;
		//}
		//std::cout << "DONE!" << std::endl;
		for (size_t i = 0; i < Particles.size(); i++)
		{
			std::cout << "Acceleration particle # << i<< = " << " ( " << a[i][0] << " , " << a[i][1] << " , " << a[i][2] << " )" << std::endl;
		}

		//integrator here

		return 0;
}
