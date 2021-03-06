#pragma once
#include <vector>

void CaclDensity(std::vector<
					std::vector<double>>&Particles,
				 std::vector<
						std::vector <
							std::vector<
								std::vector<double>>>>& density,
				 float mass,
	             double H, 
				 size_t dim);


void CalcPotential(std::vector<
					std::vector<
						std::vector<
							std::vector<double>>>>& rho,
				   size_t dim);