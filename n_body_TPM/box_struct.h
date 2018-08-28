#pragma once
#include <vector>


struct BOX
{
	std::vector<double> pos; // position of center

	double mass;
	
	double H; // dimention of the box
	double density = mass / (H*H*H);

	BOX() : mass(1.0), H(1.0), density(0) {
		std::fill(pos.begin(), pos.end(), 0.);
			}
};