#pragma once
#include <vector>

template <typename T>
class Particle {
public:

private:
	double mass;
	double r;               //radius 
	std::vector<double> pos, vel, accel;
	double x, y, z;
	pos{ x,y z };  //coordinates of particle
	double vx, vy, vz;
	vel{ vx,vy,vz }; //velocities of particle
	double ax, ay, az;
	accel{ ax,ay,az }; //acceleration of particle
 };
