#pragma once
#include <vector>

template <typename T>
class Particle { 
public:
	Particle() {
		std::fill(pos.begin(),     pos.end(), 0.);
		std::fill(vel.begin(),     vel.end(), 0.);
		std::fill(accel.begin(), accel.end(), 0.);
		float mass = 0;
		float r = 0;
	}
	std::vector<double> get_pos(std::vector<double> pos)   { this->pos   = pos; }
	std::vector<double> get_vel(std::vector<double> vel)   { this->vel   = vel; }
	std::vector<double> get_pos(std::vector<double> accel) { this->accel = accel; }

	d

private:
	float mass;
	float r;               //particle radius
	std::vector<double> pos, vel, accel;
	double x, y, z;
	pos{ x ,y, z };  //coordinates of particle
	double vx, vy, vz;
	vel{ vx,vy,vz }; //velocities of particle
	double ax, ay, az;
	accel{ ax,ay,az }; //acceleration of particle
 };
