#pragma once
#include <vector>
//
//class Particle { 
//public:
//	Particle() : mass(0), radius(0), color(0) {
//		std::fill(r.begin(), r.end(), 0.),
//		std::fill(v.begin(), v.end(), 0.),
//		std::fill(a.begin(), a.end(), 0.);
//	}
//	std::vector<double> get_r( std::vector<double> R)   { this->r   = R; }
//	std::vector<double> get_v( std::vector<double> V)   { this->v   = V; }
//	std::vector<double> get_a( std::vector<double> A)   { this->a   = A; }
//	 
//	const std::vector<double>& R()                    const   { return r; }
//	const std::vector<double>& V()                    const   { return v; }
//	const std::vector<double>& A()                    const   { return a; }
//
//	std::vector<double>& R()							  { return r; }
//	std::vector<double>& V()							  { return v; }
//	std::vector<double>& A()							  { return a; }
//
//	char& Color()						  { return color; }
//	
//	float& Radius()						  { return radius; }
//	float& Mass() { return mass; }
//
//private:
//	float mass;
//	float radius;               //particle radius
//	char color;
//	std::vector<double> r;
//	std::vector<double> v;
//	std::vector<double> a;
//	
// };
struct Particle
{
	std::vector<double> r, v, a;

	double mass;
	float radius;
	char color;
	Particle() : mass(1.0), radius(1.0), color(0) {
		std::fill(r.begin(), r.end(), 0.),
					std::fill(v.begin(), v.end(), 0.),
					std::fill(a.begin(), a.end(), 0.);
	}
 };