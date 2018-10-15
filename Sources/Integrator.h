#pragma once
#include <vector>

inline double Max_num(double f, double m);
inline double Min_num(double f, double m);
inline double Get_Step(std::vector<std::vector<double>> &a, double mass);
void Step_PM(std::vector<std::vector<double>>& Particles, std::vector<std::vector<double>>& vel, 
	std::vector<std::vector<double>>& a, double mass, double dt);
double StepP(std::vector<std::vector<double>>& Particles,
	std::vector<std::vector<double>>& vel,
	std::vector<std::vector<double>>& a,
	double mass, size_t dim, double H);