#pragma once
#include <limits>
#include "ContDynamicalSystem.h"

struct rk45opts {
	double epsilon = 1e-5;
	double h0=0.2; 
	double hMax= std::numeric_limits<double>::infinity(); 
	double hMin= 1e-14; 
	int maxSteps = 1000000;

};

class rk45
{
public: 
	explicit rk45(rk45opts opts = {});


	traj integrate(const ContDynamicalSystem& sys, const State& x0, double t0, double t1) const;


private: 
	rk45opts opts;

};

