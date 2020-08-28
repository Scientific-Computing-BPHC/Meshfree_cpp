#include "core.hpp"


inline double deg2rad(double radians) {
    return radians * (180.0 / M_PI);
}

double calculateTheta(Config configData)
{
	double theta = deg2rad(configData.core.aoa);
	return theta;
}

void getInitialPrimitive(Config configData, double primal[4])
{
	primal[0] = configData.core.rho_inf;
	double mach = configData.core.mach;
	double machcos = mach * cos(calculateTheta(configData));
	double machsin = mach * sin(calculateTheta(configData));
	primal[1] = machcos;
	primal[2] = machsin;
	primal[3] = configData.core.pr_inf;
}