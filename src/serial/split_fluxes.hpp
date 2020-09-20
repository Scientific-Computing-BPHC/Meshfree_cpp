#ifndef SPLIT_FLUXES_HPP
#define SPLIT_FLUXES_HPP

#include "point.hpp"

void flux_Gyn(double Gyn[4], double nx, double ny, double u1, double u2, double rho, double pr);

void flux_Gyp(double Gyp[4], double nx, double ny, double u1, double u2, double rho, double pr);

void flux_Gxp(double Gxp[4], double nx, double ny, double u1, double u2, double rho, double pr);

void flux_Gxn(double Gxn[4], double nx, double ny, double u1, double u2, double rho, double pr);

void flux_Gx(double Gx[4], double nx, double ny, double u1, double u2, double rho, double pr);

void flux_Gy(double Gy[4], double nx, double ny, double u1, double u2, double rho, double pr);
#endif