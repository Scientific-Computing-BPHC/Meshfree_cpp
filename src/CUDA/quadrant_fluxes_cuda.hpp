#ifndef QUADRANT_FLUXES_HPP
#define QUADRANT_FLUXES_HPP

#include "point.hpp"
#include "cmath"

__device__ void flux_quad_GxII(double G[4], double nx, double ny, double u1, double u2, double rho, double pr);

__device__ void flux_quad_GxI(double G[4], double nx, double ny, double u1, double u2, double rho, double pr);

__device__ void flux_quad_GxIII(double G[4], double nx, double ny, double u1, double u2, double rho, double pr);

__device__ void flux_quad_GxIV(double G[4], double nx, double ny, double u1, double u2, double rho, double pr);

#endif