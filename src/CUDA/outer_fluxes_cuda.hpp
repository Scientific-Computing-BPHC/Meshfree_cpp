#ifndef OUTER_FLUXES_HPP
#define OUTER_FLUXES_HPP

#include "point.hpp"
#include "split_fluxes_cuda.hpp"

__device__ void outer_dGx_pos(Point* globaldata, int idx, double Gxp[4], double* result, Config configData);

__device__ void outer_dGx_neg(Point* globaldata, int idx, double Gxn[4], double* result, Config configData);

__device__ void outer_dGy_pos(Point* globaldata, int idx, double Gyp[4], double* result, Config configData);

#endif