#ifndef INTERIOR_FLUXES_HPP
#define INTERIOR_FLUXES_HPP

#include "point.hpp"
#include "split_fluxes_cuda.hpp"

__device__ void interior_dGx_pos(Point* globaldata, int idx, double Gxp[4], Config configData);

__device__ void interior_dGx_neg(Point* globaldata, int idx, double Gxn[4], Config configData);

__device__ void interior_dGy_pos(Point* globaldata, int idx, double Gyp[4], Config configData);

__device__ void interior_dGy_neg(Point* globaldata, int idx,  double Gyn[4], Config configData);

#endif