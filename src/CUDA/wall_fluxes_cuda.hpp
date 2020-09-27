#ifndef WALL_FLUXES_HPP
#define WALL_FLUXES_HPP

#include "point.hpp"
#include "split_fluxes_cuda.hpp"

__device__ void wall_dGx_pos(Point* globaldata, int idx, double Gxp[4], double* result, double* sig_del_x_del_f, double* sig_del_y_del_f,  double* qtilde_i, double* qtilde_k, double* phi_i, double* phi_k, Config configData);

__device__ void wall_dGx_neg(Point* globaldata, int idx, double Gxn[4], double* result, double* sig_del_x_del_f, double* sig_del_y_del_f,  double* qtilde_i, double* qtilde_k, double* phi_i, double* phi_k, Config configData);

__device__ void wall_dGy_neg(Point* globaldata, int idx, double Gyn[4], double* result, double* sig_del_x_del_f, double* sig_del_y_del_f,  double* qtilde_i, double* qtilde_k, double* phi_i, double* phi_k, Config configData);
#endif