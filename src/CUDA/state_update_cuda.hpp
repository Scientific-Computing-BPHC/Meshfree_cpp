#ifndef STATE_UPDATE_HPP
#define STATE_UPDATE_HPP

#include "point.hpp"
#include "cmath"

__global__ void call_func_delta_cuda(Point* globaldata, int numPoints, double cfl, dim3 thread_dim);
__global__ void state_update_cuda(Point* globaldata, int numPoints, Config configData, int iter, double res_old[1], int rk, int rks, double* res_sqr, dim3 thread_dim);
__device__ void state_update_wall(Point* globaldata, int idx, double max_res, double* res_sqr, double U[4], double Uold[4], int rk, int euler);
__device__ void state_update_outer(Point* globaldata, int idx, double Mach, double gamma, double pr_inf, double rho_inf, double theta, double max_res, double* res_sqr, double U[4], double Uold[4], int rk, int euler);
__device__ void state_update_interior(Point* globaldata, int idx, double max_res, double* res_sqr, double U[4], double Uold[4], int rk, int euler);
__device__ void track_sig_res_sqr(double sig_res_sqr[1], int iter, int rk, int idx);

#endif