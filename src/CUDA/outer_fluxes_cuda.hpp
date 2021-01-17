#ifndef OUTER_FLUXES_HPP
#define OUTER_FLUXES_HPP

/* Internal Headers */
#include "point.hpp"
#include "split_fluxes_cuda.hpp"

__device__ void outer_dGx_pos(Point* globaldata, int idx, double Gxp[4], double* result, double* sig_del_x_del_f, double* sig_del_y_del_f, \
double* qtilde_i, double* qtilde_k, Config configData, int* xpos_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

__device__ void outer_dGx_neg(Point* globaldata, int idx, double Gxn[4], double* result, double* sig_del_x_del_f, double* sig_del_y_del_f, \
double* qtilde_i, double* qtilde_k, Config configData, int* xneg_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

__device__ void outer_dGy_pos(Point* globaldata, int idx, double Gyp[4], double* result, double* sig_del_x_del_f, double* sig_del_y_del_f, \
double* qtilde_i, double* qtilde_k, Config configData, int* ypos_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

#endif