#ifndef OUTER_FLUXES_HPP
#define OUTER_FLUXES_HPP

#include "point.hpp"
#include "split_fluxes.hpp"

void outer_dGx_pos(Point* globaldata, int idx, double Gxp[4], Config configData, int* xpos_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

void outer_dGx_neg(Point* globaldata, int idx, double Gxn[4], Config configData, int* xneg_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

void outer_dGy_pos(Point* globaldata, int idx, double Gyp[4], Config configData, int* ypos_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

#endif