#ifndef WALL_FLUXES_HPP
#define WALL_FLUXES_HPP

#include "point.hpp"
#include "split_fluxes.hpp"

void wall_dGx_pos(Point* globaldata, int idx, double Gxp[4], Config configData, int* xpos_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

void wall_dGx_neg(Point* globaldata, int idx, double Gxn[4], Config configData, int* xneg_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

void wall_dGy_neg(Point* globaldata, int idx, double Gyn[4], Config configData, int* yneg_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2);
#endif