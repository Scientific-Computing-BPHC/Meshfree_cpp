#ifndef FLUX_RESIDUAL_HPP
#define FULX_RESIDUAL_HPP

void cal_flux_residual(Point* globaldata, int numPoints, Config configData, int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn, \
double* flux_res, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

void wallindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData, \
int* xpos_conn, int* xneg_conn, int* yneg_conn, double* flux_res, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

void outerindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData, \
int* xpos_conn, int* xneg_conn, int* ypos_conn, double* flux_res, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

void interiorindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData, \
int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn, double* flux_res, double* q, double* max_q, double* min_q, double* dq1, double* dq2);

#endif