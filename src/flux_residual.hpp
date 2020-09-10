#ifndef FLUX_RESIDUAL_HPP
#define FULX_RESIDUAL_HPP

void cal_flux_residual(Point* globaldata, int numPoints, Config configData);

void wallindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData);

void outerindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData);

void interiorindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData);

#endif