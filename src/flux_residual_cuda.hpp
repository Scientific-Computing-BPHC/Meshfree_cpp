#ifndef FLUX_RESIDUAL_HPP
#define FULX_RESIDUAL_HPP

__global__ void cal_flux_residual_cuda(Point* globaldata, int numPoints, Config configData, dim3 thread_dim);

__device__ void wallindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData);

__device__ void outerindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData);

__device__ void interiorindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData);

#endif