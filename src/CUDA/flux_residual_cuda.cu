#include "point.hpp"
#include "flux_residual_cuda.hpp"
#include "wall_fluxes_cuda.hpp"
#include "outer_fluxes_cuda.hpp"
#include "interior_fluxes_cuda.hpp"

template <class Type>
bool isNan(Type var)
{
    if(var!=var) return true;
    return false;
}

__global__ void cal_flux_residual_cuda(Point* globaldata, int numPoints, Config configData, dim3 thread_dim)
{
    int bx = blockIdx.x;
    int threadx = threadIdx.x;
    int idx = bx*thread_dim.x + threadx;

	double Gxp[4] = {0}, Gxn[4] = {0}, Gyp[4] = {0}, Gyn[4] = {0};

	if(idx < numPoints)
	{

		if (globaldata[idx].flag_1 == 0)
			wallindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, configData);
		else if (globaldata[idx].flag_1 == 2)
			outerindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, configData);
		else if (globaldata[idx].flag_1 == 1)
			interiorindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, configData);

	}
}

__device__ void wallindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData)
{


	wall_dGx_pos(globaldata, idx, Gxp, configData);
	wall_dGx_neg(globaldata, idx, Gxn, configData);
	wall_dGy_neg(globaldata, idx, Gyn, configData);

	double Gtemp[4] = {0};

	for(int i=0; i<4; i++)
	{
		Gtemp[i] = globaldata[idx].delta * (Gxp[i] + Gxn[i] + Gyn[i]) * 2;
	}

	for(int i=0; i<4; i++)
	{
		globaldata[idx].flux_res[i] = Gtemp[i];
	}
}

__device__ void outerindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData)
{

	outer_dGx_pos(globaldata, idx, Gxp, configData);
	outer_dGx_neg(globaldata, idx, Gxn, configData);
	outer_dGy_pos(globaldata, idx, Gyp, configData);

	double Gtemp[4] = {0};

	for(int i=0; i<4; i++)
		Gtemp[i] = globaldata[idx].delta * (Gxp[i] + Gxn[i] + Gyp[i]);

	for(int i=0; i<4; i++)
		globaldata[idx].flux_res[i] = Gtemp[i];

}

__device__ void interiorindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData)
{
	interior_dGx_pos(globaldata, idx, Gxp, configData);
	interior_dGx_neg(globaldata, idx, Gxn, configData);
	interior_dGy_pos(globaldata, idx, Gyp, configData);
	interior_dGy_neg(globaldata, idx, Gyn, configData); 

	double Gtemp[4] = {0};

	for(int i=0; i<4; i++)
		Gtemp[i] = globaldata[idx].delta * (Gxp[i] + Gxn[i] + Gyp[i] + Gyn[i]);

	for(int i=0; i<4; i++)
		globaldata[idx].flux_res[i] = Gtemp[i];
}
