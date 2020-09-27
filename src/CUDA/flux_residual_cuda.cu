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

	__shared__ double Gxp[4*32];
	__shared__ double Gxn[4*32];
	__shared__ double Gyp[4*32];
	__shared__ double Gyn[4*32];

	__shared__ double result[4*32];
	__shared__ double sig_del_x_del_f[4*32];
	__shared__ double sig_del_y_del_f[4*32];
	__shared__ double qtilde_i[4*32];
	__shared__ double qtilde_k[4*32];

	__shared__ double phi_i[4*32];
	__shared__ double phi_k[4*32];

	if(idx < numPoints)
	{

		if (globaldata[idx].flag_1 == 0)
			wallindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, phi_i, phi_k, configData);
		else if (globaldata[idx].flag_1 == 2)
			outerindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, phi_i, phi_k, configData);
		else if (globaldata[idx].flag_1 == 1)
			interiorindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, phi_i, phi_k, configData);

	}
}

__device__ void wallindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double* result, double* sig_del_x_del_f, double* sig_del_y_del_f, double* qtilde_i, double* qtilde_k, double* phi_i, double* phi_k, Config configData)
{

	wall_dGx_pos(globaldata, idx, Gxp, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, phi_i, phi_k, configData);
	wall_dGx_neg(globaldata, idx, Gxn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, phi_i, phi_k, configData);
	wall_dGy_neg(globaldata, idx, Gyn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, phi_i, phi_k, configData);

	//double Gtemp[4] = {0};

	for(int i=0; i<4; i++)
		globaldata[idx].flux_res[i] = globaldata[idx].delta * (Gxp[i + 4*threadIdx.x] + Gxn[i + 4*threadIdx.x] + Gyn[i+ 4*threadIdx.x]) * 2;

	// for(int i=0; i<4; i++)
	// {
	// 	globaldata[idx].flux_res[i] = Gtemp[i];
	// }
}

__device__ void outerindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double* result, double* sig_del_x_del_f, double* sig_del_y_del_f, double* qtilde_i, double* qtilde_k, double* phi_i, double* phi_k, Config configData)
{

	outer_dGx_pos(globaldata, idx, Gxp, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, phi_i, phi_k, configData);
	outer_dGx_neg(globaldata, idx, Gxn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, phi_i, phi_k, configData);
	outer_dGy_pos(globaldata, idx, Gyp, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, phi_i, phi_k, configData);

	//double Gtemp[4] = {0};

	for(int i=0; i<4; i++)
		globaldata[idx].flux_res[i] = globaldata[idx].delta * (Gxp[i + 4*threadIdx.x] + Gxn[i+ 4*threadIdx.x] + Gyp[i+ 4*threadIdx.x]);

	// for(int i=0; i<4; i++)
	// 	globaldata[idx].flux_res[i] = Gtemp[i];

}

__device__ void interiorindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double* result, double* sig_del_x_del_f, double* sig_del_y_del_f,  double* qtilde_i, double* qtilde_k, double* phi_i, double* phi_k, Config configData)
{
	interior_dGx_pos(globaldata, idx, Gxp, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, phi_i, phi_k, configData);
	interior_dGx_neg(globaldata, idx, Gxn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, phi_i, phi_k, configData);
	interior_dGy_pos(globaldata, idx, Gyp, result, sig_del_x_del_f, sig_del_y_del_f,  qtilde_i, qtilde_k, phi_i, phi_k, configData);
	interior_dGy_neg(globaldata, idx, Gyn, result, sig_del_x_del_f, sig_del_y_del_f,  qtilde_i, qtilde_k, phi_i, phi_k, configData); 

	// double Gtemp[4] = {0};

	for(int i=0; i<4; i++)
		globaldata[idx].flux_res[i] = globaldata[idx].delta * (Gxp[i + 4*threadIdx.x] + Gxn[i + 4*threadIdx.x] + Gyp[i+ 4*threadIdx.x] + Gyn[i+ 4*threadIdx.x]);

	// for(int i=0; i<4; i++)
	// 	globaldata[idx].flux_res[i] = Gtemp[i];
}
