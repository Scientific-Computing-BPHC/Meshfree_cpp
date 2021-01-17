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

__global__ void cal_flux_residual_cuda(Point* globaldata, int numPoints, Config configData, dim3 thread_dim, int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn, double* flux_res, \
	double* q, double* max_q, double* min_q, double* dq1, double* dq2)
{
    int bx = blockIdx.x;
    int threadx = threadIdx.x;
    int idx = bx*thread_dim.x + threadx;

	__shared__ double Gxp[4*64];
	__shared__ double Gxn[4*64];
	__shared__ double Gyp[4*64];
	__shared__ double Gyn[4*64];

	__shared__ double result[4*64];
	__shared__ double sig_del_x_del_f[4*64];
	__shared__ double sig_del_y_del_f[4*64];
	__shared__ double qtilde_i[4*64];
	__shared__ double qtilde_k[4*64];

	if(idx < numPoints)
	{

		if (globaldata[idx].flag_1 == 0)
			wallindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData,\
				xpos_conn, xneg_conn, yneg_conn, flux_res, q, max_q, min_q, dq1, dq2);
			
		else if (globaldata[idx].flag_1 == 2)
			outerindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, \
			xpos_conn, xneg_conn, ypos_conn, flux_res, q, max_q, min_q, dq1, dq2);
			
		else if (globaldata[idx].flag_1 == 1)
			interiorindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, \
				xpos_conn, xneg_conn, ypos_conn, yneg_conn, flux_res, q, max_q, min_q, dq1, dq2);

	}
}

__device__ void wallindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], \
	double* result, double* sig_del_x_del_f, double* sig_del_y_del_f,  double* qtilde_i, double* qtilde_k, Config configData, \
	int* xpos_conn, int* xneg_conn, int* yneg_conn, double* flux_res, double* q, double* max_q, double* min_q, double* dq1, double* dq2)
{

	wall_dGx_pos(globaldata, idx, Gxp, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, xpos_conn, q, max_q, min_q, dq1, dq2);
	wall_dGx_neg(globaldata, idx, Gxn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, xneg_conn, q, max_q, min_q, dq1, dq2);
	wall_dGy_neg(globaldata, idx, Gyn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, yneg_conn, q, max_q, min_q, dq1, dq2);

	//double Gtemp[4] = {0};

	for(int i=0; i<4; i++)
		flux_res[idx*4 + i] = globaldata[idx].delta * (Gxp[i + 4*threadIdx.x] + Gxn[i + 4*threadIdx.x] + Gyn[i+ 4*threadIdx.x]) * 2;

	// if(idx ==0)
	// {
	// 	printf("\n");
	// 	for(int index = 0; index<4; index++)
	// 	{
	// 		printf("%.17f   ", flux_res[idx*4 + index]);
	// 	}
	// }
}

__device__ void outerindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], \
	double* result, double* sig_del_x_del_f, double* sig_del_y_del_f, double* qtilde_i, double* qtilde_k, Config configData, \
	int* xpos_conn, int* xneg_conn, int* ypos_conn, double* flux_res, double* q, double* max_q, double* min_q, double* dq1, double* dq2)
{

	outer_dGx_pos(globaldata, idx, Gxp, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, xpos_conn, q, max_q, min_q, dq1, dq2);
	outer_dGx_neg(globaldata, idx, Gxn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, xneg_conn, q, max_q, min_q, dq1, dq2);
	outer_dGy_pos(globaldata, idx, Gyp, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, ypos_conn, q, max_q, min_q, dq1, dq2);

	//double Gtemp[4] = {0};

	for(int i=0; i<4; i++)
		flux_res[idx*4 + i] = globaldata[idx].delta * (Gxp[i + 4*threadIdx.x] + Gxn[i+ 4*threadIdx.x] + Gyp[i+ 4*threadIdx.x]);

	// for(int i=0; i<4; i++)
	// 	globaldata[idx].flux_res[i] = Gtemp[i];

}

__device__ void interiorindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], \
	double* result, double* sig_del_x_del_f, double* sig_del_y_del_f,  double* qtilde_i, double* qtilde_k, Config configData, \
	int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn, double* flux_res, double* q, double* max_q, double* min_q, double* dq1, double* dq2)
{
	interior_dGx_pos(globaldata, idx, Gxp, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, xpos_conn, q, max_q, min_q, dq1, dq2);
	interior_dGx_neg(globaldata, idx, Gxn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, xneg_conn, q, max_q, min_q, dq1, dq2);
	interior_dGy_pos(globaldata, idx, Gyp, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, ypos_conn, q, max_q, min_q, dq1, dq2);
	interior_dGy_neg(globaldata, idx, Gyn, result, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k, configData, yneg_conn, q, max_q, min_q, dq1, dq2); 

	// double Gtemp[4] = {0};

	for(int i=0; i<4; i++)
		flux_res[idx*4 + i] = globaldata[idx].delta * (Gxp[i + 4*threadIdx.x] + Gxn[i + 4*threadIdx.x] + Gyp[i+ 4*threadIdx.x] + Gyn[i+ 4*threadIdx.x]);

	// for(int i=0; i<4; i++)
	// 	globaldata[idx].flux_res[i] = Gtemp[i];
}
