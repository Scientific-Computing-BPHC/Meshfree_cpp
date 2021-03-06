#include "outer_fluxes_cuda.hpp"
#include "limiters_cuda.hpp"
#include "quadrant_fluxes_cuda.hpp"
#include "split_fluxes_cuda.hpp"

__device__ void outer_dGx_pos(Point* globaldata, int idx, double Gxp[4], double* result,  double* sig_del_x_del_f, double* sig_del_y_del_f, \
	double* qtilde_i, double* qtilde_k, Config configData, int* xpos_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2)
{
	double power = configData.core.power;
    int limiter_flag = configData.core.limiter_flag;
    double vl_const = configData.core.vl_const;
    double gamma = configData.core.gamma;
	double phi_i[4] ={0}, phi_k[4] = {0}, G_i[4] = {0}, G_k[4] = {0};
	
	double sig_del_x_sqr = 0.0;
	double sig_del_y_sqr = 0.0;
	double sig_del_x_del_y = 0.0;

	for(int i=0; i<4; i++)
	{
		sig_del_x_del_f[i + 4*threadIdx.x] = 0.0;
		sig_del_y_del_f[i + 4*threadIdx.x] = 0.0;
	}

	double x_i = globaldata[idx].x;
	double y_i = globaldata[idx].y;

	double nx = globaldata[idx].nx;
	double ny = globaldata[idx].ny;

	double tx = ny;
	double ty = -nx;

	for(int i=0; i<30; i++)
	{
		int conn = xpos_conn[idx*30 + i];
		if(conn == -1) break;

		conn = conn - 1;

		double delta_x, delta_y, delta_s_weights, delta_n_weights;
		connectivity_stats(x_i, y_i, nx, ny, power, globaldata[conn].x, globaldata[conn].y, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y, delta_x, delta_y, delta_s_weights, delta_n_weights);
		
		calculate_qtile(qtilde_i, qtilde_k, globaldata, idx, conn, delta_x, delta_y, vl_const, gamma, limiter_flag, phi_i, phi_k, q, max_q, min_q, dq1, dq2);

		qtilde_to_primitive(result, qtilde_i, gamma);

		flux_quad_GxIII(G_i, nx, ny, result[0 + 4*threadIdx.x], result[1 + 4*threadIdx.x], result[2 + 4*threadIdx.x], result[3 + 4*threadIdx.x]);

		qtilde_to_primitive(result, qtilde_k, gamma);
        flux_quad_GxIII(G_k, nx, ny, result[0 + 4*threadIdx.x], result[1 + 4*threadIdx.x], result[2 + 4*threadIdx.x], result[3 + 4*threadIdx.x]);

        update_delf(sig_del_x_del_f, sig_del_y_del_f, G_k, G_i, delta_s_weights, delta_n_weights);
     }

    double det = sig_del_x_sqr * sig_del_y_sqr - sig_del_x_del_y * sig_del_x_del_y;
    double one_by_det = 1.0/det;
    for(int iter =0; iter<4; iter++)
    	Gxp[iter + 4*threadIdx.x] = (sig_del_x_del_f[iter + 4*threadIdx.x]*sig_del_y_sqr - sig_del_y_del_f[iter + 4*threadIdx.x]*sig_del_x_del_y)*one_by_det;

	
}

__device__ void outer_dGx_neg(Point* globaldata, int idx, double Gxn[4], double* result,  double* sig_del_x_del_f, double* sig_del_y_del_f, \
	double* qtilde_i, double* qtilde_k, Config configData, int* xneg_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2)
{
	double sig_del_x_sqr = 0.0;
	double sig_del_y_sqr = 0.0;
	double sig_del_x_del_y = 0.0;

	double power = configData.core.power;
    int limiter_flag = configData.core.limiter_flag;
    double vl_const = configData.core.vl_const;
    double gamma = configData.core.gamma;

    double phi_i[4] ={0}, phi_k[4] = {0}, G_i[4] = {0}, G_k[4] = {0};

	for(int i=0; i<4; i++)
	{
		sig_del_x_del_f[i + 4*threadIdx.x] = 0.0;
		sig_del_y_del_f[i + 4*threadIdx.x] = 0.0;
	}

	double x_i = globaldata[idx].x;
	double y_i = globaldata[idx].y;

	double nx = globaldata[idx].nx;
	double ny = globaldata[idx].ny;

	double tx = ny;
	double ty = -nx;

	for(int i=0; i<30; i++)
	{
		int conn = xneg_conn[idx*30 + i];
		if(conn == -1) break;

		conn = conn - 1;

		double delta_x, delta_y, delta_s_weights, delta_n_weights;
		connectivity_stats(x_i, y_i, nx, ny, power, globaldata[conn].x, globaldata[conn].y, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y, delta_x, delta_y, delta_s_weights, delta_n_weights);

		calculate_qtile(qtilde_i, qtilde_k, globaldata, idx, conn, delta_x, delta_y, vl_const, gamma, limiter_flag, phi_i, phi_k, q, max_q, min_q, dq1, dq2);

		qtilde_to_primitive(result, qtilde_i, gamma);

		flux_quad_GxIV(G_i, nx, ny, result[0 + 4*threadIdx.x], result[1 + 4*threadIdx.x], result[2 + 4*threadIdx.x], result[3 + 4*threadIdx.x]);

		qtilde_to_primitive(result, qtilde_k, gamma);
        flux_quad_GxIV(G_k, nx, ny, result[0 + 4*threadIdx.x], result[1 + 4*threadIdx.x], result[2 + 4*threadIdx.x], result[3 + 4*threadIdx.x]);

        update_delf(sig_del_x_del_f, sig_del_y_del_f, G_k, G_i, delta_s_weights, delta_n_weights);
     }

    double det = sig_del_x_sqr * sig_del_y_sqr - sig_del_x_del_y * sig_del_x_del_y;
    double one_by_det = 1.0/det;
    for(int iter =0; iter<4; iter++)
    	Gxn[iter+ 4*threadIdx.x] = (sig_del_x_del_f[iter + 4*threadIdx.x]*sig_del_y_sqr - sig_del_y_del_f[iter + 4*threadIdx.x]*sig_del_x_del_y)*one_by_det;


}

__device__ void outer_dGy_pos(Point* globaldata, int idx, double Gyp[4], double* result,  double* sig_del_x_del_f, double* sig_del_y_del_f, \
	double* qtilde_i, double* qtilde_k, Config configData, int* ypos_conn, double* q, double* max_q, double* min_q, double* dq1, double* dq2)
{
	double sig_del_x_sqr = 0.0;
	double sig_del_y_sqr = 0.0;
	double sig_del_x_del_y = 0.0;

	double power = configData.core.power;
    int limiter_flag = configData.core.limiter_flag;
    double vl_const = configData.core.vl_const;
    double gamma = configData.core.gamma;

	double phi_i[4] ={0}, phi_k[4] = {0}, G_i[4] = {0}, G_k[4] = {0};
	
	for(int i=0; i<4; i++)
	{
		sig_del_x_del_f[i + 4*threadIdx.x] = 0.0;
		sig_del_y_del_f[i + 4*threadIdx.x] = 0.0;
	}

	double x_i = globaldata[idx].x;
	double y_i = globaldata[idx].y;

	double nx = globaldata[idx].nx;
	double ny = globaldata[idx].ny;

	double tx = ny;
	double ty = -nx;

	for(int i=0; i<30; i++)
	{
		int conn = ypos_conn[idx*30 + i];
		if(conn == -1) break;

		conn = conn - 1;

		double delta_x, delta_y, delta_s_weights, delta_n_weights;
		connectivity_stats(x_i, y_i, nx, ny, power, globaldata[conn].x, globaldata[conn].y, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y, delta_x, delta_y, delta_s_weights, delta_n_weights);
		
		calculate_qtile(qtilde_i, qtilde_k, globaldata, idx, conn, delta_x, delta_y, vl_const, gamma, limiter_flag, phi_i, phi_k, q, max_q, min_q, dq1, dq2);

		qtilde_to_primitive(result, qtilde_i, gamma);

		flux_Gyp(G_i, nx, ny, result[0 + 4*threadIdx.x], result[1 + 4*threadIdx.x], result[2 + 4*threadIdx.x], result[3 + 4*threadIdx.x]);

		qtilde_to_primitive(result, qtilde_k, gamma);

        flux_Gyp(G_k, nx, ny, result[0 + 4*threadIdx.x], result[1 + 4*threadIdx.x], result[2 + 4*threadIdx.x], result[3 + 4*threadIdx.x]);

        update_delf(sig_del_x_del_f, sig_del_y_del_f, G_k, G_i, delta_s_weights, delta_n_weights);
     }

    double det = sig_del_x_sqr * sig_del_y_sqr - sig_del_x_del_y * sig_del_x_del_y;
    double one_by_det = 1.0/det;
    for(int iter =0; iter<4; iter++)
    	Gyp[iter+ 4*threadIdx.x] = (sig_del_y_del_f[iter + 4*threadIdx.x]*sig_del_x_sqr - sig_del_x_del_f[iter + 4*threadIdx.x]*sig_del_x_del_y)*one_by_det;

	
}