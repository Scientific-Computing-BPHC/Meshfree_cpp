#include "point.hpp"
#include "flux_residual.hpp"
#include "wall_fluxes.hpp"

void cal_flux_residual(Point* globaldata, int numPoints, Config configData, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double main_store[62])
{
	double power = main_store[52];
	double limiter_flag = main_store[54];
	double vl_const = main_store[55];
	double gamma = main_store[58];

	for(int idx=0; idx<numPoints; idx++)
	{
		if (globaldata[idx].flag_1 == 0)
			wallindices_flux_residual(globaldata, gamma, idx, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const);
		else if (globaldata[idx].flag_1 == 2)
			outerindices_flux_residual(globaldata, gamma, idx, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const);
		else if (globaldata[idx].flag_1 == 1)
			interiorindices_flux_residual(globaldata, gamma, idx, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const);

	}
}

void wallindices_flux_residual(Point* globaldata, double gamma, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const)
{
	wall_dGx_pos(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gxp);
	//wall_dGx_neg(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gxn);
	//wall_dGy_neg(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gyn);

	for(int i=0; i<4; i++)
		Gxp[i] = (Gxp[i] + Gxn[i] + Gyn[i]) * 2;

	for(int i=0; i<4; i++)
		globaldata[idx].flux_res[i] = Gxp[i];
}

void outerindices_flux_residual(Point* globaldata, double gamma, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const)
{

}

void interiorindices_flux_residual(Point* globaldata, double gamma, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const)
{

}