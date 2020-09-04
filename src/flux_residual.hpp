#ifndef FLUX_RESIDUAL_HPP
#define FULX_RESIDUAL_HPP

void cal_flux_residual(Point* globaldata, int numPoints, Config configData, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double main_store[62]);

void wallindices_flux_residual(Point* globaldata, double gamma, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const);

void outerindices_flux_residual(Point* globaldata, double gamma, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const);

void interiorindices_flux_residual(Point* globaldata, double gamma, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const);

template <class Type>
bool isNan(Type var);


#endif