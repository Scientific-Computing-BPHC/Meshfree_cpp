#ifndef INTERIOR_FLUXES_HPP
#define INTERIOR_FLUXES_HPP

#include "point.hpp"
#include "split_fluxes.hpp"

void interior_dGx_pos(Point* globaldata, int idx, double gamma, double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const, double Gxp[4]);

void interior_dGx_neg(Point* globaldata, int idx, double gamma, double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const, double Gxn[4]);

void interior_dGy_pos(Point* globaldata, int idx, double gamma, double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const, double Gyp[4]);

void interior_dGy_neg(Point* globaldata, int idx, double gamma, double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const, double Gyn[4]);

#endif