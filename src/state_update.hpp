#ifndef STATE_UPDATE_HPP
#define STATE_UPDATE_HPP

#include "point.hpp"
#include "cmath"

void func_delta(Point* globaldata, int numPoints, double cfl);
void state_update(Point* globaldata, int numPoints, Config configData, int iter, double res_old[1], int rk, double U[4], double Uold[4], double main_store[62]);
void state_update_wall(Point* globaldata, int idx, double max_res, double sig_res_sqr[1], double U[4], double Uold[4], int rk);
void state_update_outer(Point* globaldata, int idx, double Mach, double gamma, double pr_inf, double rho_inf, double theta, double max_res, double sig_res_sqr[1], double U[4], double Uold[4], int rk);
void state_update_interior(Point* globaldata, int idx, double max_res, double sig_res_sqr[1], double U[4], double Uold[4], int rk);


#endif