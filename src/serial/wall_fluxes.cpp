#include "wall_fluxes.hpp"
#include "limiters.hpp"
#include "quadrant_fluxes.hpp"
#include "split_fluxes.hpp"

template <class Type>
bool isNan(Type var)
{
    if(var!=var) return true;
    return false;
}

void wall_dGx_pos(Point* globaldata, int idx, double Gxp[4], Config configData)
{
	
    double power = configData.core.power;
    int limiter_flag = configData.core.limiter_flag;
    double vl_const = configData.core.vl_const;
    double gamma = configData.core.gamma;

    double phi_i[4] ={0}, phi_k[4] = {0}, G_i[4] = {0}, G_k[4] = {0}, result[4] = {0}, qtilde_i[4] = {0}, qtilde_k[4] = {0}, sig_del_x_del_f[4] ={0}, sig_del_y_del_f[4] = {0};

    double sig_del_x_sqr = 0.0;
	double sig_del_y_sqr = 0.0;
	double sig_del_x_del_y = 0.0;

	for(int i=0; i<4; i++)
	{
		sig_del_x_del_f[i] = 0.0;
		sig_del_y_del_f[i] = 0.0;
	}

	double x_i = globaldata[idx].x;
	double y_i = globaldata[idx].y;

	double nx = globaldata[idx].nx;
	double ny = globaldata[idx].ny;

	double tx = ny;
	double ty = -nx;

	for(int i=0; i<20; i++)
	{
		int conn = globaldata[idx].xpos_conn[i];
		if(conn == 0) break;

		conn = conn - 1;

		double delta_x, delta_y, delta_s_weights, delta_n_weights;
		std::tie(delta_x, delta_y, delta_s_weights, delta_n_weights, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y) = connectivity_stats(x_i, y_i, nx, ny, power, globaldata[conn].x, globaldata[conn].y, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y);

		calculate_qtile(qtilde_i, qtilde_k, globaldata, idx, conn, delta_x, delta_y, vl_const, gamma, limiter_flag, phi_i, phi_k);

		qtilde_to_primitive(result, qtilde_i, gamma);

		flux_quad_GxII(G_i, nx, ny, result[0], result[1], result[2], result[3]);

		qtilde_to_primitive(result, qtilde_k, gamma);
        flux_quad_GxII(G_k, nx, ny, result[0], result[1], result[2], result[3]);

        update_delf(sig_del_x_del_f, sig_del_y_del_f, G_k, G_i, delta_s_weights, delta_n_weights);
    }
	
	cout<<endl<<"Qtilde_i,k, Idx: "<<idx<<endl;
	for(int index = 0; index<4; index++)
	{
		cout<<std::fixed<<std::setprecision(17)<<qtilde_i[index]<<"   ";
	}
	for(int index = 0; index<4; index++)
	{
		cout<<std::fixed<<std::setprecision(17)<<qtilde_k[index]<<"   ";
	}


    double det = sig_del_x_sqr * sig_del_y_sqr - sig_del_x_del_y * sig_del_x_del_y;
    double one_by_det = 1.0/det;
    for(int iter =0; iter<4; iter++)
    {
    	Gxp[iter] = (sig_del_x_del_f[iter]*sig_del_y_sqr - sig_del_y_del_f[iter]*sig_del_x_del_y)*one_by_det;

    }
	
}

void wall_dGx_neg(Point* globaldata, int idx, double Gxn[4], Config configData)
{
    double power = configData.core.power;
    int limiter_flag = configData.core.limiter_flag;
    double vl_const = configData.core.vl_const;
    double gamma = configData.core.gamma;

    double phi_i[4] ={0}, phi_k[4] = {0}, G_i[4] = {0}, G_k[4] = {0}, result[4] = {0}, qtilde_i[4] = {0}, qtilde_k[4] = {0}, sig_del_x_del_f[4] ={0}, sig_del_y_del_f[4] = {0};

    double sig_del_x_sqr = 0.0;
	double sig_del_y_sqr = 0.0;
	double sig_del_x_del_y = 0.0;

	for(int i=0; i<4; i++)
	{
		sig_del_x_del_f[i] = 0.0;
		sig_del_y_del_f[i] = 0.0;
	}

	double x_i = globaldata[idx].x;
	double y_i = globaldata[idx].y;

	double nx = globaldata[idx].nx;
	double ny = globaldata[idx].ny;

	double tx = ny;
	double ty = -nx;

	for(int i=0; i<20; i++)
	{
		int conn = globaldata[idx].xneg_conn[i];
		if(conn == 0) break;

		conn = conn - 1;

		double delta_x, delta_y, delta_s_weights, delta_n_weights;
		std::tie(delta_x, delta_y, delta_s_weights, delta_n_weights, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y) = connectivity_stats(x_i, y_i, nx, ny, power, globaldata[conn].x, globaldata[conn].y, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y);

		calculate_qtile(qtilde_i, qtilde_k, globaldata, idx, conn, delta_x, delta_y, vl_const, gamma, limiter_flag, phi_i, phi_k);

		qtilde_to_primitive(result, qtilde_i, gamma);

		flux_quad_GxI(G_i, nx, ny, result[0], result[1], result[2], result[3]);

		qtilde_to_primitive(result, qtilde_k, gamma);
        flux_quad_GxI(G_k, nx, ny, result[0], result[1], result[2], result[3]);

        update_delf(sig_del_x_del_f, sig_del_y_del_f, G_k, G_i, delta_s_weights, delta_n_weights);

     }

    double det = sig_del_x_sqr * sig_del_y_sqr - sig_del_x_del_y * sig_del_x_del_y;
    double one_by_det = 1.0/det;
    for(int iter =0; iter<4; iter++)
    {
    	Gxn[iter] = (sig_del_x_del_f[iter]*sig_del_y_sqr - sig_del_y_del_f[iter]*sig_del_x_del_y)*one_by_det;
    }

	
}

void wall_dGy_neg(Point* globaldata, int idx, double Gyn[4], Config configData)
{
	double power = configData.core.power;
    int limiter_flag = configData.core.limiter_flag;
    double vl_const = configData.core.vl_const;
    double gamma = configData.core.gamma;

    double phi_i[4] ={0}, phi_k[4] = {0}, G_i[4] = {0}, G_k[4] = {0}, result[4] = {0}, qtilde_i[4] = {0}, qtilde_k[4] = {0}, sig_del_x_del_f[4] ={0}, sig_del_y_del_f[4] = {0};

    double sig_del_x_sqr = 0.0;
	double sig_del_y_sqr = 0.0;
	double sig_del_x_del_y = 0.0;

	for(int i=0; i<4; i++)
	{
		sig_del_x_del_f[i] = 0.0;
		sig_del_y_del_f[i] = 0.0;
	}

	double x_i = globaldata[idx].x;
	double y_i = globaldata[idx].y;

	double nx = globaldata[idx].nx;
	double ny = globaldata[idx].ny;

	double tx = ny;
	double ty = -nx;

	for(int i=0; i<20; i++)
	{
		int conn = globaldata[idx].yneg_conn[i];
		if(conn == 0) break;

		conn = conn - 1;

		double delta_x, delta_y, delta_s_weights, delta_n_weights;
		std::tie(delta_x, delta_y, delta_s_weights, delta_n_weights, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y) = connectivity_stats(x_i, y_i, nx, ny, power, globaldata[conn].x, globaldata[conn].y, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y);

		calculate_qtile(qtilde_i, qtilde_k, globaldata, idx, conn, delta_x, delta_y, vl_const, gamma, limiter_flag, phi_i, phi_k);

		qtilde_to_primitive(result, qtilde_i, gamma);

		flux_Gyn(G_i, nx, ny, result[0], result[1], result[2], result[3]);

		qtilde_to_primitive(result, qtilde_k, gamma);

        flux_Gyn(G_k, nx, ny, result[0], result[1], result[2], result[3]);

        update_delf(sig_del_x_del_f, sig_del_y_del_f, G_k, G_i, delta_s_weights, delta_n_weights);
     }

    double det = sig_del_x_sqr * sig_del_y_sqr - sig_del_x_del_y * sig_del_x_del_y;
    double one_by_det = 1.0/det;
    for(int iter =0; iter<4; iter++)
    {
    	Gyn[iter] = (sig_del_y_del_f[iter]*sig_del_x_sqr - sig_del_x_del_f[iter]*sig_del_x_del_y)*one_by_det;
    }

	
}