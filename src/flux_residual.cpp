#include "point.hpp"
#include "flux_residual.hpp"
#include "wall_fluxes.hpp"
#include "outer_fluxes.hpp"
#include "interior_fluxes.hpp"

template <class Type>
bool isNan(Type var)
{
    if(var!=var) return true;
    return false;
}

void cal_flux_residual(Point* globaldata, int numPoints, Config configData, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double main_store[62])
{
	double power = main_store[52];
	double limiter_flag = main_store[54];
	double vl_const = main_store[55];
	double gamma = main_store[58];

	for(int idx=0; idx<numPoints; idx++)
	{

		// if (idx == 668)
		// {
		// 	cout<<"\n PHI_K ISSS: \n";
		// 	for(int k=0; k<4; k++) cout<<phi_k[k] <<"  ";
		// 	cout<<"\n";	
		// }

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
	
	// if (idx == 0)
	// {
	// 	cout<<"\n Inside outerr PHI_K ISSS: \n";
	// 	for(int k=0; k<4; k++) cout<<phi_k[k] <<"  ";
	// 	cout<<"\n";	
	// }

	// if(idx == 0)
	// {
	// 	cout<<"Lesse how other variables are doing, and then exit at idx = 0, iter = 0 just intoooo"<<endl;
	// 	for(int l=0; l<4; l++)
	// 	{
	// 		cout<<"Gamma: "<<gamma<<"Phi:  "<<phi_i[l]<<" and k "<<phi_k[l]<<" G_i:  "<<G_i[l]<<" G_k: "<<G_k[l]<<" result:  "<<result[l]<<" qtilde_i: "<<qtilde_i[l]<<" qtilde_k:  "<<qtilde_k[l]<<" sigx: "<<sig_del_x_del_f[l]<<" sigy: "<<sig_del_y_del_f[l]<<"  "<<power<<"  "<<limiter_flag<<"  "<<vl_const<<endl;
	// 	}
	// }

	wall_dGx_pos(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gxp);
	wall_dGx_neg(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gxn);
	wall_dGy_neg(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gyn);

	// if(idx == 0)
	// {
	// 	cout<<"Lesse how other variables are doing, and then exit at idx = 0, iter = 0"<<endl;
	// 	for(int l=0; l<4; l++)
	// 	{
	// 		cout<<"Gamma: "<<gamma<<"Phi:  "<<phi_i[l]<<" and k "<<phi_k[l]<<" G_i:  "<<G_i[l]<<" G_k: "<<G_k[l]<<" result:  "<<result[l]<<" qtilde_i: "<<qtilde_i[l]<<" qtilde_k:  "<<qtilde_k[l]<<" sigx: "<<sig_del_x_del_f[l]<<" sigy: "<<sig_del_y_del_f[l]<<"  "<<power<<"  "<<limiter_flag<<"  "<<vl_const<<endl;
	// 	}
	// }

	// if(idx == 668)
	// {
	// 	cout<<"Lesse how G is"<<endl;
	// 	for(int l=0; l<4; l++)
	// 	{
	// 		cout<<Gxp[l]<<"  "<<Gxn[l]<<"  "<<Gyp[l]<<"  "<<endl;
	// 	}
	// 	//exit(0);
	// }

	debug_Gs_again(idx, Gxp, Gxn, Gyp, Gyn);

	for(int i=0; i<4; i++)
	{
		Gxp[i] = (Gxp[i] + Gxn[i] + Gyn[i]) * 2;
		if(isNan(Gxp[i]))
		{
			cout<<"Nan in wall flux Gxp: "<<endl;
			cout<<"Idx is: "<<idx<<endl;
			cout<<"i is: "<<i<<endl;
			cout<<"Gxp"<<Gxp[i]<<endl;
			cout<<"Gxn"<<Gxn[i]<<endl;
			cout<<"Gyn"<<Gyn[i]<<endl;
			exit(0);
		}
	}

	for(int i=0; i<4; i++)
	{
		globaldata[idx].flux_res[i] = Gxp[i];
		if(isNan(globaldata[idx].flux_res[i]))
		{
			cout<<"Nan in wall flux: "<<endl;
			cout<<"Idx is: "<<idx<<endl;
			cout<<"i is: "<<i<<endl;
			exit(0);
		}
	}
}

void outerindices_flux_residual(Point* globaldata, double gamma, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const)
{

	// if (idx == 1356)
	// {
	// 	cout<<"\n Inside outerr PHI_K ISSS: \n";
	// 	for(int k=0; k<4; k++) cout<<phi_k[k] <<"  ";
	// 	cout<<"\n";	
	// }

	// if(idx == 1356)
	// {
	// 	cout<<"Lesse how other variables are doing"<<endl;
	// 	for(int l=0; l<4; l++)
	// 	{
	// 		cout<<"Gamma: "<<gamma<<"Phi:  "<<phi_i[l]<<" and k "<<phi_k[l]<<"  "<<G_i[l]<<"  "<<G_k[l]<<"  "<<qtilde_i[l]<<"  "<<qtilde_k[l]<<"  "<<endl;
	// 	}
	// }

	outer_dGx_pos(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gxp);
	outer_dGx_neg(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gxn);
	outer_dGy_pos(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gyp);

	// if(idx == 1356)
	// {
	// 	cout<<"Lesse how G is"<<endl;
	// 	for(int l=0; l<4; l++)
	// 	{
	// 		cout<<Gxp[l]<<"  "<<Gxn[l]<<"  "<<Gyp[l]<<"  "<<endl;
	// 	}
	// 	//exit(0);
	// }

	debug_Gs_again(idx, Gxp, Gxn, Gyp, Gyn);

	for(int i=0; i<4; i++)
		Gxp[i] = (Gxp[i] + Gxn[i] + Gyp[i]);

	for(int i=0; i<4; i++)
		globaldata[idx].flux_res[i] = Gxp[i];

}

void interiorindices_flux_residual(Point* globaldata, double gamma, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const)
{
	interior_dGx_pos(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gxp);
	interior_dGx_neg(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gxn);
	interior_dGy_pos(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gyp);
	interior_dGy_neg(globaldata, idx, gamma, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, power, limiter_flag, vl_const, Gyn); 

	debug_Gs_again(idx, Gxp, Gxn, Gyp, Gyn);


	for(int i=0; i<4; i++)
		Gxp[i] = (Gxp[i] + Gxn[i] + Gyp[i] + Gyn[i]);

	for(int i=0; i<4; i++)
		globaldata[idx].flux_res[i] = Gxp[i];
}

void debug_Gs_again(int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4])
{
    std::ofstream fdebug("debug_Gs_again.txt", std::ios_base::app);

    fdebug<<"\nChecking the Outputs: \n";

    fdebug<<"Index (idx+1): "<<idx+1<<"\n";
    fdebug<<"\nGxp: ";
    for(int i=0; i<4; i++)
        fdebug<<Gxp[i]<<"  ";
    fdebug<<"\nGxn: ";
    for(int i=0; i<4; i++)
        fdebug<<Gxn[i]<<"  ";
    fdebug<<"\nGyp: ";
    for(int i=0; i<4; i++)
        fdebug<<Gyp[i]<<"  ";
    fdebug<<"\nGyn: ";
    for(int i=0; i<4; i++)
        fdebug<<Gyn[i]<<"  ";
    fdebug.close();
}