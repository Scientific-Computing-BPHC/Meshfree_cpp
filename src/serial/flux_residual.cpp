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

void cal_flux_residual(Point* globaldata, int numPoints, Config configData, int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn, double* flux_res, \
	double* q, double* max_q, double* min_q, double* dq1, double* dq2)
{

	double Gxp[4] = {0}, Gxn[4] = {0}, Gyp[4] = {0}, Gyn[4] = {0};
	for(int idx=0; idx<numPoints; idx++)
	{
		if (globaldata[idx].flag_1 == 0)
			wallindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, configData, \
				xpos_conn, xneg_conn, yneg_conn, flux_res, q, max_q, min_q, dq1, dq2);

		else if (globaldata[idx].flag_1 == 2)
			outerindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, configData, \
				xpos_conn, xneg_conn, ypos_conn, flux_res, q, max_q, min_q, dq1, dq2);

		else if (globaldata[idx].flag_1 == 1)
			interiorindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, configData, \
				xpos_conn, xneg_conn, ypos_conn, yneg_conn, flux_res, q, max_q, min_q, dq1, dq2);

	}
}

void wallindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData, \
	int* xpos_conn, int* xneg_conn, int* yneg_conn, double* flux_res, double* q, double* max_q, double* min_q, double* dq1, double* dq2)
{
	cout<<endl<<"Indic: "<<idx<<endl;
	cout<<endl<<"Delta: "<<globaldata[idx].delta<<endl;
	cout<<endl<<"Printing dq1, dq2: "<<endl;

	for(int index = 0; index<4; index++)
	{
		cout<<std::fixed<<std::setprecision(17)<<dq1[idx*4 + index]<<"   ";
	}

	for(int index = 0; index<4; index++)
	{
		cout<<std::fixed<<std::setprecision(17)<<dq2[idx*4 + index]<<"   ";
	}


	wall_dGx_pos(globaldata, idx, Gxp,  configData, xpos_conn, q, max_q, min_q, dq1, dq2);
	wall_dGx_neg(globaldata, idx, Gxn,  configData, xneg_conn, q, max_q, min_q, dq1, dq2);
	wall_dGy_neg(globaldata, idx, Gyn,  configData, yneg_conn, q, max_q, min_q, dq1, dq2);


	for(int i=0; i<4; i++)
		flux_res[idx*4 + i] = globaldata[idx].delta * (Gxp[i] + Gxn[i] + Gyn[i]) * 2;
	
	cout<<endl<<"Printing Gs: "<<endl;
	for(int index = 0; index<4; index++)
	{
		cout<<std::fixed<<std::setprecision(17)<<Gxp[index]<<"   ";
	}
	cout<<endl;
	for(int index = 0; index<4; index++)
	{
		cout<<std::fixed<<std::setprecision(17)<<Gxn[index]<<"   ";
	}
	cout<<endl;
	for(int index = 0; index<4; index++)
	{
		cout<<std::fixed<<std::setprecision(17)<<Gyn[index]<<"   ";
	}

}

void outerindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData, \
	int* xpos_conn, int* xneg_conn, int* ypos_conn, double* flux_res, double* q, double* max_q, double* min_q, double* dq1, double* dq2)
{

	outer_dGx_pos(globaldata, idx, Gxp, configData, xpos_conn, q, max_q, min_q, dq1, dq2);
	outer_dGx_neg(globaldata, idx, Gxn, configData, xneg_conn, q, max_q, min_q, dq1, dq2);
	outer_dGy_pos(globaldata, idx, Gyp, configData, ypos_conn, q, max_q, min_q, dq1, dq2);


	if(idx == 1356)
	{
	    cout<<endl<<"For 1356, printing Gs: "<<endl;
        for(int index = 0; index<4; index++)
        {
            cout<<std::fixed<<std::setprecision(17)<<Gxp[index]<<"   ";
        }
        cout<<endl;
        for(int index = 0; index<4; index++)
        {
            cout<<std::fixed<<std::setprecision(17)<<Gxn[index]<<"   ";
        }
        cout<<endl;
        for(int index = 0; index<4; index++)
        {
            cout<<std::fixed<<std::setprecision(17)<<Gyp[index]<<"   ";
        }
        // cout<<endl;
        // for(int index = 0; index<4; index++)
        // {
        //     cout<<std::fixed<<std::setprecision(17)<<Gyn[index]<<"   ";
        // }
    }

	for(int i=0; i<4; i++)
		flux_res[idx*4 + i] = globaldata[idx].delta * (Gxp[i] + Gxn[i] + Gyp[i]);

}

void interiorindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData, \
	int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn, double* flux_res, double* q, double* max_q, double* min_q, double* dq1, double* dq2)
{
	interior_dGx_pos(globaldata, idx, Gxp, configData, xpos_conn, q, max_q, min_q, dq1, dq2);
	interior_dGx_neg(globaldata, idx, Gxn, configData, xneg_conn, q, max_q, min_q, dq1, dq2);
	interior_dGy_pos(globaldata, idx, Gyp, configData, ypos_conn, q, max_q, min_q, dq1, dq2);
	interior_dGy_neg(globaldata, idx, Gyn, configData, yneg_conn, q, max_q, min_q, dq1, dq2); 


	for(int i=0; i<4; i++)
		flux_res[idx*4 + i] = globaldata[idx].delta * (Gxp[i] + Gxn[i] + Gyp[i] + Gyn[i]);

}
