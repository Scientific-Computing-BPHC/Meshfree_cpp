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

void cal_flux_residual(Point* globaldata, int numPoints, Config configData)
{

	double Gxp[4] = {0}, Gxn[4] = {0}, Gyp[4] = {0}, Gyn[4] = {0};

	for(int idx=0; idx<numPoints; idx++)
	{
		

		if (globaldata[idx].flag_1 == 0)
			wallindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, configData);
		else if (globaldata[idx].flag_1 == 2)
			outerindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, configData);
		else if (globaldata[idx].flag_1 == 1)
			interiorindices_flux_residual(globaldata, idx, Gxp, Gxn, Gyp, Gyn, configData);

	}
}

void wallindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData)
{


	wall_dGx_pos(globaldata, idx, Gxp, configData);
	wall_dGx_neg(globaldata, idx, Gxn, configData);
	wall_dGy_neg(globaldata, idx, Gyn, configData);

	for(int i=0; i<4; i++)
	{
		Gxp[i] = (Gxp[i] + Gxn[i] + Gyn[i]) * 2;
	}

	for(int i=0; i<4; i++)
	{
		globaldata[idx].flux_res[i] = Gxp[i];
	}
}

void outerindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData)
{

	outer_dGx_pos(globaldata, idx, Gxp, configData);
	outer_dGx_neg(globaldata, idx, Gxn, configData);
	outer_dGy_pos(globaldata, idx, Gyp, configData);

	for(int i=0; i<4; i++)
		Gxp[i] = (Gxp[i] + Gxn[i] + Gyp[i]);

	for(int i=0; i<4; i++)
		globaldata[idx].flux_res[i] = Gxp[i];

}

void interiorindices_flux_residual(Point* globaldata, int idx, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], Config configData)
{
	interior_dGx_pos(globaldata, idx, Gxp, configData);
	interior_dGx_neg(globaldata, idx, Gxn, configData);
	interior_dGy_pos(globaldata, idx, Gyp, configData);
	interior_dGy_neg(globaldata, idx, Gyn, configData); 

	for(int i=0; i<4; i++)
		Gxp[i] = (Gxp[i] + Gxn[i] + Gyp[i] + Gyn[i]);

	for(int i=0; i<4; i++)
		globaldata[idx].flux_res[i] = Gxp[i];
}
