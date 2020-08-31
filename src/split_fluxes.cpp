#include "split_fluxes.hpp"

void flux_Gyn(double Gyn[4], double nx, double ny, double u1, double u2, double rho, double pr)
{
	double tx = ny;
	double ty = -nx;

	double ut = u1*tx + u2*ty;
	double un = u1*nx + u2*ny;

	double beta = 0.5*rho/pr;
	double S2 = un*sqrt(beta);
	double B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
	double A2neg = 0.5*(1 - erf(S2));

	double pr_by_rho = pr/rho;
	double u_sqr = ut*ut + un*un;

	Gyn[0] = (rho*(un*A2neg - B2));

	double temp1 = pr_by_rho + un*un;
	double temp2 = temp1*A2neg -un*B2;

	temp1 = ut*un*A2neg - ut*B2;
	Gyn[1] = (rho*temp1);

	Gyn[2] = (rho*temp2);

	temp1 = (7.0*pr_by_rho) + u_sqr;
	temp2 = 0.5*un*temp1*A2neg;
	temp1 = (6.0*pr_by_rho) + u_sqr;
	Gyn[3] = (rho*(temp2 - 0.5*temp1*B2));

}