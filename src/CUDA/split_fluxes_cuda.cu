#include "split_fluxes_cuda.hpp"

__device__ void flux_Gxp(double Gxp[4], double nx, double ny, double u1, double u2, double rho, double pr)
{
	double tx = ny;
	double ty = -nx;

	double ut = u1*tx + u2*ty;
	double un = u1*nx + u2*ny;

	double beta = 0.5*rho/pr;
	double S1 = ut*sqrt(beta);
	double B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
	double A1pos = 0.5*(1 + erf(S1));

	double pr_by_rho = pr/rho;
	double u_sqr = ut*ut + un*un;

	Gxp[0] = (rho*(ut*A1pos + B1));

	double temp1 = pr_by_rho + ut*ut;
	double temp2 = temp1*A1pos + ut*B1;

	Gxp[1] = (rho*temp2);

    temp1 = ut*un*A1pos + un*B1;
    Gxp[2] = (rho*temp1);

    temp1 = (7.0*pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1pos;
    temp1 = (6.0*pr_by_rho) + u_sqr;
    Gxp[3] = (rho*(temp2 + 0.5*temp1*B1));

}

__device__ void flux_Gxn(double Gxn[4], double nx, double ny, double u1, double u2, double rho, double pr)
{
	double tx = ny;
	double ty = -nx;

	double ut = u1*tx + u2*ty;
	double un = u1*nx + u2*ny;

	double beta = 0.5*rho/pr;
	double S1 = ut*sqrt(beta);
	double B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
	double A1neg = 0.5*(1 - erf(S1));

	double pr_by_rho = pr/rho;
	double u_sqr = ut*ut + un*un;

	Gxn[0] = (rho*(ut*A1neg - B1));

	double temp1 = pr_by_rho + ut*ut;
	double temp2 = temp1*A1neg - ut*B1;

	Gxn[1] = (rho*temp2);

    temp1 = ut*un*A1neg - un*B1;
    Gxn[2] = (rho*temp1);

    temp1 = (7.0*pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1neg;
    temp1 = (6.0*pr_by_rho) + u_sqr;
    Gxn[3] = (rho*(temp2 - 0.5*temp1*B1));

}

__device__ void flux_Gyn(double Gyn[4], double nx, double ny, double u1, double u2, double rho, double pr)
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

__device__ void flux_Gyp(double Gyp[4], double nx, double ny, double u1, double u2, double rho, double pr)
{
	double tx = ny;
	double ty = -nx;

	double ut = u1*tx + u2*ty;
	double un = u1*nx + u2*ny;

	double beta = 0.5*rho/pr;
	double S2 = un*sqrt(beta);
	double B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
	double A2pos = 0.5*(1 + erf(S2));

	double pr_by_rho = pr/rho;
	double u_sqr = ut*ut + un*un;

	Gyp[0] = (rho*(un*A2pos + B2));

	double temp1 = pr_by_rho + un*un;
	double temp2 = temp1*A2pos +un*B2;

	temp1 = ut*un*A2pos + ut*B2;
	Gyp[1] = (rho*temp1);

	Gyp[2] = (rho*temp2);

	temp1 = (7.0*pr_by_rho) + u_sqr;
	temp2 = 0.5*un*temp1*A2pos;
	temp1 = (6.0*pr_by_rho) + u_sqr;
	Gyp[3] = (rho*(temp2 + 0.5*temp1*B2));

}

__device__ void flux_Gx(double Gx[4], double nx, double ny, double u1, double u2, double rho, double pr)
{
	double tx = ny;
	double ty = -nx;

	double ut = u1*tx + u2*ty;
	double un = u1*nx + u2*ny;

    Gx[0] = rho*ut;

    Gx[1] = pr + rho*ut*ut;

    Gx[2] = rho*ut*un;

    double temp1 = 0.5*(ut*ut + un*un);
    double rho_e = 2.5*pr + rho*temp1;
    Gx[3] = (pr + rho_e)*ut;
}

__device__ void flux_Gy(double Gy[4], double nx, double ny, double u1, double u2, double rho, double pr)
{
	double tx = ny;
	double ty = -nx;

	double ut = u1*tx + u2*ty;
	double un = u1*nx + u2*ny;

    Gy[0] = rho*un;

    Gy[1] = rho*ut*un;

    Gy[2] = pr + rho*un*un;

    double temp1 = 0.5*(ut*ut + un*un);
    double rho_e = 2.5*pr + rho*temp1;
    Gy[3] = (pr + rho_e)*un;
}