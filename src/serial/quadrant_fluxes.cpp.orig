#include "quadrant_fluxes.hpp"

void flux_quad_GxI(double G[4], double nx, double ny, double u1, double u2, double rho, double pr)
{
    double tx = ny;
    double ty = -nx;
    double ut = u1*tx + u2*ty;
    double un = u1*nx + u2*ny;

    double beta = 0.5*rho/pr;
    double S1 = ut*sqrt(beta);
    double S2 = un*sqrt(beta);
    double B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
    double B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
    double A1neg = 0.5*(1.0 - erf(S1));
    double A2neg = 0.5*(1.0 - erf(S2));

    double pr_by_rho = pr/rho;
    double u_sqr = ut*ut + un*un;
    G[0] = rho * A2neg* (ut*A1neg - B1);

    double temp1 = pr_by_rho + ut*ut;
    double temp2 = temp1*A1neg - ut*B1;
    G[1] = rho*A2neg*temp2;

    temp1 = ut*A1neg - B1;
    temp2 = un*A2neg - B2;
    G[2] = rho*temp1*temp2;

    temp1 = (7.0 *pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1neg;

    temp1 = (6.0 *pr_by_rho) + u_sqr;
    double temp3 = 0.5*B1*temp1;

    temp1 = ut*A1neg - B1;
    double temp4 = 0.5*rho*un*B2*temp1;
    G[3] = rho*A2neg*(temp2 - temp3) - temp4;
}

void flux_quad_GxII(double G[4], double nx, double ny, double u1, double u2, double rho, double pr)
{
	double tx = ny;
	double ty = -nx;
	double ut = u1*tx + u2*ty;
	double un = u1*nx + u2*ny;

	double beta = 0.5*rho/pr;
	double S1 = ut*sqrt(beta);
    double S2 = un*sqrt(beta);
    double B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
    double B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
    double A1pos = 0.5*(1.0 + erf(S1));
    double A2neg = 0.5*(1.0 - erf(S2));

    double pr_by_rho = pr/rho;
    double u_sqr = ut*ut + un*un;
    G[0] = rho * A2neg* (ut*A1pos + B1);

    double temp1 = pr_by_rho + ut*ut;
    double temp2 = temp1*A1pos + ut*B1;
    G[1] = rho*A2neg*temp2;

    temp1 = ut*A1pos + B1;
    temp2 = un*A2neg - B2;
    G[2] = rho*temp1*temp2;

    temp1 = (7.0 *pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1pos;

    temp1 = (6.0 * pr_by_rho) + u_sqr;
    double temp3 = 0.5*B1*temp1;

    temp1 = ut*A1pos + B1;
    double temp4 = 0.5*rho*un*B2*temp1;
    G[3] = rho*A2neg*(temp2 + temp3) - temp4;
}

void flux_quad_GxIII(double G[4], double nx, double ny, double u1, double u2, double rho, double pr)
{
    double tx = ny;
    double ty = -nx;
    double ut = u1*tx + u2*ty;
    double un = u1*nx + u2*ny;

    double beta = 0.5*rho/pr;
    double S1 = ut*sqrt(beta);
    double S2 = un*sqrt(beta);
    double B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
    double B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
    double A1pos = 0.5*(1.0 + erf(S1));
    double A2pos = 0.5*(1.0 + erf(S2));

    double pr_by_rho = pr/rho;
    double u_sqr = ut*ut + un*un;
    G[0] = rho * A2pos* (ut*A1pos + B1);

    double temp1 = pr_by_rho + ut*ut;
    double temp2 = temp1*A1pos + ut*B1;
    G[1] = rho*A2pos*temp2;

    temp1 = ut*A1pos + B1;
    temp2 = un*A2pos + B2;
    G[2] = rho*temp1*temp2;

    temp1 = (7.0 *pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1pos;

    temp1 = (6.0 *pr_by_rho) + u_sqr;
    double temp3 = 0.5*B1*temp1;

    temp1 = ut*A1pos - B1;
    double temp4 = 0.5*rho*un*B2*temp1;
    G[3] = rho*A2pos*(temp2 + temp3) + temp4;
}

void flux_quad_GxIV(double G[4], double nx, double ny, double u1, double u2, double rho, double pr)
{
    double tx = ny;
    double ty = -nx;
    double ut = u1*tx + u2*ty;
    double un = u1*nx + u2*ny;

    double beta = 0.5*rho/pr;
    double S1 = ut*sqrt(beta);
    double S2 = un*sqrt(beta);
    double B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
    double B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
    double A1neg = 0.5*(1.0 - erf(S1));
    double A2pos = 0.5*(1.0 + erf(S2));

    double pr_by_rho = pr/rho;
    double u_sqr = ut*ut + un*un;
    G[0] = rho * A2pos* (ut*A1neg - B1);

    double temp1 = pr_by_rho + ut*ut;
    double temp2 = temp1*A1neg - ut*B1;
    G[1] = rho*A2pos*temp2;

    temp1 = ut*A1neg - B1;
    temp2 = un*A2pos + B2;
    G[2] = rho*temp1*temp2;

    temp1 = (7.0 *pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1neg;

    temp1 = (6.0 *pr_by_rho) + u_sqr;
    double temp3 = 0.5*B1*temp1;

    temp1 = ut*A1neg - B1;
    double temp4 = 0.5*rho*un*B2*temp1;
    G[3] = rho*A2pos*(temp2 - temp3) + temp4;
}