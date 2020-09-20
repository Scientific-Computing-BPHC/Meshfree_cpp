#include "state_update_cuda.hpp"

__device__ inline void primitive_to_conserved(double globaldata_prim[4], double nx, double ny, double U[4]);
__device__ inline void conserved_vector_Ubar(double globaldata_prim[4], double nx, double ny, double Mach, double gamma, double pr_inf, double rho_inf, double theta, double Ubar[4]);

__global__ void call_func_delta_cuda(Point* globaldata, int numPoints, double cfl, dim3 thread_dim)
{
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int idx = bx*thread_dim.x + tx;
    
    if(idx < numPoints)
	{
		double min_delt = 1.0;
		for(int i=0; i<20; i++)
		{
			int conn = globaldata[idx].conn[i];
			if (conn == 0) break;

            conn = conn -1; // To account for the indexing difference b/w Julia and C++

			double x_i = globaldata[idx].x;
			double y_i = globaldata[idx].y;
			double x_k = globaldata[conn].x;
			double y_k = globaldata[conn].y;

			double dist = hypot((x_k - x_i), (y_k - y_i));
			double mod_u = hypot(globaldata[conn].prim[1], globaldata[conn].prim[2]);
			double delta_t = dist/(mod_u + 3*sqrt(globaldata[conn].prim[3]/globaldata[conn].prim[0]));
			delta_t *= cfl;
			if (min_delt > delta_t)
				min_delt = delta_t;
		}
		globaldata[idx].delta = min_delt;
		for(int i=0; i<4; i++)
			globaldata[idx].prim_old[i] = globaldata[idx].prim[i];
	}
}

__global__ void state_update_cuda(Point* globaldata, int numPoints, Config configData, int iter, double res_old[1], int rk, int rks, double* res_sqr, dim3 thread_dim)
{
    int bx = blockIdx.x;
    int threadx = threadIdx.x;
    int idx = bx*thread_dim.x + threadx;
    
    double max_res = 0.0;

	double Mach = configData.core.mach;
	double gamma = configData.core.gamma;
	double pr_inf = configData.core.pr_inf;
	double rho_inf = configData.core.rho_inf;
	double theta = configData.core.aoa * (M_PI)/180.0;

    int euler = configData.core.euler;

    double U[4], Uold[4] = {0};

	if(idx < numPoints)
	{
		if(globaldata[idx].flag_1 == 0)
		{
			for(int i=0; i<4; i++)
			{
				U[i] = 0.0;
            }
			state_update_wall(globaldata, idx, max_res, res_sqr, U, Uold, rk, euler);
		}
		else if(globaldata[idx].flag_1 == 2)
		{
			for(int i=0; i<4; i++)
			{
				U[i] = 0.0;
            }
			state_update_outer(globaldata, idx, Mach, gamma, pr_inf, rho_inf, theta, max_res, res_sqr, U, Uold, rk, euler);
		}
		else if(globaldata[idx].flag_1 == 1)
		{
			for(int i=0; i<4; i++)
			{
				U[i] = 0.0;
            }
			state_update_interior(globaldata, idx, max_res, res_sqr, U, Uold, rk, euler);
		}
    }
}

__device__ void state_update_wall(Point* globaldata, int idx, double max_res, double* res_sqr, double U[4], double Uold[4], int rk, int euler)
{
    double nx = globaldata[idx].nx;
    double ny = globaldata[idx].ny;

    primitive_to_conserved(globaldata[idx].prim, nx, ny, U);
    primitive_to_conserved(globaldata[idx].prim_old, nx, ny, Uold);

    double temp = U[0];

    for (int iter=0; iter<4; iter++)
    {
        U[iter] = U[iter] - 0.5 * euler * globaldata[idx].flux_res[iter];
    }

    if (rk == 2)
    {
        for (int iter=0; iter<4; iter++)
            U[iter] = U[iter] * ((double)1.0)/3.0 + Uold[iter] * ((double)2.0)/3.0;
    }

    U[2] = 0.0;
    double U2_rot = U[1];
    double U3_rot = U[2];
    U[1] = U2_rot*ny + U3_rot*nx;
    U[2] = U3_rot*ny - U2_rot*nx;
    res_sqr[idx] = (U[0] - temp)*(U[0] - temp);

    Uold[0] = U[0];
    temp = 1.0 / U[0];
    Uold[1] = U[1]*temp;
    Uold[2] = U[2]*temp;
    Uold[3] = (0.4*U[3]) - ((0.2 * temp) * (U[1] * U[1] + U[2] * U[2]));
    for(int i=0; i<4; i++)
    {
    	globaldata[idx].prim[i] = Uold[i];
    }

}

__device__ void state_update_outer(Point* globaldata, int idx, double Mach, double gamma, double pr_inf, double rho_inf, double theta, double max_res, double* res_sqr, double U[4], double Uold[4], int rk, int euler)
{
    double nx = globaldata[idx].nx;
    double ny = globaldata[idx].ny;

    conserved_vector_Ubar(globaldata[idx].prim, nx, ny, Mach, gamma, pr_inf, rho_inf, theta, U);
    conserved_vector_Ubar(globaldata[idx].prim_old, nx, ny, Mach, gamma, pr_inf, rho_inf, theta, Uold);

    double temp = U[0];
    for (int iter=0; iter<4; iter++)
        U[iter] = U[iter] - 0.5 * euler * globaldata[idx].flux_res[iter];
    if (rk == 2)
    {
        for (int iter=0; iter<4; iter++)
            U[iter] = U[iter] * ((double)1.0)/3.0 + Uold[iter] * ((double)2.0)/3.0;
    }
    //U[2] = 0.0;
    double U2_rot = U[1];
    double U3_rot = U[2];
    U[1] = U2_rot*ny + U3_rot*nx;
    U[2] = U3_rot*ny - U2_rot*nx;
    res_sqr[idx] = (U[0] - temp)*(U[0] - temp);

    Uold[0] = U[0];
    temp = 1.0 / U[0];
    Uold[1] = U[1]*temp;
    Uold[2] = U[2]*temp;
    Uold[3] = (0.4*U[3]) - ((0.2 * temp) * (U[1] * U[1] + U[2] * U[2]));
    for(int i=0; i<4; i++)
    {
    	globaldata[idx].prim[i] = Uold[i];
    }

}

__device__ void state_update_interior(Point* globaldata, int idx, double max_res, double* res_sqr, double U[4], double Uold[4], int rk, int euler)
{
    double nx = globaldata[idx].nx;
    double ny = globaldata[idx].ny;

    primitive_to_conserved(globaldata[idx].prim, nx, ny, U);
    primitive_to_conserved(globaldata[idx].prim_old, nx, ny, Uold);

    double temp = U[0];
    for (int iter=0; iter<4; iter++)
        U[iter] = U[iter] - 0.5 * euler * globaldata[idx].flux_res[iter];
    if (rk == 2)
    {
        for (int iter=0; iter<4; iter++)
            U[iter] = U[iter] * ((double)1.0)/3.0 + Uold[iter] * ((double)2.0)/3.0;
    }

    double U2_rot = U[1];
    double U3_rot = U[2];
    U[1] = U2_rot*ny + U3_rot*nx;
    U[2] = U3_rot*ny - U2_rot*nx;
    res_sqr[idx] = (U[0] - temp)*(U[0] - temp);

    Uold[0] = U[0];
    temp = 1.0 / U[0];
    Uold[1] = U[1]*temp;
    Uold[2] = U[2]*temp;
    Uold[3] = (0.4*U[3]) - ((0.2 * temp) * (U[1] * U[1] + U[2] * U[2]));
    for(int i=0; i<4; i++)
    {
    	globaldata[idx].prim[i] = Uold[i];
    }

}

__device__ inline void primitive_to_conserved(double globaldata_prim[4], double nx, double ny, double U[4])
{
	double rho = globaldata_prim[0];
    U[0] = rho;
    double temp1 = rho * globaldata_prim[1];
    double temp2 = rho * globaldata_prim[2];
    U[1] = temp1*ny - temp2*nx;
    U[2] = temp1*nx + temp2*ny;
    U[3] = 2.5*globaldata_prim[3] + 0.5*(temp1*temp1 + temp2*temp2)/rho;
}

__device__ inline void conserved_vector_Ubar(double globaldata_prim[4], double nx, double ny, double Mach, double gamma, double pr_inf, double rho_inf, double theta, double Ubar[4])
{
	double u1_inf = Mach*cos(theta);
    double u2_inf = Mach*sin(theta);

    double tx = ny;
    double ty = -nx;

    double u1_inf_rot = u1_inf*tx + u2_inf*ty;
    double u2_inf_rot = u1_inf*nx + u2_inf*ny;

    double temp1 = (u1_inf_rot * u1_inf_rot + u2_inf_rot*u2_inf_rot);
    double e_inf = (pr_inf/(rho_inf*(gamma-1))) + 0.5 * (temp1);

    double beta = (0.5 * rho_inf)/pr_inf;
    double S2 = u2_inf_rot * sqrt(beta);
    double B2_inf = exp(-S2*S2)/(2.0*sqrt(M_PI*beta));
    double A2n_inf = 0.5 * (1 - erf(S2));

    double rho = globaldata_prim[0];
    double u1 = globaldata_prim[1];
    double u2 = globaldata_prim[2];
    double pr = globaldata_prim[3];

    double u1_rot = u1*tx + u2*ty;
    double u2_rot = u1*nx + u2*ny;

    temp1 = (u1_rot*u1_rot + u2_rot*u2_rot);
    double e = (pr/(rho*(gamma-1))) + 0.5*(temp1);

    beta = (rho)/(2.0*pr);
    S2 = u2_rot*sqrt(beta);
    double B2 = exp(-S2*S2)/(2.0*sqrt(M_PI*beta));
    double A2p = 0.5*(1.0 + erf(S2));

    Ubar[0] = (rho_inf*A2n_inf) + (rho*A2p);

    Ubar[1] = (rho_inf*u1_inf_rot*A2n_inf) + (rho*u1_rot*A2p);

    temp1 = rho_inf*(u2_inf_rot*A2n_inf - B2_inf);
    double temp2 = rho*(u2_rot*A2p + B2);
    Ubar[2] = (temp1 + temp2);

    temp1 = (rho_inf*A2n_inf* e_inf - 0.5*rho_inf*u2_inf_rot*B2_inf);
    temp2 = (rho*A2p*e + 0.5*rho*u2_rot*B2);

    Ubar[3] = (temp1 + temp2);
}

