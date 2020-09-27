#include "point.hpp"
#include "limiters_cuda.hpp"

__device__ inline void update_qtildes(double qtilde[4], double q[4], double dq1[4], double dq2[4], double delta_x, double delta_y);
__device__ inline void update_qtildes(double qtilde[4], double q[4], double dq1[4], double dq2[4], double delta_x, double delta_y, double phi[4]);

__device__ void venkat_limiter(double qtilde[4], double vl_const, Point* globaldata, int index, double gamma, double phi[4])
{
	double ds = globaldata[index].short_distance;
	double epsi = vl_const * ds;
	epsi = epsi * epsi * epsi;
	double del_pos = 0.0;
	double del_neg = 0.0;
	VLBroadcaster(globaldata[index].q, qtilde, globaldata[index].max_q, globaldata[index].min_q, phi, epsi, del_pos, del_neg);
}

__device__ void VLBroadcaster(double q[4], double qtilde[4], double max_q[4], double min_q[4], double phi[4], double epsi, double del_pos, double del_neg)
{
	for(int i=0; i<4; i++)
	{
		del_neg = qtilde[i + 4*threadIdx.x] - q[i];
		if(abs(del_neg) <= 1e-5)
			phi[i] = 1.0;
		else if (abs(del_neg) > 1e-5)
		{
			if (del_neg > 0)
				del_pos = max_q[i] - q[i];
			else if (del_neg < 0)
				del_pos = min_q[i] - q[i];

			double num = (del_pos*del_pos) + (epsi*epsi);
			num = (num*del_neg) + 2 * (del_neg*del_neg*del_pos);

			double den = (del_pos*del_pos) + (2*del_neg*del_neg);
			den = den + (del_neg*del_pos) + (epsi*epsi);
			den = den*del_neg;

			double temp = num/den;
			if (temp<1.0)
				phi[i] = temp;
			else
				phi[i] = 1.0;
		}

	}
}

__device__ void connectivity_stats(double x_i, double y_i, double nx, double ny, double power, double conn_x, double conn_y, double& sig_del_x_sqr, double& sig_del_y_sqr, double& sig_del_x_del_y, double& delta_x, double& delta_y, double& delta_s_weights, double& delta_n_weights)
{
    double x_k = conn_x;
    double y_k = conn_y;
    
    delta_x = x_k - x_i;
    delta_y = y_k - y_i;
    
    double delta_s = delta_x*ny - delta_y*nx;
    double delta_n = delta_x*nx + delta_y*ny;

    double dist = sqrt(delta_s*delta_s + delta_n*delta_n);
    double weights = pow(dist, power);
    
    delta_s_weights = delta_s*weights;
    delta_n_weights = delta_n*weights;
    
    sig_del_x_sqr += (delta_s*delta_s_weights);
    sig_del_y_sqr += (delta_n*delta_n_weights);
	sig_del_x_del_y += (delta_s*delta_n_weights);
	
 }

 __device__ void calculate_qtile(double qtilde_i[4], double qtilde_k[4], Point* globaldata, int idx, int conn, double delta_x, double delta_y, double vl_const, double gamma, int limiter_flag, double phi_i[4], double phi_k[4])
{
	update_qtildes(qtilde_i, globaldata[idx].q, globaldata[idx].dq1, globaldata[idx].dq2, delta_x, delta_y);
	update_qtildes(qtilde_k, globaldata[conn].q, globaldata[conn].dq1, globaldata[conn].dq2, delta_x, delta_y);

	if(limiter_flag == 1)
	{
		venkat_limiter(qtilde_i, vl_const, globaldata, idx, gamma, phi_i);
		venkat_limiter(qtilde_k, vl_const, globaldata, conn, gamma, phi_k);
		update_qtildes(qtilde_i, globaldata[idx].q, globaldata[idx].dq1, globaldata[idx].dq2, delta_x, delta_y, phi_i);
		update_qtildes(qtilde_k, globaldata[conn].q, globaldata[conn].dq1, globaldata[conn].dq2, delta_x, delta_y, phi_k);
	}
}



__device__ inline void update_qtildes(double qtilde[4], double q[4], double dq1[4], double dq2[4], double delta_x, double delta_y)
{
	for(int iter=0; iter<4; iter++)
	{
		qtilde[iter + 4*threadIdx.x] = q[iter] - 0.5 * (delta_x * dq1[iter] + delta_y * dq2[iter]);
	}
}

__device__ inline void update_qtildes(double qtilde[4], double q[4], double dq1[4], double dq2[4], double delta_x, double delta_y, double phi[4])
{
	for(int iter=0; iter<4; iter++)
	{
		qtilde[iter + 4*threadIdx.x] = q[iter] - 0.5 * phi[iter] * (delta_x * dq1[iter] + delta_y * dq2[iter]);
	}
}

__device__ void update_delf(double sig_del_x_del_f[4], double sig_del_y_del_f[4], double G_k[4], double G_i[4], double delta_s_weights, double delta_n_weights)
{
	for(int iter=0; iter<4; iter++)
	{
		double intermediate_var = G_k[iter] - G_i[iter];
		sig_del_x_del_f[iter + 4*threadIdx.x] += (intermediate_var * delta_s_weights);
		sig_del_y_del_f[iter + 4*threadIdx.x] += (intermediate_var * delta_n_weights);
	}
}

__device__ void qtilde_to_primitive(double result[4], double qtilde[4], double gamma)
{
    double beta = -qtilde[3 + 4*threadIdx.x]*0.5;
    double temp = 0.5/beta;
    double u1 = qtilde[1 + 4*threadIdx.x]*temp;
    double u2 = qtilde[2 + 4*threadIdx.x]*temp;

    double temp1 = qtilde[0 + 4*threadIdx.x] + beta*(u1*u1 + u2*u2);
    double temp2 = temp1 - (log(beta)/(gamma-1));
    double rho = exp(temp2);
    double pr = rho*temp;
    result[0 + 4*threadIdx.x] = u1;
    result[1 + 4*threadIdx.x] = u2;
    result[2 + 4*threadIdx.x] = rho;
    result[3 + 4*threadIdx.x] = pr;
}
