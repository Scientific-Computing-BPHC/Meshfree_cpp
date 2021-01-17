/* External Headers */
#include<thrust/reduce.h>
#include <thrust/system/cuda/execution_policy.h>

/* Internal Headers */
#include "core_cuda.hpp"
#include "point.hpp"
#include "state_update_cuda.hpp"
#include "flux_residual_cuda.hpp"
#include "utils.hpp"

__device__ inline void q_var_derivatives_get_sum_delq_innerloop(Point* globaldata, int idx, int conn, double weights, double delta_x, double delta_y, \
    double qi_tilde[4], double qk_tilde[4], double sig_del_x_del_q[4], double sig_del_y_del_q[4], double* q, double* dq1, double* dq2);

template <class Type>
bool isNan(Type var)
{
    if(var!=var) return true;
    return false;
}

double calculateTheta(Config configData)
{
    return (configData.core.aoa * (M_PI)/180.0);
}

void getInitialPrimitive(Config configData, double primal[4])
{
	primal[0] = configData.core.rho_inf;
	double mach = configData.core.mach;
	double machcos = mach * cos(calculateTheta(configData));
	double machsin = mach * sin(calculateTheta(configData));
	primal[1] = machcos;
	primal[2] = machsin;
	primal[3] = configData.core.pr_inf;

}

void placeNormals(Point* globaldata, int idx, Config configData, long long interior, long long wall, long long outer)
{
	int flag = globaldata[idx].flag_1;
    if (flag == wall || flag == outer)
    {
        xy_tuple currpt = getxy(globaldata[idx]);
        int leftpt_tmp = globaldata[idx].left;
        leftpt_tmp = leftpt_tmp - 1; // To account for indexing
        xy_tuple leftpt = getxy(globaldata[leftpt_tmp]);
        int rightpt_tmp = globaldata[idx].right;
        rightpt_tmp = rightpt_tmp - 1; // To account for indexing
        xy_tuple rightpt = getxy(globaldata[rightpt_tmp]);
        xy_tuple normals = calculateNormals(leftpt, rightpt, std::get<0>(currpt), std::get<1>(currpt));
        setNormals(globaldata, idx, normals);
     }
    else if (flag == interior)
     	setNormals(globaldata, idx, std::make_tuple(0.0, 1.0));
    else
    	cout<<"Illegal Point Type"<<endl;
}

xy_tuple calculateNormals(xy_tuple left, xy_tuple right, double mx, double my)
{
	double lx = std::get<0>(left);
	double ly = std::get<1>(left);
	double rx = std::get<0>(right);
	double ry = std::get<1>(right);
	double nx1 = my - ly;
	double nx2 = ry - my;
	double ny1 = mx - lx;
	double ny2 = rx - mx;
	double nx = 0.5*(nx1 + nx2);
	double ny = 0.5*(ny1 + ny2);
	double det = hypot(nx, ny);
	nx = -nx/det;
	ny = ny/det;
	return std::make_tuple(nx, ny);
}

void calculateConnectivity(Point* globaldata, int idx, int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn, int* connec)
{
	Point ptInterest = globaldata[idx];
	double currx = ptInterest.x;
    double curry = ptInterest.y;
    double nx = ptInterest.nx;
    double ny = ptInterest.ny;
    int flag = ptInterest.flag_1;
    double tx = ny;
    double ty = -nx;
    int xpos_nbhs = 0;
    int xneg_nbhs = 0;
    int ypos_nbhs = 0;
    int yneg_nbhs = 0; 
    int xpos_conn_tmp[20] = {0};
    int ypos_conn_tmp[20] = {0};
    int xneg_conn_tmp[20] = {0};
    int yneg_conn_tmp[20] = {0};
    
    /* Start Connectivity Generation */
    for (int i=0; i<20; i++)
    {
    	int itm = connec[idx*20 +i];
    	if (itm==0) 
    	{
    		break;
    	}
        itm = itm -1; // to account for indexing
    	double itmx = globaldata[itm].x;
    	double itmy = globaldata[itm].y;
    	double delta_x = itmx - currx;
    	double delta_y = itmy - curry;
    	double delta_s = delta_x*tx + delta_y*ty;
    	double delta_n = delta_x*nx + delta_y*ny;
        itm = itm + 1; // to reaccount for indexing when we add the point below xpos_conn[xpos_nbhs] = itm;
    	if(delta_s <= 0.0)
    	{
    		xpos_conn_tmp[xpos_nbhs] = itm;
            xpos_nbhs+=1;
    	}
    	if(delta_s >= 0.0)
    	{
    		xneg_conn_tmp[xneg_nbhs] = itm;
            xneg_nbhs+=1;
    	}
    	if(flag==1)
    	{
    		if(delta_n<=0.0)
    		{
    			ypos_conn_tmp[ypos_nbhs] = itm;
                ypos_nbhs+=1;
    		}
    		if(delta_n>=0.0)
    		{	
    			yneg_conn_tmp[yneg_nbhs] = itm;
                yneg_nbhs+=1;
    		}
    	}
    	else if (flag==0)
    	{
    		yneg_conn_tmp[yneg_nbhs] = itm;
            yneg_nbhs+=1;
    	}
    	else if (flag==2)
    	{
    		ypos_conn_tmp[ypos_nbhs] = itm;
            ypos_nbhs+=1;
    	}
    }
    /* End Connectivity Generation */

    for(int i=0; i<20; i++)
    {
    	xpos_conn[idx*20 + i] = xpos_conn_tmp[i];
    	xneg_conn[idx*20 + i] = xneg_conn_tmp[i];
    	ypos_conn[idx*20 + i] = ypos_conn_tmp[i];
    	yneg_conn[idx*20 + i] = yneg_conn_tmp[i];
    }
    globaldata[idx].xpos_nbhs = xpos_nbhs;
    globaldata[idx].xneg_nbhs = xneg_nbhs;
    globaldata[idx].ypos_nbhs = ypos_nbhs;
    globaldata[idx].yneg_nbhs = yneg_nbhs;	
}

void fpi_solver(int iter, Point* globaldata_d, Config configData, double* res_old_d, double* res_sqr_d, int numPoints, \
    TempqDers* tempdq_d, cudaStream_t stream, double res_old[1], double* res_sqr, unsigned int mem_size_C, unsigned int mem_size_D, \
    int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn, \
    int* connec_d, double* prim_d, double* flux_res_d, double* q_d, double* dq1_d, double* dq2_d, double* max_q_d, double* min_q_d, double* prim_old_d)
{
    if (iter == 0) cout<<"\nStarting FuncDelta"<<endl;
    int block_size = configData.core.threadsperblock;
    dim3 threads(block_size);
    dim3 grid((numPoints / threads.x +1));
    int rks = configData.core.rks;
    double cfl = configData.core.cfl;
    double power = configData.core.power;

    call_func_delta_cuda<<<grid, threads, 0, stream>>>(globaldata_d, numPoints, cfl, threads, connec_d, prim_d, prim_old_d);
    cudaDeviceSynchronize();

    for(int rk=0; rk<rks; rk++)
    {
        call_rem_fpi_solver_cuda(globaldata_d, numPoints, power, tempdq_d, block_size, configData, res_old_d, res_sqr_d, iter, rk, rks, threads, \
            grid, stream, res_old, res_sqr, mem_size_C, mem_size_D, xpos_conn, xneg_conn, ypos_conn, yneg_conn, \
            connec_d, prim_d, flux_res_d, q_d, dq1_d, dq2_d, max_q_d, min_q_d, prim_old_d);
    }
}

void call_rem_fpi_solver_cuda(Point* globaldata_d, int numPoints, double power, TempqDers* tempdq_d, int block_size, Config configData, double* res_old_d, \
    double* res_sqr_d, int iter, int rk, int rks, dim3 threads, dim3 grid, cudaStream_t stream, double res_old[1], \
    double* res_sqr, unsigned int mem_size_C, unsigned int mem_size_D, int* xpos_conn_d, int* xneg_conn_d, int* ypos_conn_d, int* yneg_conn_d,\
    int* connec_d, double* prim_d, double* flux_res_d, double* q_d, double* dq1_d, double* dq2_d, double* max_q_d, double* min_q_d, double* prim_old_d)
{

    /* Make the kernel calls */
    q_variables_cuda<<<grid, threads, 0, stream>>>(globaldata_d, numPoints, threads, prim_d, q_d);
    q_var_derivatives_cuda<<<grid, threads, 0, stream>>>(globaldata_d, numPoints, power, threads, q_d, dq1_d, dq2_d, connec_d, max_q_d, min_q_d);
    for(int inner_iters=0; inner_iters<2; inner_iters++) // Basically, three inner iters
    {
        q_var_derivatives_innerloop_cuda<<<grid, threads, 0, stream>>>(globaldata_d, numPoints, power, tempdq_d, threads, connec_d, q_d, dq1_d, dq2_d);
        q_var_derivatives_update_innerloop_cuda<<<grid, threads, 0, stream>>>(globaldata_d, tempdq_d, threads, dq1_d, dq2_d);
    }
    cal_flux_residual_cuda<<<grid, threads, 0, stream>>>(globaldata_d, numPoints, configData, threads, xpos_conn_d, xneg_conn_d, ypos_conn_d, yneg_conn_d, flux_res_d,\
    q_d, max_q_d, min_q_d, dq1_d, dq2_d);
    checkCudaErrors(cudaMemcpyAsync(res_old_d, res_old, mem_size_C, cudaMemcpyHostToDevice, stream)); 
    checkCudaErrors(cudaMemcpyAsync(res_sqr_d, res_sqr, mem_size_D, cudaMemcpyHostToDevice, stream)); 
    
    state_update_cuda<<<grid, threads, 0, stream>>>(globaldata_d, numPoints, configData, iter, res_old_d, rk, rks, res_sqr_d, threads, prim_d, prim_old_d, flux_res_d);

    cudaDeviceSynchronize();
    checkCudaErrors(cudaMemcpyAsync(res_old, res_old_d, mem_size_C, cudaMemcpyDeviceToHost, stream));

    double sig_res_sqr = thrust::reduce(thrust::cuda::par.on(stream), res_sqr_d, res_sqr_d + numPoints, (double) 0.0, thrust::plus<double>());
    double res_new = sqrt(sig_res_sqr)/numPoints;
	double residue = 0.0;
	if(iter<=1)   
	{
		res_old[0] = res_new;
		residue = 0.0;
	}
	else
		residue = log10(res_new/res_old[0]);
	if(rk == rks-1)
		cout<<std::fixed<<std::setprecision(17)<<"\n Residue: "<<iter+1<<" "<<residue<<endl;
}

__global__ void q_variables_cuda(Point* globaldata, int numPoints, dim3 thread_dim, double* prim, double* q)
{
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int idx = bx*thread_dim.x + tx;
    double q_result[4] = {0};
    if(idx < numPoints)
    {
        double rho = prim[idx*4 + 0];
        double u1 = prim[idx*4 + 1];
        double u2 = prim[idx*4 + 2];
        double pr = prim[idx*4 + 3];
        double beta = 0.5 * (rho/pr);
        double two_times_beta = 2.0 * beta;
        q_result[0] = log(rho) + log(beta) * 2.5 - (beta * ((u1 * u1) + (u2 * u2)));
        q_result[1] = (two_times_beta * u1);
        q_result[2] = (two_times_beta * u2);
        q_result[3] = -two_times_beta;
        for(int i=0; i<4; i++)
        {
            q[idx*4 + i] = q_result[i];
        }
    }
}

__global__ void q_var_derivatives_cuda(Point* globaldata, int numPoints, double power, dim3 thread_dim, double* q, double* dq1, double* dq2, int* connec, double* max_q, double* min_q)
{
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int idx = bx*thread_dim.x + tx;
    double sig_del_x_del_q[4], sig_del_y_del_q[4], min_q_tmp[4], max_q_tmp[4];
    if(idx < numPoints)
    {
        double x_i = globaldata[idx].x;
        double y_i = globaldata[idx].y;
        double sig_del_x_sqr = 0.0;
        double sig_del_y_sqr = 0.0;
        double sig_del_x_del_y = 0.0;
        #pragma unroll
        for(int i=0; i<4; i++)
        {
            sig_del_x_del_q[i] = 0.0;
            sig_del_y_del_q[i] = 0.0;
        }
        #pragma unroll
        for(int i=0; i<4; i++)
        {
            max_q_tmp[i] = q[idx*4 + i];
            min_q_tmp[i] = q[idx*4 + i];
        }
        #pragma unroll
        for(int i=0; i<20; i++)
        {
            int conn = connec[idx*20 + i];
            if(conn == 0) 
            {
                break;
            }
            conn = conn - 1; // To account for the indexing difference
            double x_k = globaldata[conn].x;
            double y_k = globaldata[conn].y;
            double delta_x = x_k - x_i;
            double delta_y = y_k - y_i;
            double dist = hypot(delta_x, delta_y);
            double weights = pow(dist, power);
            sig_del_x_sqr += ((delta_x * delta_x) * weights);
            sig_del_y_sqr += ((delta_y * delta_y) * weights);
            sig_del_x_del_y += ((delta_x * delta_y) * weights);
            #pragma unroll
            for(int iter=0; iter<4; iter++)
            {
                double intermediate_var = weights * (q[conn*4 + iter] - q[idx*4 + iter]);
                sig_del_x_del_q[iter] = sig_del_x_del_q[iter] + (delta_x * intermediate_var);
                sig_del_y_del_q[iter] = sig_del_y_del_q[iter] + (delta_y * intermediate_var);

            }
            #pragma unroll
            for(int j=0; j<4; j++)
            {
                if (max_q_tmp[j] < q[conn*4 + j])
                {
                    max_q_tmp[j] = q[conn*4 + j];
                }
                if(min_q_tmp[j] > q[conn*4 + j])
                {
                    min_q_tmp[j] = q[conn*4 + j];
                }
            }
        }
        #pragma unroll
        for(int i=0; i<4; i++)
        {
            max_q[idx*4 + i] = max_q_tmp[i];
            min_q[idx*4 + i] = min_q_tmp[i];
        }
        double det = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
        double one_by_det = 1.0/det;
        #pragma unroll
        for(int iter=0; iter<4; iter++)
        {
            dq1[idx*4 + iter] = one_by_det * (sig_del_x_del_q[iter] * sig_del_y_sqr - sig_del_y_del_q[iter] * sig_del_x_del_y);
            dq2[idx*4 + iter] = one_by_det * (sig_del_y_del_q[iter] * sig_del_x_sqr - sig_del_x_del_q[iter] * sig_del_x_del_y);
        }
    }
}

__global__ void q_var_derivatives_innerloop_cuda(Point* globaldata, int numPoints, double power, TempqDers* tempdq, dim3 thread_dim, int* connec, double* q, double* dq1, double* dq2)
{   
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int idx = bx*thread_dim.x + tx;
    double sig_del_x_del_q[4], sig_del_y_del_q[4], qi_tilde[4] ={0}, qk_tilde[4] = {0};
    if(idx <numPoints)
    {
        double x_i = globaldata[idx].x;
        double y_i = globaldata[idx].y;
        double sig_del_x_sqr = 0.0;
        double sig_del_y_sqr = 0.0;
        double sig_del_x_del_y = 0.0;
        #pragma unroll
        for(int i=0; i<4; i++)
        {
            sig_del_x_del_q[i] = 0.0;
            sig_del_y_del_q[i] = 0.0;
        }
        #pragma unroll
        for(int i=0; i<20; i++)
        {
            int conn = connec[idx*20 + i];
            if(conn == 0) break;
            conn = conn - 1;
            double x_k = globaldata[conn].x;
            double y_k = globaldata[conn].y;
            double delta_x = x_k - x_i;
            double delta_y = y_k - y_i;
            double dist = hypot(delta_x, delta_y);
            double weights = pow(dist, power);
            sig_del_x_sqr += ((delta_x * delta_x) * weights);
            sig_del_y_sqr += ((delta_y * delta_y) * weights);
            sig_del_x_del_y += ((delta_x * delta_y) * weights);

            q_var_derivatives_get_sum_delq_innerloop(globaldata, idx, conn, weights, delta_x, delta_y, qi_tilde, qk_tilde, sig_del_x_del_q, sig_del_y_del_q, q, dq1, dq2);
        }
        double det = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
        double one_by_det = 1.0/det;
        #pragma unroll
        for(int iter =0; iter<4; iter++)
        {   
            tempdq[idx].dq1[iter] = one_by_det * (sig_del_x_del_q[iter] * sig_del_y_sqr - sig_del_y_del_q[iter] * sig_del_x_del_y);
            tempdq[idx].dq2[iter] = one_by_det * (sig_del_y_del_q[iter] * sig_del_x_sqr - sig_del_x_del_q[iter] * sig_del_x_del_y);
        }
    }
}

__device__ inline void q_var_derivatives_get_sum_delq_innerloop(Point* globaldata, int idx, int conn, double weights, double delta_x, double delta_y, \
    double qi_tilde[4], double qk_tilde[4], double sig_del_x_del_q[4], double sig_del_y_del_q[4], double* q, double* dq1, double* dq2)
{
   #pragma unroll 
   for(int iter=0; iter<4; iter++)
    {

        qi_tilde[iter] = q[idx*4 + iter] - 0.5 * (delta_x * dq1[idx*4 + iter] + delta_y * dq2[idx*4 + iter]);
        qk_tilde[iter] = q[conn*4 + iter] - 0.5 * (delta_x * dq1[conn*4 + iter] + delta_y * dq2[conn*4 + iter]);
        double intermediate_var = weights * (qk_tilde[iter] - qi_tilde[iter]);
        sig_del_x_del_q[iter] += (delta_x * intermediate_var);
        sig_del_y_del_q[iter] += (delta_y * intermediate_var);
    }
}

__global__ void q_var_derivatives_update_innerloop_cuda(Point* globaldata, TempqDers* tempdq, dim3 thread_dim, double* dq1, double* dq2)
{
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int idx = bx*thread_dim.x + tx;
    #pragma unroll 
    for(int iter=0; iter<4; iter++)
    {
        dq1[idx*4 + iter] = tempdq[idx].dq1[iter];
        dq2[idx*4 + iter] = tempdq[idx].dq2[iter];
    }
}



