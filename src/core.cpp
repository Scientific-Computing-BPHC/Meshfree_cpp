#include "core.hpp"
#include "point.hpp"
#include "state_update.hpp"
#include "flux_residual.hpp"

inline void q_var_derivatives_update(double sig_del_x_sqr, double sig_del_y_sqr, double sig_del_x_del_y, double sig_del_x_del_q[4], double sig_del_y_del_q[4], double dq1_store[4], double dq2_store[2]);
inline void q_var_derivatives_get_sum_delq_innerloop(Point* globaldata, int idx, int conn, double weights, double delta_x, double delta_y, double qi_tilde[4], double qk_tilde[4], double sig_del_x_del_q[4], double sig_del_y_del_q[4]);
inline void q_var_derivatives_update_innerloop(double dq1[4], double dq2[4], int idx, double tempdq[][2][4]);


inline double deg2rad(double radians) {
    return radians * (180.0 / M_PI);
}

double calculateTheta(Config configData)
{
	double theta = deg2rad(configData.core.aoa);
	return theta;
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
        xy_tuple leftpt = getxy(globaldata[leftpt_tmp]);
        int rightpt_tmp = globaldata[idx].right;
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
	ny = -ny/det;

	return std::make_tuple(nx, ny);
}

void calculateConnectivity(Point* globaldata, int idx)
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
    int xpos_conn[20] = {0};
    int ypos_conn[20] = {0};
    int xneg_conn[20] = {0};
    int yneg_conn[20] = {0};

    // /* Start Connectivity Generation */
    for (int i=0; i<20; i++)
    {
    	int itm = ptInterest.conn[i];
    	if (itm==0) 
    	{
    		//cout<<"\n Breaking"<<endl;
    		break;
    	}

    	//cout<< "\n Unbroken \n";
    	double itmx = globaldata[itm].x;
    	double itmy = globaldata[itm].y;

    	double delta_x = itmx - currx;
    	double delta_y = itmy - curry;

    	double delta_s = delta_x*tx + delta_y*ty;
    	double delta_n = delta_x*nx + delta_y*ny;

    	if(delta_s <= 0.0)
    	{
    		xpos_nbhs+=1;
    		xpos_conn[xpos_nbhs] = itm;
    	}

    	if(delta_s >= 0.0)
    	{
    		xneg_nbhs+=1;
    		xneg_conn[xpos_nbhs] = itm;
    	}

    	if(flag==1)
    	{
    		if(delta_n<=0.0)
    		{
    			ypos_nbhs+=1;
    			ypos_conn[ypos_nbhs] = itm;
    		}

    		if(delta_n>=0.0)
    		{
    			yneg_nbhs+=1;
    			yneg_conn[yneg_nbhs] = itm;
    		}
    	}

    	else if (flag==0)
    	{
    		yneg_nbhs+=1;
    		yneg_conn[yneg_nbhs] = itm;
    	}

    	else if (flag==2)
    	{
    		ypos_nbhs+=1;
    		ypos_conn[ypos_nbhs] = itm;
    	}
    }
    /* End Connectivity Generation */

    //cout<<"\n Just Checking: "<<xpos_conn[0]<<endl;

    for(int i=0; i<20; i++)
    {
    	globaldata[idx].xpos_conn[i] = xpos_conn[i];
    	globaldata[idx].xneg_conn[i] = xneg_conn[i];
    	globaldata[idx].ypos_conn[i] = ypos_conn[i];
    	globaldata[idx].yneg_conn[i] = yneg_conn[i];
    }

    globaldata[idx].xpos_nbhs = xpos_nbhs;
    globaldata[idx].xneg_nbhs = xneg_nbhs;
    globaldata[idx].ypos_nbhs = ypos_nbhs;
    globaldata[idx].yneg_nbhs = yneg_nbhs;	

}

void fpi_solver(int iter, Point* globaldata, Config configData, double res_old[1], int numPoints, double main_store[62], double tempdq[][2][4])
{
    if (iter == 0)
        cout<<"\nStarting FuncDelta"<<endl;

    double power = main_store[52];
    double cfl = main_store[53];

    func_delta(globaldata, numPoints, cfl);

    double phi_i[4], phi_k[4], G_i[4], G_k[4], result[4], qtilde_i[4], qtilde_k[4];
    double Gxp[4], Gxn[4], Gyp[4], Gyn[4], sig_del_x_del_f[4], sig_del_y_del_f[4];

    for(int i=0; i<4; i++)
        phi_i[i] = main_store[i];
    for(int i=4; i<8; i++)
        phi_k[i-4] = main_store[i];
    for(int i=8; i<12; i++)
        G_i[i-8] = main_store[i];
    for(int i=12; i<16; i++)
        G_k[i-12] = main_store[i];
    for(int i=16; i<20; i++)
        result[i-16] = main_store[i];
    for(int i=20; i<24; i++)
        qtilde_i[i-20] = main_store[i];
    for(int i=24; i<28; i++)
        qtilde_k[i-24] = main_store[i];
    for(int i=28; i<32; i++)
        Gxp[i-28] = main_store[i];
    for(int i=32; i<36; i++)
        Gxn[i-32] = main_store[i];
    for(int i=36; i<40; i++)
        Gyp[i-36] = main_store[i];
    for(int i=40; i<44; i++)
        Gyn[i-40] = main_store[i];
    for(int i=44; i<48; i++)
        sig_del_x_del_f[i-44] =main_store[i];
    for(int i=48; i<52; i++)
        sig_del_y_del_f[i-48] = main_store[i];

    cout<<"\nIteration Number: "<<iter+1<<endl;

    for(int rk=0; rk<4; rk++)
    {
        
        cout<<"\n rk: "<<rk<<"\n";

        q_variables(globaldata, numPoints, result);

        //debug_globaldata(globaldata, numPoints, iter, rk);

        q_var_derivatives(globaldata, numPoints, power, tempdq, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k);



        for(int inner_iters=0; inner_iters<3; inner_iters++)
        {
            q_var_derivatives_innerloop(globaldata, numPoints, power, tempdq, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k);
        }

        cout<<"\nCalculating Flux Residual\n";

        cal_flux_residual(globaldata, numPoints, configData, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k,
            result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, main_store);

        printDebug(globaldata, numPoints, configData, iter, res_old, rk, sig_del_x_del_f, sig_del_y_del_f, main_store);

        cout<<"\nDone Calculating Flux Residual\n";

        state_update(globaldata, numPoints, configData, iter, res_old, rk, sig_del_x_del_f, sig_del_y_del_f, main_store);
    }
}

void q_variables(Point* globaldata, int numPoints, double q_result[4])
{
    cout<<"\nInside q_variables\n";
    //cout<<"\nSecond Check for No. Of Points: "<<numPoints<<endl;

    for(int idx=0; idx<numPoints; idx++)
    {
        double rho = globaldata[idx].prim[0];
        double u1 = globaldata[idx].prim[1];
        double u2 = globaldata[idx].prim[2];
        double pr = globaldata[idx].prim[3];
        double beta = 0.5 * ((double)rho/pr);
        double two_times_beta = 2.0 * beta;
        q_result[0] = log(rho) + log(beta) * 2.5 - (beta * ((u1 * u1) + (u2 * u2)));
        q_result[1] = (two_times_beta * u1);
        q_result[2] = (two_times_beta * u2);
        q_result[3] = -two_times_beta;
        for(int i=0; i<4; i++)
        {
            globaldata[idx].q[i] = q_result[i];
            //cout<<"q_result: "<<i<<" "<<q_result[i]<<endl;
        }
    }
    cout<<"\nGoing outta q_variables\n";
}

void q_var_derivatives(Point* globaldata, int numPoints, double power, double tempdq[][2][4], double sig_del_x_del_q[4], double sig_del_y_del_q[4], double max_q[4], double min_q[4])
{
    cout<<"\nInside q_var_derivatives\n";

    for(int idx=0; idx<numPoints; idx++)
    {
        double x_i = globaldata[idx].x;
        double y_i = globaldata[idx].y;
        double sig_del_x_sqr = 0.0;
        double sig_del_y_sqr = 0.0;
        double sig_del_x_del_y = 0.0;

        for(int i=0; i<4; i++)
        {
            sig_del_x_del_q[i] = 0.0;
            sig_del_y_del_q[i] = 0.0;
        }

        for(int i=0; i<20; i++)
        {
            int conn = globaldata[idx].conn[i];
            if(conn == 0) break;
            double x_k = globaldata[conn].x;
            double y_k = globaldata[conn].y;

            double delta_x = x_k - x_i;
            double delta_y = y_k - y_i;

            double dist = hypot(delta_x, delta_y);
            double weights = pow(dist, power);
            sig_del_x_sqr += (delta_x * delta_x) * weights;
            sig_del_y_sqr += (delta_y * delta_y) * weights;
            sig_del_x_del_y += (delta_x * delta_y) * weights;

            for(int iter=0; iter<4; iter++)
            {
                double intermediate_var = weights * (globaldata[conn].q[iter] - globaldata[idx].q[iter]);
                sig_del_x_del_q[iter] += delta_x * intermediate_var;
                sig_del_y_del_q[iter] += delta_y * intermediate_var;
            }

            for(int i=0; i<4; i++)
            {
                if (max_q[i] < globaldata[conn].q[i])
                    max_q[i] = globaldata[conn].q[i];
                if(min_q[i] > globaldata[conn].q[i])
                    min_q[i] = globaldata[conn].q[i];
            }

            for(int i=0; i<4; i++)
            {
                globaldata[idx].max_q[i] = max_q[i];
                globaldata[idx].min_q[i] = min_q[i];
            }

        }

        q_var_derivatives_update(sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y, sig_del_x_del_q, sig_del_y_del_q, max_q, min_q);

        for(int i=0; i<4; i++)
        {
            globaldata[idx].dq1[i] = max_q[i];
            globaldata[idx].dq2[i] = min_q[i];

            //cout<<"\n Q_Ders: "<<max_q[i]<<"\t"<<min_q[i];
        }

    }
    cout<<"\nGoing outta q_var_derivatives\n";
}

inline void q_var_derivatives_update(double sig_del_x_sqr, double sig_del_y_sqr, double sig_del_x_del_y, double sig_del_x_del_q[4], double sig_del_y_del_q[4], double dq1_store[4], double dq2_store[2])
{
    double det = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
    double one_by_det = 1.0/det;
    for(int iter=0; iter<4; iter++)
    {
        dq1_store[iter] = one_by_det * (sig_del_x_del_q[iter] * sig_del_y_sqr - sig_del_y_del_q[iter] * sig_del_x_del_y);
        dq2_store[iter] = one_by_det * (sig_del_y_del_q[iter] * sig_del_x_sqr - sig_del_x_del_q[iter] * sig_del_x_del_y);
    }
}

void q_var_derivatives_innerloop(Point* globaldata, int numPoints, double power, double tempdq[][2][4], double sig_del_x_del_q[4], double sig_del_y_del_q[4], double qi_tilde[4], double qk_tilde[4])
{   
    cout<<"\nInside q_var_derivatives Inner\n";

    for(int idx=0; idx<numPoints; idx++)
    {
        double x_i = globaldata[idx].x;
        double y_i = globaldata[idx].y;
        double sig_del_x_sqr = 0.0;
        double sig_del_y_sqr = 0.0;
        double sig_del_x_del_y = 0.0;

        for(int i=0; i<4; i++)
        {
            sig_del_x_del_q[i] = 0.0;
            sig_del_y_del_q[i] = 0.0;
        }

        for(int i=0; i<20; i++)
        {
            //cout<<i<<" ";
            int conn = globaldata[idx].conn[i];
            if(conn == 0) break;
            double x_k = globaldata[conn].x;
            double y_k = globaldata[conn].y;

            double delta_x = x_k - x_i;
            double delta_y = y_k - y_i;

            double dist = hypot(delta_x, delta_y);
            double weights = pow(dist, power);
            sig_del_x_sqr += (delta_x * delta_x) * weights;
            sig_del_y_sqr += (delta_y * delta_y) * weights;
            sig_del_x_del_y += (delta_x * delta_y) * weights;

            q_var_derivatives_get_sum_delq_innerloop(globaldata, idx, conn, weights, delta_x, delta_y, qi_tilde, qk_tilde, sig_del_x_del_q, sig_del_y_del_q);

            double det = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
            double one_by_det = 1.0/det;

            for(int iter =0; iter<4; iter++)
            {
                tempdq[idx][0][iter] = one_by_det * (sig_del_x_del_q[iter] * sig_del_y_sqr - sig_del_y_del_q[iter] * sig_del_x_del_y);
                tempdq[idx][1][iter] = one_by_det * (sig_del_y_del_q[iter] * sig_del_x_sqr - sig_del_x_del_q[iter] * sig_del_x_del_y);
            }

            for(int idx=0; idx<numPoints; idx++)
            {
                q_var_derivatives_update_innerloop(qi_tilde, qk_tilde, idx, tempdq);
                for(int i=0; i<4; i++)
                {    
                    globaldata[idx].dq1[i] = qi_tilde[i];
                    globaldata[idx].dq2[i] = qk_tilde[i];

                    cout<<"\n Q_Ders: "<<qi_tilde[i]<<"\t"<<qk_tilde[i];
                }
            }

        }
    }
    cout<<"\nOutta q_var_derivatives Inner\n";
}

inline void q_var_derivatives_get_sum_delq_innerloop(Point* globaldata, int idx, int conn, double weights, double delta_x, double delta_y, double qi_tilde[4], double qk_tilde[4], double sig_del_x_del_q[4], double sig_del_y_del_q[4])
{
    for(int iter=0; iter<4; iter++)
    {
        qi_tilde[iter] = globaldata[idx].q[iter] - 0.5 * (delta_x * globaldata[idx].dq1[iter] + delta_y * globaldata[idx].dq2[iter]);
        qk_tilde[iter] = globaldata[conn].q[iter] - 0.5 * (delta_x * globaldata[conn].dq1[iter] + delta_y * globaldata[conn].dq2[iter]);

        double intermediate_var = weights * (qk_tilde[iter] - qi_tilde[iter]);
        sig_del_x_del_q[iter] += delta_x * intermediate_var;
        sig_del_y_del_q[iter] += delta_y * intermediate_var;
    }
}

inline void q_var_derivatives_update_innerloop(double dq1[4], double dq2[4], int idx, double tempdq[][2][4])
{
    for(int iter=0; iter<4; iter++)
    {
        dq1[iter] = tempdq[idx][0][iter];
        dq2[iter] = tempdq[idx][1][iter];
    }
}

void printDebug(Point* globaldata, int numPoints, Config configData, int iter, double res_old[1], int rk, double sig_del_x_del_f[4], double sig_del_y_del_f[4], double main_store[62])
{
    std::ofstream fdebug("debug_output.txt", std::ios_base::app);

    fdebug<<"\nChecking the Outputs: \n";
    fdebug<<"numPoints:"<<numPoints<<"\n";

    fdebug<<"Iteration: "<<iter+1<<"\n";
    fdebug<<"Res_Old: "<<res_old[0]<<"\n";
    fdebug<<"rk: "<<rk<<"\n";
    for(int i=0; i<4; i++)
        fdebug<<"sig_del_x_del_f "<<i<<": "<<sig_del_x_del_f[i]<<"\n";
    for(int i=0; i<4; i++)
        fdebug<<"sig_del_y_del_f "<<i<<": "<<sig_del_y_del_f[i]<<"\n";
    for(int i=0; i<62; i++)
        fdebug<<"main_store "<<i<<": "<<main_store[i]<<"\n";
    fdebug.close();
}


void debug_globaldata(Point* globaldata, int numPoints, int iter, int rk)
{
    std::ofstream fdebug("debug_globaldata.txt", std::ios_base::app);
    fdebug<<"Iteration: "<<iter+1<<" And rk: "<<rk<<"\n";
    for(int i=0; i<numPoints; i++)
    {
        fdebug<<"Point:  "<<i<<"\n";
            for(int j=0; j<4; j++)
                fdebug<<globaldata[i].q[j]<<"\t";
    }
    fdebug.close();
}
