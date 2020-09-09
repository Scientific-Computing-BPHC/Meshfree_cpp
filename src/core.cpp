#include "core.hpp"
#include "point.hpp"
#include "state_update.hpp"
#include "flux_residual.hpp"

inline void q_var_derivatives_update(double sig_del_x_sqr, double sig_del_y_sqr, double sig_del_x_del_y, double sig_del_x_del_q[4], double sig_del_y_del_q[4], double dq1_store[4], double dq2_store[2]);
inline void q_var_derivatives_get_sum_delq_innerloop(Point* globaldata, int idx, int conn, double weights, double delta_x, double delta_y, double qi_tilde[4], double qk_tilde[4], double sig_del_x_del_q[4], double sig_del_y_del_q[4]);
inline void q_var_derivatives_update_innerloop(double dq1[4], double dq2[4], int idx, double tempdq[][2][4]);

template <class Type>
bool isNan(Type var)
{
    if(var!=var) return true;
    return false;
}



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

        itm = itm -1; // to account for indexing

    	//cout<< "\n Unbroken \n";
    	double itmx = globaldata[itm].x;
    	double itmy = globaldata[itm].y;

    	double delta_x = itmx - currx;
    	double delta_y = itmy - curry;

    	double delta_s = delta_x*tx + delta_y*ty;
    	double delta_n = delta_x*nx + delta_y*ny;

        itm = itm + 1; // to reaccount for indexing when we add the point below xpos_conn[xpos_nbhs] = itm;

    	if(delta_s <= 0.0)
    	{
    		
    		xpos_conn[xpos_nbhs] = itm;
            xpos_nbhs+=1;
    	}

    	if(delta_s >= 0.0)
    	{
    		
    		xneg_conn[xneg_nbhs] = itm;
            xneg_nbhs+=1;
    	}

    	if(flag==1)
    	{
    		if(delta_n<=0.0)
    		{
    			
    			ypos_conn[ypos_nbhs] = itm;
                ypos_nbhs+=1;
    		}

    		if(delta_n>=0.0)
    		{
    			
    			yneg_conn[yneg_nbhs] = itm;
                yneg_nbhs+=1;
    		}
    	}

    	else if (flag==0)
    	{
    		
    		yneg_conn[yneg_nbhs] = itm;
            yneg_nbhs+=1;
    	}

    	else if (flag==2)
    	{
    		
    		ypos_conn[ypos_nbhs] = itm;
            ypos_nbhs+=1;
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

    int rks = configData.core.rks;
    int euler = configData.core.euler;

    func_delta(globaldata, numPoints, cfl);

    //debug_globaldata(globaldata, numPoints, iter, 0);

    // double phi_i[4], phi_k[4], G_i[4], G_k[4], result[4], qtilde_i[4], qtilde_k[4];
    // double Gxp[4], Gxn[4], Gyp[4], Gyn[4], sig_del_x_del_f[4], sig_del_y_del_f[4];

    // for(int i=0; i<4; i++)
    //     phi_i[i] = main_store[i];
    // for(int i=4; i<8; i++)
    //     phi_k[i-4] = main_store[i];
    // for(int i=8; i<12; i++)
    //     G_i[i-8] = main_store[i];
    // for(int i=12; i<16; i++)
    //     G_k[i-12] = main_store[i];
    // for(int i=16; i<20; i++)
    //     result[i-16] = main_store[i];
    // for(int i=20; i<24; i++)
    //     qtilde_i[i-20] = main_store[i];
    // for(int i=24; i<28; i++)
    //     qtilde_k[i-24] = main_store[i];
    // for(int i=28; i<32; i++)
    //     Gxp[i-28] = main_store[i];
    // for(int i=32; i<36; i++)
    //     Gxn[i-32] = main_store[i];
    // for(int i=36; i<40; i++)
    //     Gyp[i-36] = main_store[i];
    // for(int i=40; i<44; i++)
    //     Gyn[i-40] = main_store[i];
    // for(int i=44; i<48; i++)
    //     sig_del_x_del_f[i-44] =main_store[i];
    // for(int i=48; i<52; i++)
    //     sig_del_y_del_f[i-48] = main_store[i];

    double *phi_i, *phi_k, *G_i, *G_k, *result, *qtilde_i, *qtilde_k;
    double *Gxp, *Gxn, *Gyp, *Gyn, *sig_del_x_del_f, *sig_del_y_del_f;

    int i = 0;
    phi_i = &(main_store[i]);
    phi_k = &(main_store[i+4]);
    G_i = &(main_store[i+8]);
    G_k = &(main_store[i+12]);
    result = &(main_store[i+16]);
    qtilde_i = &(main_store[i+20]);
    qtilde_k = &(main_store[i+24]);
    Gxp = &(main_store[i+28]);
    Gxn = &(main_store[i+32]);
    Gyp = &(main_store[i+36]);
    Gyn = &(main_store[i+40]);
    sig_del_x_del_f = &(main_store[i+44]);
    sig_del_y_del_f = &(main_store[i+48]);

    //cout<<"Intial Check: "<<main_store[23]<<endl;

    cout<<"\nIteration Number: (iter+1)  "<<iter+1<<endl;

    for(int rk=0; rk<rks; rk++)
    {
        
        cout<<"\n rk: (rk+1)  "<<rk+1<<"\n";

        int deb = 1;
        if(deb)
        {
            int idx = 0;
            {
                for(int i=0; i<4; i++)
                {
                    cout<<"Prim: "<<i<<" "<<globaldata[idx].prim[i]<<endl;
                }
            }
        }

        q_variables(globaldata, numPoints, result);

        if(deb)
        {
            int idx = 0;
            {
                for(int i=0; i<4; i++)
                {
                    cout<<"Result: "<<i<<" "<<result[i]<<", dq1: "<<i<<" "<<globaldata[idx].dq1[i]<<", dq2: "<<i<<" "<<globaldata[idx].dq2[i]<<endl;
                }
            }
        }

        //debug_globaldata(globaldata, numPoints, iter, rk);

        q_var_derivatives(globaldata, numPoints, power, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k);

        //debug_main_store_3(main_store);

        for(int inner_iters=0; inner_iters<3; inner_iters++)
        {
            q_var_derivatives_innerloop(globaldata, numPoints, power, tempdq, sig_del_x_del_f, sig_del_y_del_f, qtilde_i, qtilde_k);
        }

        //debug_main_store_3(main_store);

        cout<<"\nCalculating Flux Residual\n";

        //debug_globaldata(globaldata, numPoints, iter, rk, main_store);


        // cout<<"\n \n \n";
        // cout<<"This is iter: "<<iter+1<<endl;
        // cout<<"This is rk (rk+1): "<<rk<<endl;
        // cout<<"Let's see the qtildes: "<<endl;
        // for(int m=0; m<4; m++)
        //     cout<<qtilde_i[m]<<"  "<<qtilde_k[m]<<"\n";
        // cout<<"\n \n \n";

        //debug_Gs_and_qtildes(iter, rk, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k, result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, main_store);

        cal_flux_residual(globaldata, numPoints, configData, Gxp, Gxn, Gyp, Gyn, phi_i, phi_k, G_i, G_k,
            result, qtilde_i, qtilde_k, sig_del_x_del_f, sig_del_y_del_f, main_store);

        // if(rk==1)
        // {
        //     cout<<"Exiting for rk==1, inside core";
        //     exit(0);
        // }



        //printDebug(globaldata, numPoints, configData, iter, res_old, rk, sig_del_x_del_f, sig_del_y_del_f, main_store);

        cout<<"\nDone Calculating Flux Residual\n";

        state_update(globaldata, numPoints, configData, iter, res_old, rk, sig_del_x_del_f, sig_del_y_del_f, main_store, euler);

        //debug_globaldata(globaldata, numPoints, iter, rk, main_store);
        //printDebug(globaldata, numPoints, configData, iter, res_old, rk, sig_del_x_del_f, sig_del_y_del_f, main_store);

        // if(iter == 0 && rk == 3)
        // {
        //     cout<<"\n Stopping inside fpi solver, 4th rk and 1st iteration done \n";
        //     exit(0);

        // }

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
        double beta = 0.5 * (rho/pr);

        if(isNan(rho) || isNan(u1) || isNan (u2) || isNan(pr))
        {
            cout<<"Ah we have a Nan inside q_variables here";
            exit(0);
        }

        if(beta!=beta)
        {
            cout<<"\nNAN encountered at Beta\n";
            exit(0);
        }
        double two_times_beta = 2.0 * beta;
        q_result[0] = log(rho) + log(beta) * 2.5 - (beta * ((u1 * u1) + (u2 * u2)));
        if(q_result[0]!=q_result[0])
        {
            cout<<"\nNAN encountered at q_result[0]\n";
            cout<<"Let's find out who the culprit is: "<<endl;
            cout<<"Rho: "<<rho<<endl;
            cout<<"Beta: "<<beta<<endl;
            cout<<"u1: "<<u1<<endl;
            cout<<"u2: "<<u2<<endl;
            cout<<"The Point is: "<<idx<<endl;
            exit(0);
        }
        q_result[1] = (two_times_beta * u1);
        q_result[2] = (two_times_beta * u2);
        q_result[3] = -two_times_beta;
        for(int i=0; i<4; i++)
        {
            globaldata[idx].q[i] = q_result[i];
            //cout<<"q_result: "<<i<<" "<<q_result[i]<<endl;
        }
        for(int i=0; i<4; i++)
        {
            if(q_result[i]!=q_result[i])
            {
                cout<<"\nNAN encountered at q_result_"<<i<<"\n";
                exit(0);
            }
        }
    }
    cout<<"\nGoing outta q_variables\n";
}

void q_var_derivatives(Point* globaldata, int numPoints, double power, double sig_del_x_del_q[4], double sig_del_y_del_q[4], double max_q[4], double min_q[4])
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

        for(int i=0; i<4; i++)
        {
            max_q[i] = globaldata[idx].q[i];
            min_q[i] = globaldata[idx].q[i];
        }

        for(int i=0; i<20; i++)
        {
            //cout<<"\n Count "<<i<<endl;
            int conn = globaldata[idx].conn[i];
            if(idx == 48737) cout<<"Conn: "<<conn<<endl;
            if(conn == 0) 
            {
                //cout<<"BROKENNNNNN Lesse i"<<i<<endl;
                //cout<<"Just to check if exiting";
                //exit(0);
                break;
            }
            if(conn!=conn)
            {
                    cout<<"\nNAN encountered at conn\n";
                    exit(0);
            }

            conn = conn - 1; // To account for the indexing difference

            double x_k = globaldata[conn].x;
            double y_k = globaldata[conn].y;

            double delta_x = x_k - x_i;
            double delta_y = y_k - y_i;

            //cout<<"Deltaaa x: "<<delta_x<<"Deltaa y:"<<delta_y<<endl;

            double dist = hypot(delta_x, delta_y);
            double weights = pow(dist, power);
            sig_del_x_sqr += ((delta_x * delta_x) * weights);
            sig_del_y_sqr += ((delta_y * delta_y) * weights);
            sig_del_x_del_y += ((delta_x * delta_y) * weights);

            //cout<<"\nLesseee these values:\n"<<endl;
            //cout<<sig_del_x_sqr<<"\t"<<sig_del_y_sqr<<"\t"<<sig_del_x_del_y<<endl;

            if(sig_del_x_sqr!=sig_del_x_sqr)
                {
                    cout<<"\nNAN encountered at sig_del_x_sqr\n";
                    exit(0);
                }
            if(sig_del_y_sqr!=sig_del_y_sqr)
                {
                    cout<<"\nNAN encountered at sig_del_y_sqr\n";
                    exit(0);
                }
            if(sig_del_x_del_y!=sig_del_x_del_y)
                {
                    cout<<"\nNAN encountered at sig_del_x_del_y\n";
                    exit(0);
                }
            for(int iter=0; iter<4; iter++)
            {
                double intermediate_var = weights * (globaldata[conn].q[iter] - globaldata[idx].q[iter]);
                sig_del_x_del_q[iter] += (delta_x * intermediate_var);
                sig_del_y_del_q[iter] += (delta_y * intermediate_var);
                if(intermediate_var!=intermediate_var)
                {
                    cout<<"\nNAN encountered at intermediate_var\n";
                    exit(0);
                }

            }

            for(int j=0; j<4; j++)
            {
                //cout<<"Yoo max: "<<max_q[j]<<endl;
                if (max_q[j] < globaldata[conn].q[j])
                {
                    max_q[j] = globaldata[conn].q[j];
                    //cout<<"\t"<<max_q[j]<<"\t";
                    if(globaldata[conn].q[j]!=globaldata[conn].q[j])
                    {
                        cout<<"\nNan\n";
                        exit(0);
                    }
                }
                if(min_q[j] > globaldata[conn].q[j])
                {
                    min_q[j] = globaldata[conn].q[j];
                    //cout<<"\t"<<min_q[j]<<"\t";
                    if(globaldata[conn].q[j]!=globaldata[conn].q[j])
                    {
                        cout<<"\nNan\n";
                        exit(0);
                    }
                }
            }
        }



        for(int i=0; i<4; i++)
        {
            globaldata[idx].max_q[i] = max_q[i];
            //cout<<"\n";
            //cout<<"\t"<<max_q[j]<<"\t";
            //cout<<"\n And \n";
            globaldata[idx].min_q[i] = min_q[i];
            //cout<<"\t"<<min_q[j]<<"\t";
        }

        for(int i=0; i<4; i++)
        {    
            if(max_q[i]!=max_q[i])
            {
                cout<<"\nNAN encountered at max_q[i]\n";
                cout<<"\t"<<max_q[i];
                exit(0);
            }

            if(min_q[i]!=min_q[i])
            {
                cout<<"\nNAN encountered at min_q[i]\n";
                exit(0);
            }
        }

        

        //cout<<"Lesse what we're passing"<<"\n"<<endl;
        //cout<<sig_del_x_sqr<<"\t"<<sig_del_y_sqr<<"\t"<<sig_del_x_del_y<<endl;

        //cout<<"\nLesse what would happen to det"<<"\n"<<endl;
        //cout<<sig_del_x_sqr*sig_del_y_sqr<<"\t"<<sig_del_x_del_y*sig_del_x_del_y<<endl;

        if(idx==48737) cout<<"Before update min_q[0]: "<<min_q[0]<<endl;

        q_var_derivatives_update(sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y, sig_del_x_del_q, sig_del_y_del_q, max_q, min_q);

        if(idx==48737) cout<<"After update min_q[0]: "<<min_q[0]<<endl;
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
    // cout<<"\nThis is sqr:"<<(sig_del_x_sqr * sig_del_y_sqr);
    // cout<<"\nSub from sqr: "<<(sig_del_x_del_y * sig_del_x_del_y);
    // cout<<"\nThis is det: "<<det<<endl;
    double one_by_det = 1.0/det;

    if(one_by_det!=one_by_det || det!=det)
    {
        cout<<"\nNAN encountered at one_by_det or det\n";
        exit(0);
    }

    for(int iter=0; iter<4; iter++)
    {
        if(sig_del_x_sqr != sig_del_x_sqr)
        {
            cout<<"sig_del_x_sqr_ Why god, why?"<<endl;
            exit(0);
        }

        if(sig_del_y_sqr != sig_del_y_sqr)
        {
            cout<<"sig_del_y_sqr_ Why god, why?"<<endl;
            exit(0);
        }


        if(sig_del_x_del_q[iter] != sig_del_x_del_q[iter])
        {
            cout<<"sig_del_x_q_ Why god, why?"<<endl;
            exit(0);
        }


        if(sig_del_y_del_q[iter] != sig_del_y_del_q[iter])
        {
            cout<<"sig_del_y_q_ Why god, why?"<<endl;
            exit(0);
        }


        if(sig_del_x_del_y != sig_del_x_del_y)
        {
            cout<<"sig_del_x_y_ Why god, why?"<<endl;
            exit(0);
        }


    }

    for(int iter=0; iter<4; iter++)
    {
        dq1_store[iter] = one_by_det * (sig_del_x_del_q[iter] * sig_del_y_sqr - sig_del_y_del_q[iter] * sig_del_x_del_y);
        dq2_store[iter] = one_by_det * (sig_del_y_del_q[iter] * sig_del_x_sqr - sig_del_x_del_q[iter] * sig_del_x_del_y);
    }


    for(int iter=0; iter<4; iter++)
    {
        if(dq1_store[iter] != dq1_store[iter])
        {
            cout<<"Why god, why dq1? + "<<iter<<endl;
            cout<<"One by det: "<<one_by_det<<endl;
            exit(0);
        }

        if(dq2_store[iter] != dq2_store[iter])
        {
            cout<<"Why god, why? dq2 + "<<iter<<endl;
            cout<<"One by det: "<<one_by_det<<endl;
            exit(0);
        }
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
            //cout<<i<<" "<<endl;
            int conn = globaldata[idx].conn[i];
            if(conn == 0) break;

            conn = conn - 1;

            double x_k = globaldata[conn].x;
            double y_k = globaldata[conn].y;

            if(idx==48737) cout<<"conn (+1): "<<conn+1<<" y_k: "<<y_k<<endl;

            double delta_x = x_k - x_i;
            double delta_y = y_k - y_i;

            double dist = hypot(delta_x, delta_y);
            double weights = pow(dist, power);
            sig_del_x_sqr += ((delta_x * delta_x) * weights);
            sig_del_y_sqr += ((delta_y * delta_y) * weights);
            sig_del_x_del_y += ((delta_x * delta_y) * weights);


            if(isNan(sig_del_x_sqr) || isNan(sig_del_y_sqr) || isNan(sig_del_x_del_y))
            {
                cout<<"Ah we have a Nan inside q_var_derivatives_innerloop here";
                exit(0);
            }  

            int deb = 0;
            if(deb)
            {
                cout<<"sig_del_x_sqr_: "<<sig_del_x_sqr<<endl;
                cout<<"sig_del_y_sqr_: "<<sig_del_y_sqr<<endl;
                cout<<"sig_del_x_y_: "<<sig_del_x_del_y<<endl;
                double tmp = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
                cout<<"My man:_1 "<<tmp<<endl;
                if(!tmp) cout<<"Nbh is "<<i<<" "<<" Point is "<<idx<<" Check karo "<<endl;
            }

            if(idx==48737) cout<<"\n Sig dely delq BEFOREEEE[0]: "<<sig_del_y_del_q[0]<<endl;

            q_var_derivatives_get_sum_delq_innerloop(globaldata, idx, conn, weights, delta_x, delta_y, qi_tilde, qk_tilde, sig_del_x_del_q, sig_del_y_del_q);
        }

        int deb = 0;
        if(deb)
        {
            cout<<"sig_del_x_sqr_: "<<sig_del_x_sqr<<endl;
            cout<<"sig_del_y_sqr_: "<<sig_del_y_sqr<<endl;
            cout<<"sig_del_x_y_: "<<sig_del_x_del_y<<endl;
            cout<<"My man_2: "<<(sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
        }

        if(idx == 48737)
        {
            cout<<"\n Sig x sqr: "<<sig_del_x_sqr<<endl;
            cout<<"\n Sig y sqr: "<<sig_del_y_sqr<<endl;
            cout<<"\n Sig delx dely: "<<sig_del_x_del_y<<endl;
            cout<<"\n Sig delx delq[0]: "<<sig_del_x_del_q[0]<<endl;
            cout<<"\n Sig dely delq[0]: "<<sig_del_y_del_q[0]<<endl;
        }

        double det = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
        double one_by_det = 1.0/det;

        if(isNan(det) || isNan(one_by_det))
        {
            cout<<"Ah we have a Nan inside q_var_derivatives_innerloop one by det or det here";
            exit(0);
        } 

        for(int iter =0; iter<4; iter++)
        {
            
            if(isNan(sig_del_x_del_q[iter]) || isNan(sig_del_y_del_q[iter]))
            {
                cout<<"Ah we have a Nan inside q_var_derivatives_innerloop here sig_del";
                exit(0);
            } 

            tempdq[idx][0][iter] = one_by_det * (sig_del_x_del_q[iter] * sig_del_y_sqr - sig_del_y_del_q[iter] * sig_del_x_del_y);
            tempdq[idx][1][iter] = one_by_det * (sig_del_y_del_q[iter] * sig_del_x_sqr - sig_del_x_del_q[iter] * sig_del_x_del_y);

            if(isNan(tempdq[idx][0][iter]) || isNan(tempdq[idx][1][iter]))
            {
                cout<<"Ah we have a Nan inside q_var_derivatives_innerloop tempdq here\n";
                if(isNan(sig_del_x_del_q[iter])) cout<<"\n because of sig x q";
                if(isNan(sig_del_y_del_q[iter])) cout<<"\n because of sig y q";
                if(isNan(sig_del_x_sqr)) cout<<"\n because of sig x sq";
                if(isNan(sig_del_y_sqr)) cout<<"\n because of sig y sq";
                if(isNan(one_by_det)) cout<<"\n because of one_by_det";
                cout<<(sig_del_x_del_q[iter] * sig_del_y_sqr - sig_del_y_del_q[iter] * sig_del_x_del_y)<<endl;
                cout<<"one by det"<<one_by_det<<endl;
                cout<<"Ah damn\n";
                cout <<"Hi"<<endl;
                exit(0);
            } 
        }
    }

    cout<<"\nTempdq[1][0] for idx 48737: "<<tempdq[48737][1][0];


    for(int k=0; k<numPoints; k++)
    {
        q_var_derivatives_update_innerloop(qi_tilde, qk_tilde, k, tempdq);
        for(int j=0; j<4; j++)
        {    
            globaldata[k].dq1[j] = qi_tilde[j];
            globaldata[k].dq2[j] = qk_tilde[j];

            //cout<<"\n Q_Ders: "<<qi_tilde[i]<<"\t"<<qk_tilde[i];
            if(isNan(qi_tilde[j]) || isNan(qk_tilde[j]))
            {
                cout<<"Ah we have a Nan inside q_var_derivatives_innerloop here, at the end";
                exit(0);
            } 
        }

        if(k == 48737)
        {
            cout<<"\n DEBUG: qtilde_i[0]: "<<qi_tilde[0]<<endl;
            cout<<"\n DEBUG: qtilde_k[0]: "<<qk_tilde[0]<<endl;
        }
    }

    cout<<"\nOutta q_var_derivatives Inner\n";
}

inline void q_var_derivatives_get_sum_delq_innerloop(Point* globaldata, int idx, int conn, double weights, double delta_x, double delta_y, double qi_tilde[4], double qk_tilde[4], double sig_del_x_del_q[4], double sig_del_y_del_q[4])
{
    for(int iter=0; iter<4; iter++)
    {
        if(isNan(qi_tilde[iter]) || isNan(qk_tilde[iter]))
        {
            cout<<"Ah we have a Nan inside q_var_derivatives_get_suminnerloop before only here";
            exit(0);
        } 

        qi_tilde[iter] = globaldata[idx].q[iter] - 0.5 * (delta_x * globaldata[idx].dq1[iter] + delta_y * globaldata[idx].dq2[iter]);
        qk_tilde[iter] = globaldata[conn].q[iter] - 0.5 * (delta_x * globaldata[conn].dq1[iter] + delta_y * globaldata[conn].dq2[iter]);

        if(isNan(qi_tilde[iter]) || isNan(qk_tilde[iter]))
        {
            cout<<"Ah we have a Nan inside q_var_derivatives_get_suminnerloop here";
            exit(0);
        }

        if(idx == 48737)
        {   
            cout<<"\n Conn: (+1) "<<conn+1<<endl;
            cout<<"Weights: "<<weights<<endl;
            cout<<"Delta_y: "<<delta_y<<endl;
            cout<<"q_i[0]: "<<qi_tilde[0]<<endl;
            cout<<"q_k[0]: "<<qk_tilde[0]<<endl;
        } 

        double intermediate_var = weights * (qk_tilde[iter] - qi_tilde[iter]);
        sig_del_x_del_q[iter] += (delta_x * intermediate_var);
        sig_del_y_del_q[iter] += (delta_y * intermediate_var);
        if(idx == 48737 && iter == 0)
        {
            cout<<"Check intermediate sig y q [0]: "<<sig_del_y_del_q[0]<<endl;
        }
        if(isNan (sig_del_x_del_q[iter]) || isNan(sig_del_y_del_q[iter]))
        {
            cout<<"Ah we have a Nan inside q_var_derivatives_get_suminnerloop here at the end";
            exit(0);
        } 
    }
}

inline void q_var_derivatives_update_innerloop(double dq1[4], double dq2[4], int idx, double tempdq[][2][4])
{
    for(int iter=0; iter<4; iter++)
    {
        dq1[iter] = tempdq[idx][0][iter];
        dq2[iter] = tempdq[idx][1][iter];

        if(isNan(dq1[iter]) || isNan(dq2[iter]))
        {
            cout<<"Ah we have a Nan inside q_var_derivatives_update_innerloop dq";
            exit(0);
        } 
    }

}

void printDebug(Point* globaldata, int numPoints, Config configData, int iter, double res_old[1], int rk, double sig_del_x_del_f[4], double sig_del_y_del_f[4], double main_store[62])
{
    std::ofstream fdebug("debug_output.txt", std::ios_base::app);

    fdebug<<"\nChecking the Outputs: \n";
    fdebug<<"numPoints:"<<numPoints<<"\n";

    fdebug<<"Iteration: "<<iter+1<<"\n";
    fdebug<<"Res_Old: "<<res_old[0]<<"\n";
    fdebug<<"rk (rk+1): "<<rk+1<<"\n";
    for(int i=0; i<4; i++)
        fdebug<<"sig_del_x_del_f "<<i<<": "<<sig_del_x_del_f[i]<<"\n";
    for(int i=0; i<4; i++)
        fdebug<<"sig_del_y_del_f "<<i<<": "<<sig_del_y_del_f[i]<<"\n";
    for(int i=0; i<62; i++)
        fdebug<<"main_store "<<i<<": "<<main_store[i]<<"\n";
    fdebug.close();
}


void debug_globaldata(Point* globaldata, int numPoints, int iter, int rk, double main_store[62])
{
    std::ofstream fdebug("debug_globaldata_dq.txt", std::ios_base::app);
    fdebug<<"Iteration: "<<iter+1<<" And rk: "<<rk<<"\n";
    for(int i=0; i<numPoints; i++)
    {
        fdebug<<"\nPoint:  "<<i<<"  ";
            for(int j=0; j<4; j++)
                fdebug<<std::setprecision(17)<<globaldata[i].dq1[j]<<" and "<<globaldata[i].dq2[j]<<" ,  ";
        //fdebug<<globaldata[i].x<<", ";
    }
    fdebug.close();
}

void debug_main_store_3(double main_store[62])
{
    std::ofstream fdebug("debug_main_store_3.txt", std::ios_base::app);
    fdebug<<"\n--------------------------------------------------------\n";
    for(int i=0; i<62; i++)
        fdebug<<"main_store "<<i<<": "<<main_store[i]<<"\n";
    fdebug.close();
}

void debug_Gs_and_qtildes(int iter, int rk, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double main_store[62])
{
    std::ofstream fdebug("debug_Gs.txt", std::ios_base::app);

    fdebug<<"\nChecking the Outputs: \n";

    fdebug<<"Iteration: "<<iter+1<<"\n";
    fdebug<<"rk (rk+1): "<<rk+1<<"\n";
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
    fdebug<<"\nPhi_i: ";
    for(int i=0; i<4; i++)
        fdebug<<phi_i[i]<<"  ";
    fdebug<<"\nPhi_k: ";
    for(int i=0; i<4; i++)
        fdebug<<phi_k[i]<<"  ";
    fdebug<<"\nG_i: ";
    for(int i=0; i<4; i++)
        fdebug<<G_i[i]<<"  ";
    fdebug<<"\nG_k: ";
    for(int i=0; i<4; i++)
        fdebug<<G_k[i]<<"  ";
    fdebug<<"\nResult: ";
    for(int i=0; i<4; i++)
        fdebug<<result[i]<<"  ";
    fdebug<<"\nqtilde_i: ";
    for(int i=0; i<4; i++)
        fdebug<<qtilde_i[i]<<"  ";
    fdebug<<"\nqtilde_k: ";
    for(int i=0; i<4; i++)
        fdebug<<qtilde_k[i]<<"  ";
    fdebug<<"\nsig_del_x_del_f: ";
    for(int i=0; i<4; i++)
        fdebug<<sig_del_x_del_f[i]<<"  ";
    fdebug<<"\nsig_del_y_del_f: ";
    for(int i=0; i<4; i++)
        fdebug<<sig_del_y_del_f[i]<<"  ";
    for(int i=0; i<62; i++)
        fdebug<<"\nmain_store "<<i<<": "<<main_store[i]<<"\n";
    fdebug.close();
}


