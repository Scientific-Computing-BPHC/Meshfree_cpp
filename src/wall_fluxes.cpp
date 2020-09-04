#include "wall_fluxes.hpp"
#include "limiters.hpp"
#include "quadrant_fluxes.hpp"
#include "split_fluxes.hpp"

template <class Type>
bool isNan(Type var)
{
    if(var!=var) return true;
    return false;
}

void wall_dGx_pos(Point* globaldata, int idx, double gamma, double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const, double Gxp[4])
{
	double sig_del_x_sqr = 0.0;
	double sig_del_y_sqr = 0.0;
	double sig_del_x_del_y = 0.0;

	for(int i=0; i<4; i++)
	{
		sig_del_x_del_f[i] = 0.0;
		sig_del_y_del_f[i] = 0.0;
	}

	double x_i = globaldata[idx].x;
	double y_i = globaldata[idx].y;

	double nx = globaldata[idx].nx;
	double ny = globaldata[idx].ny;

	double tx = ny;
	double ty = -nx;

	for(int i=0; i<20; i++)
	{
		int conn = globaldata[idx].xpos_conn[i];
		if(conn == 0) break;

		conn = conn - 1;

		if(idx==804) cout<<"WITH EACH I: "<<i<<" Conn is: "<<conn<<" gloabl conn: "<<globaldata[conn].x<<"\t"<<globaldata[conn].y<<endl;

		double delta_x, delta_y, delta_s_weights, delta_n_weights;
		std::tie(delta_x, delta_y, delta_s_weights, delta_n_weights, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y) = connectivity_stats(x_i, y_i, nx, ny, power, globaldata[conn].x, globaldata[conn].y, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y);

		calculate_qtile(qtilde_i, qtilde_k, globaldata, idx, conn, delta_x, delta_y, vl_const, gamma, limiter_flag, phi_i, phi_k);

		qtilde_to_primitive(result, qtilde_i, gamma);

		flux_quad_GxII(G_i, nx, ny, result[0], result[1], result[2], result[3]);

		qtilde_to_primitive(result, qtilde_k, gamma);
        flux_quad_GxII(G_k, nx, ny, result[0], result[1], result[2], result[3]);

        update_delf(sig_del_x_del_f, sig_del_y_del_f, G_k, G_i, delta_s_weights, delta_n_weights);
     }

    double det = sig_del_x_sqr * sig_del_y_sqr - sig_del_x_del_y * sig_del_x_del_y;
    double one_by_det = 1.0/det;
    for(int iter =0; iter<4; iter++)
    {
    	Gxp[iter] = (sig_del_x_del_f[iter]*sig_del_y_sqr - sig_del_y_del_f[iter]*sig_del_x_del_y)*one_by_det;
    	if(isNan(Gxp[iter]))
    	{
            cout<<"Debug details, error in Gxp inside dGxpos: "<<endl;
            cout<<"i: "<<iter<<endl;
            cout<<"idx: "<<idx<<endl;
            cout<<"flag: "<<globaldata[idx].flag_1<<endl;
            cout<<"One by det"<<one_by_det<<endl;
            cout<<"Det"<<det<<endl;
            cout<<"Update delf: "<<sig_del_x_del_f[iter]<<"\t"<<sig_del_y_del_f[iter]<<endl;
            exit(0); 
    	}

    }

	
}

void wall_dGx_neg(Point* globaldata, int idx, double gamma, double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const, double Gxn[4])
{
	double sig_del_x_sqr = 0.0;
	double sig_del_y_sqr = 0.0;
	double sig_del_x_del_y = 0.0;

	for(int i=0; i<4; i++)
	{
		sig_del_x_del_f[i] = 0.0;
		sig_del_y_del_f[i] = 0.0;
	}

	double x_i = globaldata[idx].x;
	double y_i = globaldata[idx].y;

	double nx = globaldata[idx].nx;
	double ny = globaldata[idx].ny;

	double tx = ny;
	double ty = -nx;

	for(int i=0; i<20; i++)
	{
		int conn = globaldata[idx].xneg_conn[i];
		if(conn == 0) break;

		conn = conn - 1;
		//cout<<"WITH EACH I: "<<i<<" Conn is: "<<conn<<" gloabl conn: "<<globaldata[conn].x<<"\t"<<globaldata[conn].y<<endl;

		double delta_x, delta_y, delta_s_weights, delta_n_weights;

		int deb=0;
		if(deb)
		{
			cout<<"\n \n \nsig_del_x_y_ before limiterssssssssssss: "<<sig_del_x_del_y<<endl;	
		}



		std::tie(delta_x, delta_y, delta_s_weights, delta_n_weights, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y) = connectivity_stats(x_i, y_i, nx, ny, power, globaldata[conn].x, globaldata[conn].y, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y);



     	
        if(deb)
            {	
                cout<<"_del_x__: "<<delta_x<<endl;
                cout<<"_del_y__: "<<delta_y<<endl;
                cout<<"sig_del_x_y_: "<<sig_del_x_del_y<<endl;
                cout<<"weights: "<<delta_s_weights<<"\t"<<delta_n_weights<<endl;
                cout<<"gloabl conn: "<<globaldata[conn].x<<"\t"<<globaldata[conn].y<<endl;
                cout<<"\nlet's see what conn is and then let's see what it's co-ords are"<<endl;
                cout<<conn<<endl;
                cout<<globaldata[conn].x<<"\t"<<globaldata[conn].y<<endl;
                cout<<"Wokay"<<endl;
            }
        //exit(0);


		calculate_qtile(qtilde_i, qtilde_k, globaldata, idx, conn, delta_x, delta_y, vl_const, gamma, limiter_flag, phi_i, phi_k);

		qtilde_to_primitive(result, qtilde_i, gamma);

		flux_quad_GxI(G_i, nx, ny, result[0], result[1], result[2], result[3]);

		qtilde_to_primitive(result, qtilde_k, gamma);
        flux_quad_GxI(G_k, nx, ny, result[0], result[1], result[2], result[3]);

        if(deb)
            {
                cout<<"sig_del_x_sqr_: "<<sig_del_x_sqr<<endl;
                cout<<"sig_del_y_sqr_: "<<sig_del_y_sqr<<endl;
                cout<<"sig_del_x_y_: "<<sig_del_x_del_y<<endl;
                double tmp = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
                cout<<"My man inside wall dGx:_1 "<<tmp<<endl;
                if(!tmp) cout<<"Nbh is "<<i<<" "<<" Point is "<<idx<<" Check karo "<<endl;
            }

        if(deb)
        {
        	cout<<"\n \n \n CALLS TO UPDATE DEL_F"<<endl;
        	cout<<"Params: "<<delta_s_weights<<delta_n_weights;
        	cout<<"The Sigs and Gs, and then EXITING after 20 runs";
        	for(int a=0; a<4; a++)
        	{
        		cout<<"\n"<<sig_del_x_del_f[a]<<"    "<<sig_del_y_del_f[a]<<"    "<<G_i[i]<<"    "<<G_k[i];
        	}
        }

        update_delf(sig_del_x_del_f, sig_del_y_del_f, G_k, G_i, delta_s_weights, delta_n_weights);

         if(deb)
        {
            cout<<"sig_del_x_sqr_: "<<sig_del_x_sqr<<endl;
            cout<<"sig_del_y_sqr_: "<<sig_del_y_sqr<<endl;
            cout<<"sig_del_x_y_: "<<sig_del_x_del_y<<endl;
            double tmp = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
            cout<<"My man inside wall dGx:_2 "<<tmp<<endl;
            if(!tmp) cout<<"Nbh is "<<i<<" "<<" Point is "<<idx<<" Check karo "<<endl;
        }
     }
    //exit(0);

    double det = sig_del_x_sqr * sig_del_y_sqr - sig_del_x_del_y * sig_del_x_del_y;
    double one_by_det = 1.0/det;
    for(int iter =0; iter<4; iter++)
    {
    	Gxn[iter] = (sig_del_x_del_f[iter]*sig_del_y_sqr - sig_del_y_del_f[iter]*sig_del_x_del_y)*one_by_det;
    	if(isNan(Gxn[iter]))
    	{
            cout<<"Debug details, error in Gxn inside dGxneg: "<<endl;
            cout<<"i: "<<iter<<endl;
            cout<<"idx: "<<idx<<endl;
            cout<<"flag: "<<globaldata[idx].flag_1<<endl;
            cout<<"One by det"<<one_by_det<<endl;
            cout<<"Det"<<det<<endl;
            cout<<"Update delf: "<<sig_del_x_del_f[iter]<<"\t"<<sig_del_y_del_f[iter]<<endl;
            exit(0); 
    	}
    }

	
}

void wall_dGy_neg(Point* globaldata, int idx, double gamma, double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double power, int limiter_flag, double vl_const, double Gyn[4])
{
	double sig_del_x_sqr = 0.0;
	double sig_del_y_sqr = 0.0;
	double sig_del_x_del_y = 0.0;

	for(int i=0; i<4; i++)
	{
		sig_del_x_del_f[i] = 0.0;
		sig_del_y_del_f[i] = 0.0;
	}

	double x_i = globaldata[idx].x;
	double y_i = globaldata[idx].y;

	double nx = globaldata[idx].nx;
	double ny = globaldata[idx].ny;

	double tx = ny;
	double ty = -nx;

	for(int i=0; i<20; i++)
	{
		int conn = globaldata[idx].yneg_conn[i];
		if(conn == 0) break;

		conn = conn - 1;

		double delta_x, delta_y, delta_s_weights, delta_n_weights;
		std::tie(delta_x, delta_y, delta_s_weights, delta_n_weights, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y) = connectivity_stats(x_i, y_i, nx, ny, power, globaldata[conn].x, globaldata[conn].y, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y);

		calculate_qtile(qtilde_i, qtilde_k, globaldata, idx, conn, delta_x, delta_y, vl_const, gamma, limiter_flag, phi_i, phi_k);

		qtilde_to_primitive(result, qtilde_i, gamma);

		flux_Gyn(G_i, nx, ny, result[0], result[1], result[2], result[3]);

		qtilde_to_primitive(result, qtilde_k, gamma);

        flux_Gyn(G_k, nx, ny, result[0], result[1], result[2], result[3]);

        int deb=0;
        if(deb)
            {
                cout<<"sig_del_x_sqr_: "<<sig_del_x_sqr<<endl;
                cout<<"sig_del_y_sqr_: "<<sig_del_y_sqr<<endl;
                cout<<"sig_del_x_y_: "<<sig_del_x_del_y<<endl;
                double tmp = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
                cout<<"My man inside wall Gdy:_1 "<<tmp<<endl;
                if(!tmp) cout<<"Nbh is "<<i<<" "<<" Point is "<<idx<<" Check karo "<<endl;
            }

        update_delf(sig_del_x_del_f, sig_del_y_del_f, G_k, G_i, delta_s_weights, delta_n_weights);

        if(deb)
            {
                cout<<"sig_del_x_sqr_: "<<sig_del_x_sqr<<endl;
                cout<<"sig_del_y_sqr_: "<<sig_del_y_sqr<<endl;
                cout<<"sig_del_x_y_: "<<sig_del_x_del_y<<endl;
                double tmp = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
                cout<<"My man inside wall Gdy:_2 "<<tmp<<endl;
                if(!tmp) cout<<"Nbh is "<<i<<" "<<" Point is "<<idx<<" Check karo "<<endl;
            }
     }

    double det = sig_del_x_sqr * sig_del_y_sqr - sig_del_x_del_y * sig_del_x_del_y;
    double one_by_det = 1.0/det;
    for(int iter =0; iter<4; iter++)
    {
    	Gyn[iter] = (sig_del_y_del_f[iter]*sig_del_x_sqr - sig_del_x_del_f[iter]*sig_del_x_del_y)*one_by_det;
    	if(isNan(Gyn[iter]))
    	{
            cout<<"Debug details, error in Gyn inside dGyneg: "<<endl;
            cout<<"i: "<<iter<<endl;
            cout<<"idx: "<<idx<<endl;
            cout<<"flag: "<<globaldata[idx].flag_1<<endl;
            cout<<"One by det: "<<one_by_det<<endl;
            cout<<"Det: "<<det<<endl;
            cout<<"Update delf: "<<sig_del_x_del_f[iter]<<"\t"<<sig_del_y_del_f[iter]<<endl;
            exit(0); 
    	}
    }

	
}