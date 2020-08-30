#include "point.hpp"
#include "state_update.hpp"

void func_delta(Point* globaldata, int numPoints, double cfl)
{
	{
	for(int idx=0; idx<numPoints; idx++)
	{
		double min_delt = 1.0;
		for(int i=0; i<20; i++)
		{
			int conn = globaldata[idx].conn[i];
			if (conn == 0) break;

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
}