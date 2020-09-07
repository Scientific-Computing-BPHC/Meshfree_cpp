#ifndef CORE_HPP
#define CORE_HPP

#include<stdio.h>
#include<iostream>
#include<math.h>
#include <unistd.h> 
#include <armadillo>
#include <fstream>
#include<map>
#include<string>
#include<iomanip>
#include<limits>
#include<string>
#include<regex>
#include<sstream>
#include<tuple>

#define ARMA_DONT_PRINT_ERRORS
using namespace arma;

typedef std::map<char, float> char_float_map;
typedef std::map<char, double> char_double_map;
typedef std::vector<std::string> vec_str;
typedef std::vector<double> vec_doub;
typedef std::vector<long double> vec_ldoub;
typedef std::tuple <double, double> xy_tuple;

struct Point
{
	int localID;
	double x, y;
	int left, right;
	int flag_1, flag_2; // Int8 in the Julia code
	double short_distance;
	int nbhs;
	int conn[20];
	double nx, ny;
	// Size 4 (Pressure, vx, vy, density) x numberpts
	double prim[4];
	double flux_res[4];
	double q[4];
	// Size 2(x,y) 4(Pressure, vx, vy, density) numberpts
	double dq1[4];
	double dq2[4];
	double entropy;
	int xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs;
	int xpos_conn[20];
	int xneg_conn[20];
	int ypos_conn[20];
	int yneg_conn[20];
	double delta;
	double max_q[4];
	double min_q[4];
	double prim_old[4];

	//Point Constructor

	Point() {}

};

struct Config
{
	struct Core
	{
		int points;
		double cfl;
		int max_iters;
		double mach;
		double aoa;
		double power;
		int limiter_flag;
		double vl_const;
		int initial_conditions_flag;
		int interior_points_normal_flag;
		double shapes;
		double rho_inf;
		double pr_inf;
		int threadsperblock;
		double gamma;
		int clcd_flag;

		Core() {}

	};

	struct PointConfig
	{
		int wall, interior, outer;

		PointConfig() {}
	};

	struct Format
	{
		char* type;

		Format() {}
	};


	Core core = Core();
	PointConfig point_config = PointConfig();
	Format format = Format();

	Config() {}

	void setConfig()
	{
		core.points = 9600;
		core.cfl = 0.5;
		core.max_iters = 10;
		core.mach = 1.2;
		core.aoa = 0.0;
		core.power = 0.0;
		core.limiter_flag = 1;
		core.vl_const = 200.0;
		core.initial_conditions_flag = 0;
		core.interior_points_normal_flag = 0;
		core.shapes = 1.0;
		core.rho_inf = 1.0;
		core.pr_inf = 0.7142857142857142;
		core.threadsperblock = 128;
		core.gamma = 1.4;
		core.clcd_flag = 0;

		point_config.wall = 0;
		point_config.interior = 1;
		point_config.outer = 2;

		format.type = "quadtree";

	}

	void printConfig()
	{
		cout<<"Core points: "<<core.points<<endl;
		cout<<"Core cfl: "<<core.cfl<<endl;
		cout<<"Core max_iters: "<<core.max_iters<<endl;
		cout<<"Core mach: "<<core.mach<<endl;
		cout<<"Core aoa: "<<core.aoa<<endl;
		cout<<"Core power: "<<core.power<<endl;
		cout<<"Core limiter_flag: "<<core.limiter_flag<<endl;
		cout<<"Core vl_const: "<<core.vl_const<<endl;
		cout<<"Core initial_conditions_flag: "<<core.initial_conditions_flag<<endl;
		cout<<"Core interior_points_normal_flag: "<<core.interior_points_normal_flag<<endl;
		cout<<"Core shapes: "<<core.shapes<<endl;
		cout<<"Core rho_inf: "<<core.rho_inf<<endl;
		cout<<"Core pr_inf: "<<core.pr_inf<<endl;
		cout<<"Core threadsperblock: "<<core.threadsperblock<<endl;
		cout<<"Core gamma: "<<core.gamma<<endl;
		cout<<"Core clcd_flag: "<<core.clcd_flag<<endl;

		cout<<"Config wall: "<<point_config.wall<<endl;
		cout<<"Config interior: "<<point_config.interior<<endl;
		cout<<"Config outer: "<<point_config.outer<<endl;

		cout<<"Format type: "<<format.type<<endl;

	}
};

double calculateTheta(Config configData);

void getInitialPrimitive(Config configData, double primal[4]);

void placeNormals(Point* globaldata, int idx, Config configData, long long interior, long long wall, long long outer);

xy_tuple calculateNormals(xy_tuple left, xy_tuple right, double mx, double my);

void calculateConnectivity(Point* globaldata, int idx);

void fpi_solver(int iter, Point* globaldata, Config configData, double res_old[1], int numPoints, double main_store[62], double tempdq[][2][4]);

void q_variables(Point* globaldata, int numPoints, double q_result[4]);

void q_var_derivatives(Point* globaldata, int numPoints, double power, double sig_del_x_del_q[4], double sig_del_y_del_q[4], double max_q[4], double min_q[4]);

void q_var_derivatives_innerloop(Point* globaldata, int numPoints, double power, double tempdq[][2][4], double sig_del_x_del_q[4], double sig_del_y_del_q[4], double qi_tilde[4], double qk_tilde[4]);

void debug_globaldata(Point* globaldata, int numPoints, int iter, int rk, double main_store[62]);

void printDebug(Point* globaldata, int numPoints, Config configData, int iter, double res_old[1], int rk, double sig_del_x_del_f[4], double sig_del_y_del_f[4], double main_store[62]);

void debug_Gs_and_qtildes(int iter, int rk, double Gxp[4], double Gxn[4], double Gyp[4], double Gyn[4], double phi_i[4], double phi_k[4], double G_i[4], double G_k[4], double result[4], double qtilde_i[4], double qtilde_k[4], double sig_del_x_del_f[4], double sig_del_y_del_f[4], double main_store[62]);

void debug_main_store_3(double main_store[62]);

template <class Type>
bool isNan(Type var);

#endif