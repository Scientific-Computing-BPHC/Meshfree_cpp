#include<stdio.h>
#include<iostream>
#include<math.h>
#include <unistd.h> 
#include <armadillo>
// #include <json/json.h>
// #include <json/value.h>
#include <fstream>
#include<map>
#include<string>
#include<iomanip>
#include<limits>

using namespace arma;

typedef std::map<char, float> char_float_map;
typedef std::map<char, double> char_double_map;
//typedef std::numeric_limits< double > dbl;


struct Point
{
	int localID;
	double x, y;
	int left, right;
	int flag_1, flag_2; //Int8 in the Julia code
	double short_distance;
	int nbhs;
	int conn[20];
	double nx, ny;
	double prim[4];
	double flux_res[4];
	double q[4];
	double dq1[4];
	double dq2[4];
	double entropy;
	int xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs;
	int xpos_conn[20];
	int xneg_conn[20];
	int ypos_conn[20];
	int yneg_conn[20];
	double delta;
	double max_q[20];
	double min_q[20];
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

	};

	struct PointConfig
	{
		int wall, interior, outer;
	};

	struct Format
	{
		char* type;
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
		core.limiter_flag = 0;
		core.vl_const = 200.0;
		core.initial_conditions_flag = 0;
		core.interior_points_normal_flag = 0;
		core.shapes = 0.0;
		core.rho_inf = 0.0;
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
		cout<<"Core threadsperblock: "<<core.threadsperblock<<endl;
		cout<<"Core gamma: "<<core.gamma<<endl;
		cout<<"Core clcd_flag: "<<core.clcd_flag<<endl;

		cout<<"Config wall: "<<point_config.wall<<endl;
		cout<<"Config interior: "<<point_config.interior<<endl;
		cout<<"Config outer: "<<point_config.outer<<endl;

		cout<<"Format type: "<<format.type<<endl;

	}
};

void meshfree_solver(char* filename, int num_iters);
//void getConfig(Json::Value &config);
//void set_manual_config(std::map<char, char_float_map> &config);



int main(int argc, char **argv)
{
	printf("Meshfree AD\n");

	/* initialize random seed*/
	srand (time(NULL));
	arma_rng::set_seed_random();

	meshfree_solver(argv[1], (int)argv[2]); //Warning: Casting from char to int loses precision
	//Gotta see if maintaining a global array id efficient or if passing around by reference is efficient
	//For all we know, maintaining a global data structure instead of passing it around might be more efficient
}

#if 0

void getConfig(Json::Value &config)
{
	std::ifstream config_file("../config.json", std::ifstream::binary);
	config_file>>config;
	cout<<config;
	config["format"] = "quadtree";
}

void set_manual_config(std::map<char, char_double_map> &config)
{
	char_double_map core_config =
	{
		{'1', 9600},
        {'2', 0.5},
        {'3', 10},
        {'4', 1.2},
        {'5', 0.0},
        {'6', 0.0},
        {'7', 1},
        {'8', 200.0},
        {'9', 0},
        {'A', 0},
        {'B', 1.0},
        {'C', 1.0},
        {'D', 0.7142857142857142},
        {'E', 128},
        {'F', 1.4},
        {'G', 0}
	};

	char_double_map point_config = 
	{
		{'1', 0},
		{'2', 1},
		{'3', 2}
	}; // Gotta typecast these values back to int, or use templating appropriately later.

	char_double_map format_config =
	{
		{'1', 1} // Type 1 refers to Quadtree, Type 0 refers to non-quadtree
	};

	config = 
	{
	    { '1', core_config },
	    { '2', point_config },
	    { '3', format_config }
	};

For reference:
{
    "core":{
        "points":9600,  // 1
        "cfl": 0.5, // 2
        "max_iters": 10, // 3
        "mach": 1.2, // 4
        "aoa": 0.0, // 5
        "power": 0.0, // 6
        "limiter_flag": 1, // 7
        "vl_const": 200.0, // 8
        "initial_conditions_flag": 0, // 9
        "interior_points_normal_flag": 0, // 10
        "shapes": 1.0, // 11
        "rho_inf": 1.0, // 12
        "pr_inf": 0.7142857142857142, // 13
        "threadsperblock": 128, // 14
        "gamma": 1.4, // 15
        "clcd_flag": 0 // 16
    },
    "point":{
        "wall": 0,
        "interior": 1,
        "outer": 2
    },
    "format":{
        "type":"quadtree"
    }
}


}

#endif


void meshfree_solver(char* filename, int max_iters)
{
	//Json::Value configData = getConfig(Json::Value configData); // Assuming config file is located at "../config.json"
	//std::cout << std::setprecision(17); // Show 17 digits of precision
	//cout.precision(17);

#if 0

	std::map<char, char_double_map> configData;
	set_manual_config(configData);

	if(configData['1']['F'] == 1.45)
		cout<<"All good"<<endl;
	else
		cout<<"Values Changed"<<endl;		
	std::cout<<"Checking config read: "<<endl;
	cout<<"gamma "<<configData['1']['F']<<endl;

	std::streamsize ss = std::cout.precision();
	cout.precision(dbl::max_digits10);
	cout<<std::fixed<<"pr_inf "<<configData['1']['D']<<endl; // Shows 0.71428...29. Why am I getting the extra 9

	cout.precision(ss);
	cout<<"gamma "<<configData['1']['F']<<endl;
	cout<<"threadsperblock "<<configData['1']['E']<<endl;

	if(configData['1']['F'] == 1.39999997615814209)
		cout<<"Values Changed"<<endl;
	else
		cout<<"All good"<<endl;

	// For some reason this double works well but float doesn't
	// Please note that setprecision only changes what we display through cout, the internal representation used for comparisons anol is still ABSOLUTELY THE SAME.

	double gamma = 1.4;
	cout<<gamma<<endl;

	double pr_inf = 0.7142857142857142;
	cout<<pr_inf<<endl;

	if(gamma == 1.4)
		cout<<"All good"<<endl;
	else
		cout<<"Values Changed"<<endl;

	if(pr_inf == 0.7142857142857142)
		cout<<"All good"<<endl;
	else
		cout<<"Values Changed"<<endl;

	cout.precision(dbl::max_digits10);

	cout<<pr_inf<<endl;

	if(gamma == 1.4)
		cout<<"All good"<<endl;
	else
		cout<<"Values Changed"<<endl;

	if(pr_inf == 0.7142857142857142)
		cout<<"All good"<<endl;
	else
		cout<<"Values Changed"<<endl;

#endif

	Config configData = Config();
	configData.setConfig();
	configData.printConfig();

}