// Please repeatedly remind yourself that whatever starts with index x in Julia starts with index x-1 in C++
// Again, and again remind youself

// Learn to use memset to initialize all the points in a single shot

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

#include "utils.hpp"
#include "core.hpp"

#define ARMA_DONT_PRINT_ERRORS
using namespace arma;

typedef std::map<char, float> char_float_map;
typedef std::map<char, double> char_double_map;
typedef std::vector<std::string> vec_str;
typedef std::vector<double> vec_doub;
typedef std::vector<long double> vec_ldoub;

bool debug_mode = true;

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

void meshfree_solver(char* file_name, int num_iters);

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

void meshfree_solver(char* file_name, int max_iters)
{

	Config configData = Config();
	configData.setConfig();
	if (debug_mode)
		configData.printConfig();

	std::string format = configData.format.type;
	if (debug_mode)
		cout<<"Format: "<<format<<endl;

	cout<<"Filename: "<<file_name<<endl;
	cout<<"hi";

	int numPoints = 0;
	std::fstream datafile(file_name, ios::in);

	if(datafile.is_open())
	{
		std::string temp;
		getline(datafile, temp);
		std::stringstream num(temp);
		num >> numPoints;	
	}

	std::string new_file[numPoints], temp;
	if(datafile.is_open())
	{
		for(int i=0; i<numPoints; i++)
		{
			getline(datafile, new_file[i]);
		}
	}
	datafile.close();
	
	//cout<<new_file[1][0]<<endl;
	cout<<"No. of points: "<<numPoints<<endl;

	std::regex ws_re("\\s+"); 
	std::vector<vec_str> result;
	for(int i=0; i<numPoints; i++)  // This might include the first line as well so, you need to store that separately or just throw the 1st line away (the one that contains the file length)
									// There are 48738 lines apart from that
	{ 
    	std::vector<std::string> temp{std::sregex_token_iterator(new_file[i].begin(), new_file[i].end(), ws_re, -1), {}};
    	result.push_back(temp);
    }
#if 0

    for(int j=0; j<numPoints; j++)
    {
		for (int i=0; i<result[j].size(); i++)
			cout<<"Result: "<<j<<" " <<i<<" "<<result[j][i]<<endl;

	}
	
#endif

	std::vector<vec_doub> result_doub;
	for(int j=0; j<numPoints; j++)
    {
		vec_doub temp;
		for (int i=0; i<result[j].size(); i++)
			temp.push_back(std::stod(result[j][i]));
		result_doub.push_back(temp);

	}


#if 0
	
	checkFileRead(result_doub, numPoints);

	for(int j=0; j<20; j++)
    {
		for (int i=0; i<result_doub[j].size(); i++)
			cout<<std::fixed<<std::setprecision(20)<<"Result Doub: "<<j<<" " <<i<<" "<<result_doub[j][i]<<endl;

	}
#endif

	Point* globaldata = new Point[numPoints];
	double res_old = 0.0;
	double main_store[62] = {0};

	double defprimal[4];
	getInitialPrimitive(configData, defprimal);
	printPrimal(defprimal);

	cout<<"-----Start Read-----"<<endl; // readFile function in Julia

	//cout<<result_doub[0][0]<<endl;

	for(int idx=0; idx<numPoints; idx++)
	{

		int connectivity[20] = {0};
		for(int iter=8; iter<result_doub[idx].size(); iter++)
		{
			connectivity[iter-8] = result_doub[idx][iter];
		}

		// Assign values to each of the elements of global data
		// Make sure the types specifications are adhered to

		globaldata[idx].localID = idx;
		globaldata[idx].x = result_doub[idx][0];
		globaldata[idx].y = result_doub[idx][1];
		globaldata[idx].left = (int)result_doub[idx][2];
		globaldata[idx].right = (int)result_doub[idx][3];
		globaldata[idx].flag_1 = (int)result_doub[idx][4];
		globaldata[idx].flag_2 = (int)result_doub[idx][5];
		globaldata[idx].short_distance = result_doub[idx][6];
		globaldata[idx].nbhs = (int)result_doub[idx][7];

		for(int i=0; i<20; i++)
			globaldata[idx].conn[i] = connectivity[i];

		globaldata[idx].nx = 0.0;
		globaldata[idx].ny = 0.0;

		for(int i=0; i<4; i++)
		{
			globaldata[idx].prim[i] = defprimal[i];
			globaldata[idx].flux_res[i] = 0.0;
			globaldata[idx].q[i] = 0.0;
			globaldata[idx].dq1[i] = 0.0;
			globaldata[idx].dq2[i] = 0.0;
			globaldata[idx].max_q[i] = 0.0;
			globaldata[idx].min_q[i] = 0.0;
			globaldata[idx].prim_old[i] = 0.0;
		}

		globaldata[idx].xpos_nbhs = 0;
		globaldata[idx].xneg_nbhs = 0;
		globaldata[idx].ypos_nbhs = 0;
		globaldata[idx].yneg_nbhs = 0;

		globaldata[idx].entropy = 0.0;
		globaldata[idx].delta = 0.0;


		for(int i=0; i<20; i++)
		{
			globaldata[idx].xpos_conn[i] = 0;
			globaldata[idx].xneg_conn[i] = 0;
			globaldata[idx].ypos_conn[i] = 0;
			globaldata[idx].yneg_conn[i] = 0;
		}

	}

	cout<<"-----End Read-----"<<endl;

	// Interior, Out and Wall were defined as Int64 in Julia, so am defining them as long long

	long long interior = configData.point_config.interior;
	long long wall = configData.point_config.wall;
	long long outer = configData.point_config.outer;




}	