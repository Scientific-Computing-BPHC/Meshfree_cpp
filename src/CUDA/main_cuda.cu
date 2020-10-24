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
#include <chrono>

#include "utils.hpp"
#include "core_cuda.hpp"
#include "point.hpp"

#define ARMA_DONT_PRINT_ERRORS
using namespace arma;

typedef std::map<char, float> char_float_map;
typedef std::map<char, double> char_double_map;
typedef std::vector<std::string> vec_str;
typedef std::vector<double> vec_doub;
typedef std::vector<long double> vec_ldoub;

bool debug_mode = true;

void meshfree_solver(char* file_name, int num_iters);
void run_code(Point* globaldata, Config configData, double res_old[1], int numPoints, TempqDers* tempdq, int max_iters, int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn, \
	int* connec, double* prim, double* flux_res, double* q, double* dq1, double* dq2, double* max_q, double* min_q, double* prim_old);
void test_code(Point* globaldata, Config configData, double res_old[1], int numPoints, int max_iters, int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn,\
	int* connec, double* prim, double* flux_res, double* q, double* dq1, double* dq2, double* max_q, double* min_q, double* prim_old);

int main(int argc, char **argv)
{
	printf("\nMeshfree AD\n");

	/* initialize random seed*/
	srand (time(NULL));
	arma_rng::set_seed_random();

	meshfree_solver(argv[1], std::stoi(argv[2])); //Warning: Casting from char to int loses precision
	//Gotta see if maintaining a global array id efficient or if passing around by reference is efficient
	//For all we know, maintaining a global data structure instead of passing it around might be more efficient

	cout<<"\n Max Iters: "<<std::stoi(argv[2])<<endl;
}

void meshfree_solver(char* file_name, int max_iters)
{

	Config configData = Config();
	configData.setConfig();
	if (debug_mode)
		configData.printConfig();

	std::string format = configData.format.type;
	if (debug_mode)
		cout<<"\nFormat: "<<format<<endl;

	cout<<"\nFilename: "<<file_name<<endl;
	//cout<<"hi"<<endl;

	int numPoints = 0;
	std::fstream datafile(file_name, ios::in);

	if(datafile.is_open())
	{
		std::string temp;
		getline(datafile, temp);
		std::stringstream num(temp);
		num >> numPoints;	
	}

	std::string temp;
	std::string* new_file = new std::string[numPoints];
	if(datafile.is_open())
	{
		for(int i=0; i<numPoints; i++)
		{
			getline(datafile, new_file[i]);
		}
	}
	datafile.close();
	
	//cout<<new_file[1][0]<<endl;
	cout<<"\nNo. of points: "<<numPoints<<endl;

	std::regex ws_re("\\s+"); 
	std::vector<vec_str> result;
	for(int i=0; i<numPoints; i++)  // This might include the first line as well so, you need to store that separately or just throw the 1st line away (the one that contains the file length)
									// There are 48738 lines apart from that
	{ 
    	std::vector<std::string> temp{std::sregex_token_iterator(new_file[i].begin(), new_file[i].end(), ws_re, -1), {}};
    	result.push_back(temp);
	}
	
	// Free up the space taken by new_file
	delete[] new_file;

	std::vector<vec_doub> result_doub;
	for(int j=0; j<numPoints; j++)
    {
		vec_doub temp;
		for (int i=0; i<result[j].size(); i++)
			temp.push_back(std::stod(result[j][i]));
		result_doub.push_back(temp);

	}

	std::vector<vec_str>().swap(result); // Free up the space taken up by result

	Point* globaldata = new Point[numPoints];
	int* xpos_conn = new int[numPoints*20];
	int* xneg_conn = new int[numPoints*20];
	int* ypos_conn = new int[numPoints*20];
	int* yneg_conn = new int[numPoints*20];
	double res_old[1] = {0.0};

	int* connec = new int[numPoints*20];

	double* prim = new double[numPoints*4];
	double* flux_res = new double[numPoints*4];
	double* q = new double[numPoints*4];
	double* dq1 = new double[numPoints*4];
	double* dq2 = new double[numPoints*4];
	double* max_q = new double[numPoints*4];
	double* min_q = new double[numPoints*4];
	double* prim_old = new double[numPoints*4];

	double defprimal[4];
	getInitialPrimitive(configData, defprimal);
	//printPrimal(defprimal);

	cout<<"\n-----Start Read-----\n"; // readFile function in Julia

	//cout<<result_doub[0][0]<<endl;


/* Read File Legacy */

#if 0
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
		{
			connec[idx*20 +i] = connectivity[i];
			// (connectivity[i]!=0) cout<<"\n non-zero connectivity \n";
		}


		globaldata[idx].nx = 0.0;
		globaldata[idx].ny = 0.0;

		for(int i=0; i<4; i++)
		{
			prim[idx*4 + i] = defprimal[i];
			flux_res[idx*4 + i] = 0.0;
			q[idx*4 + i] = 0.0;
			dq1[idx*4 + i] = 0.0;
			dq2[idx*4 + i] = 0.0;
			max_q[idx*4 + i] = 0.0;
			min_q[idx*4 + i] = 0.0;
			prim_old[idx*4 + i] = 0.0;
		}

		globaldata[idx].xpos_nbhs = 0;
		globaldata[idx].xneg_nbhs = 0;
		globaldata[idx].ypos_nbhs = 0;
		globaldata[idx].yneg_nbhs = 0;

		globaldata[idx].entropy = 0.0;
		globaldata[idx].delta = 0.0;


		for(int i=0; i<20; i++)
		{
			xpos_conn[idx*20 + i] = 0;
			xneg_conn[idx*20 + i] = 0;
			ypos_conn[idx*20 + i] = 0;
			yneg_conn[idx*20 + i] = 0; 
		}

	}

	cout<<"\n-----End Read-----\n";

#endif

/* Read File Quadtree */

	for(int idx=0; idx<numPoints; idx++)
	{

		int connectivity[20] = {0};
		for(int iter=11; iter<result_doub[idx].size(); iter++)
		{
			connectivity[iter-11] = result_doub[idx][iter];
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
		globaldata[idx].short_distance = result_doub[idx][9];
		globaldata[idx].nbhs = (int)result_doub[idx][10];

		for(int i=0; i<20; i++)
		{
			connec[idx*20 +i] = connectivity[i];
			// (connectivity[i]!=0) cout<<"\n non-zero connectivity \n";
		}


		globaldata[idx].nx = result_doub[idx][6];
		globaldata[idx].ny = result_doub[idx][7];

		for(int i=0; i<4; i++)
		{
			prim[idx*4 + i] = defprimal[i];
			flux_res[idx*4 + i] = 0.0;
			q[idx*4 + i] = 0.0;
			dq1[idx*4 + i] = 0.0;
			dq2[idx*4 + i] = 0.0;
			max_q[idx*4 + i] = 0.0;
			min_q[idx*4 + i] = 0.0;
			prim_old[idx*4 + i] = 0.0;
		}

		globaldata[idx].xpos_nbhs = 0;
		globaldata[idx].xneg_nbhs = 0;
		globaldata[idx].ypos_nbhs = 0;
		globaldata[idx].yneg_nbhs = 0;

		globaldata[idx].entropy = 0.0;
		globaldata[idx].delta = 0.0;


		for(int i=0; i<20; i++)
		{
			xpos_conn[idx*20 + i] = 0;
			xneg_conn[idx*20 + i] = 0;
			ypos_conn[idx*20 + i] = 0;
			yneg_conn[idx*20 + i] = 0; 
		}

	}

	cout<<"\n-----End Read-----\n";
	if(configData.core.restart != 1) std::vector<vec_doub>().swap(result_doub);
	
	if(configData.core.restart == 1)
    {
		// Check if any of the memory deallocations with delete[] new_file, or std::vector<vec_doub>().swap(result_doub) etc. are causing issues for the restart version
        char* file_name = "/home/hari/Work/Meshfree_cpp/restart.dat";
        std::fstream datafile(file_name, ios::in);
        std::string temp;
        std::regex ws_re("\\s+"); 

        int numPoints, max_iters;
        double residue;

        if(datafile.is_open())
        {
            
            getline(datafile, temp);
            std::vector<std::string> temp_vec{std::sregex_token_iterator(temp.begin(), temp.end(), ws_re, -1), {}};
            numPoints = (int) std::stod(temp_vec[1]);
        }

        std::string new_file[numPoints];
        if(datafile.is_open())
        {
            for(int i=0; i<numPoints; i++)
            {
                getline(datafile, new_file[i]);
            }
        }
        datafile.close();

        std::vector<vec_str> result;
        for(int i=0; i<numPoints; i++)  // This might include the first line as well so, you need to store that separately or just throw the 1st line away (the one that contains the file length)
                                        // There are 48738 lines apart from that
        { 
            std::vector<std::string> temp{std::sregex_token_iterator(new_file[i].begin(), new_file[i].end(), ws_re, -1), {}};
            result.push_back(temp);
        }

        std::vector<vec_doub> result_doub;

        for(int j=0; j<numPoints; j++)
        {
            vec_doub temp;
            for (int i=1; i<result[j].size(); i++)
                temp.push_back(std::stod(result[j][i]));
            result_doub.push_back(temp);
        }

        for(int i=0; i<numPoints; i++)
        {
        	for(int j=0; j<4; j++)
        		prim[i*4 + j] = result_doub[i][5+j];
        }
    }

	// Interior, Out and Wall were defined as Int64 in Julia, so am defining them as long long

	long long interior = configData.point_config.interior;
	long long wall = configData.point_config.wall;
	long long outer = configData.point_config.outer;


	cout<<"\n-----Computing Normals-----\n";

	for(int idx=0; idx<numPoints; idx++)
		placeNormals(globaldata, idx, configData, interior, wall, outer);

	cout<<"\n-----Start Connectivity Generation-----\n";

	for(int idx=0; idx<numPoints; idx++)
		calculateConnectivity(globaldata, idx, xpos_conn, xneg_conn, ypos_conn, yneg_conn, connec);

	cout<<"\n-----Connectivity Generation Done-----\n";  

	cout<<"\n"<<max_iters+1<<endl;

	test_code(globaldata, configData, res_old, numPoints, max_iters, xpos_conn, xneg_conn, ypos_conn, yneg_conn, connec, prim, flux_res, q, dq1, dq2, max_q, min_q, prim_old); 


	cout<<"\n--------Done--------\n"<<endl;

}	


void run_code(Point* globaldata, Config configData, double res_old[1], int numPoints, TempqDers* tempdq, int max_iters, int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn, \
	int* connec, double* prim, double* flux_res, double* q, double* dq1, double* dq2, double* max_q, double* min_q, double* prim_old)
{
	auto begin = std::chrono::high_resolution_clock::now();
	// I'M ASSUMING IN FORTRAN ALSO THEY MEASURE THE COPY TIME AND NOT JUST COMPUTATION TIME. Just doesn't make sense to ignore copy time, although it won't be too huge
	cudaStream_t stream;  
    Point* globaldata_d;
    unsigned int mem_size_A = sizeof(struct Point) * numPoints;
    TempqDers* tempdq_d;
    unsigned int mem_size_B = sizeof(struct TempqDers) * numPoints;
    double* res_sqr_d, *res_old_d;
    unsigned int mem_size_C = sizeof(double);
	unsigned int mem_size_D = sizeof(double) * numPoints;
	unsigned int mem_size_E = sizeof(int) * 20 * numPoints;
	unsigned int mem_size_F = sizeof(double) * 4 * numPoints;
	int* xpos_conn_d, *xneg_conn_d, *ypos_conn_d, *yneg_conn_d;
	int* connec_d;
	double* prim_d, *flux_res_d, *q_d, *dq1_d, *dq2_d, *max_q_d, *min_q_d, *prim_old_d;

    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&globaldata_d), mem_size_A));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&tempdq_d), mem_size_B)); 
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&res_old_d), mem_size_C));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&res_sqr_d), mem_size_D));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&xpos_conn_d), mem_size_E));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&xneg_conn_d), mem_size_E));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&ypos_conn_d), mem_size_E));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&yneg_conn_d), mem_size_E));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&connec_d), mem_size_E));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&prim_d), mem_size_F));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&flux_res_d), mem_size_F));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&q_d), mem_size_F));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&dq1_d), mem_size_F));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&dq2_d), mem_size_F));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&max_q_d), mem_size_F));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&min_q_d), mem_size_F));
	checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&prim_old_d), mem_size_F));
	checkCudaErrors(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));

    double* res_sqr = new double[numPoints];

    // Copy from host to device
    checkCudaErrors(cudaMemcpyAsync(globaldata_d, globaldata, mem_size_A, cudaMemcpyHostToDevice, stream));
	checkCudaErrors(cudaMemcpyAsync(tempdq_d, tempdq, mem_size_B, cudaMemcpyHostToDevice, stream));  
	checkCudaErrors(cudaMemcpyAsync(xpos_conn_d, xpos_conn, mem_size_E, cudaMemcpyHostToDevice, stream));  
	checkCudaErrors(cudaMemcpyAsync(xneg_conn_d, xneg_conn, mem_size_E, cudaMemcpyHostToDevice, stream));  
	checkCudaErrors(cudaMemcpyAsync(ypos_conn_d, ypos_conn, mem_size_E, cudaMemcpyHostToDevice, stream));  
	checkCudaErrors(cudaMemcpyAsync(yneg_conn_d, yneg_conn, mem_size_E, cudaMemcpyHostToDevice, stream));  
	checkCudaErrors(cudaMemcpyAsync(connec_d, connec, mem_size_E, cudaMemcpyHostToDevice, stream)); 
	checkCudaErrors(cudaMemcpyAsync(prim_d, prim, mem_size_F, cudaMemcpyHostToDevice, stream)); 
	checkCudaErrors(cudaMemcpyAsync(flux_res_d, flux_res, mem_size_F, cudaMemcpyHostToDevice, stream));
	checkCudaErrors(cudaMemcpyAsync(q_d, q, mem_size_F, cudaMemcpyHostToDevice, stream));
	checkCudaErrors(cudaMemcpyAsync(dq1_d, dq1, mem_size_F, cudaMemcpyHostToDevice, stream));
	checkCudaErrors(cudaMemcpyAsync(dq2_d, dq2, mem_size_F, cudaMemcpyHostToDevice, stream));
	checkCudaErrors(cudaMemcpyAsync(max_q_d, max_q, mem_size_F, cudaMemcpyHostToDevice, stream));
	checkCudaErrors(cudaMemcpyAsync(min_q_d, min_q, mem_size_F, cudaMemcpyHostToDevice, stream));
	checkCudaErrors(cudaMemcpyAsync(prim_old_d, prim_old, mem_size_F, cudaMemcpyHostToDevice, stream));
	
	for (int i=0; i<max_iters; i++)
	{
		fpi_solver(i, globaldata_d, configData, res_old_d, res_sqr_d, numPoints, tempdq_d, stream, res_old, res_sqr, mem_size_C, mem_size_D, \
			xpos_conn_d, xneg_conn_d, ypos_conn_d, yneg_conn_d, connec_d, prim_d, flux_res_d, q_d, dq1_d, dq2_d, max_q_d, min_q_d, prim_old_d); 

	}
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
	printf("\nTime measured: %.5f seconds.\n", elapsed.count() * 1e-9);
	// Copy from device to host, and free the device memory
	checkCudaErrors(cudaMemcpyAsync(globaldata, globaldata_d, mem_size_A, cudaMemcpyDeviceToHost, stream));
	checkCudaErrors(cudaMemcpyAsync(tempdq, tempdq_d, mem_size_B, cudaMemcpyDeviceToHost, stream));
	// Do we wanna copy more stuff here? ...
    checkCudaErrors(cudaFree(globaldata_d)); 
    checkCudaErrors(cudaFree(tempdq_d));
    checkCudaErrors(cudaFree(res_old_d));
	checkCudaErrors(cudaFree(res_sqr_d));
	checkCudaErrors(cudaFree(xpos_conn_d));
	checkCudaErrors(cudaFree(xneg_conn_d));
	checkCudaErrors(cudaFree(ypos_conn_d));
	checkCudaErrors(cudaFree(yneg_conn_d));
	checkCudaErrors(cudaFree(connec_d));
	checkCudaErrors(cudaFree(prim_d));
	checkCudaErrors(cudaFree(flux_res_d));
	checkCudaErrors(cudaFree(q_d));
	checkCudaErrors(cudaFree(dq1_d));
	checkCudaErrors(cudaFree(dq2_d));
	checkCudaErrors(cudaFree(max_q_d));
	checkCudaErrors(cudaFree(min_q_d));
	checkCudaErrors(cudaFree(prim_old_d));

}


void test_code(Point* globaldata, Config configData, double res_old[1], int numPoints, int max_iters, int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn, \
int* connec, double* prim, double* flux_res, double* q, double* dq1, double* dq2, double* max_q, double* min_q, double* prim_old)
{
	cout<<"\nStarting warmup function \n";
	res_old[0] = 0.0;

	cout<<"\nStarting main function \n";
	TempqDers* tempdq = new TempqDers[numPoints];
	
	for (int i=0; i<numPoints; i++)
		tempdq[i].setTempdq();

	run_code(globaldata, configData, res_old, numPoints, tempdq, max_iters, xpos_conn, xneg_conn, ypos_conn, yneg_conn, connec, prim, flux_res, q, dq1, dq2, max_q, min_q, prim_old); 

}