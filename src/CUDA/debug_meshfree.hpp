#ifndef DEBUG_MESHFREE
#define DEBUG_MESHFREE
#include "point.hpp"

void printNormals(int idx_low, int idx_high, Point* globalData);

void printConnectivity(int idx_low, int idx_high, int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn);

__global__ void printFuncDelta(int idx_low, int idx_high, Point* globalData);

__global__ void print_qVar_cuda(int idx_low, int idx_high, double* q);

__global__ void print_qVarDer_cuda(int idx_low, int idx_high, double* dq_1, double* dq_2);

__global__ void print_qVarDer_after_innerloop_cuda(int idx_low, int idx_high, double* dq_1, double* dq_2);

#endif