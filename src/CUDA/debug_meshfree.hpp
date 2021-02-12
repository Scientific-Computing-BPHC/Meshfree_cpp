#ifndef DEBUG_MESHFREE
#define DEBUG_MESHFREE
#include "point.hpp"

void printNormals(int idx_low, int idx_high, Point* globalData);

void printConnectivity(int idx_low, int idx_high, int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn);

#endif