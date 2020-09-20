#ifndef POINT_HPP
#define POINT_HPP

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

#include "core_cuda.hpp"

#define ARMA_DONT_PRINT_ERRORS
using namespace arma;

typedef std::map<char, float> char_float_map;
typedef std::map<char, double> char_double_map;
typedef std::vector<std::string> vec_str;
typedef std::vector<double> vec_doub;
typedef std::vector<long double> vec_ldoub;
typedef std::tuple <double, double> xy_tuple;
typedef std::tuple <double, double, double, double, double, double, double> conn_tuple;

xy_tuple getxy(Point self);
void setNormals(Point* globaldata, int idx, xy_tuple n);

#endif
