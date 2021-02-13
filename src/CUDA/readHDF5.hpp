#ifndef READ_HDF5_HPP
#define READ_HDF5_HPP

/* External Headers */
#include<iostream>
#include <string>
#include <fstream>
using std::cout;
using std::endl;
using std::string;

#include <H5Cpp.h>
using namespace H5;

double* readHDF5file(string file_path, string dataset_name, string data_filename, string folder, double* result_doub, int& numPoints, int& maxCol);
void printHDF5file(double* result_doub_arr, int max_row, int max_col);
void saveHDF5file(double* result_doub_arr, int max_row, int max_col);

#endif