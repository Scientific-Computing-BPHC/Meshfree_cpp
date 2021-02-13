#include "readHDF5.hpp"

void printHDF5file(double* result_doub_arr, int max_row, int max_col)
{
    for(int i=0; i<max_row; i++)
    {
        cout<<"\n";
        for(int j=0; j<max_col; j++)
            cout<<result_doub_arr[i*max_col + j]<<"  ";
    }
}

void saveHDF5file(double* result_doub_arr, int max_row, int max_col)
{
    printf("\n \nSaving file....");
    std::ofstream fout("HDF5file.txt");
    for(int i=0; i<max_row; i++)
    {
        fout<<"\n";
        for(int j=0; j<max_col; j++)
            fout<<result_doub_arr[i*max_col + j]<<"  ";
    }
    fout.close();
}

double* readHDF5file(string file_path, string dataset_name, string data_filename, string folder, double* result_doub_arr, int& numPoints, int& maxCol)
{
    const H5std_string FILE_NAME(file_path);
    const H5std_string DATASET_NAME(dataset_name);

    const int RANK_OUT = 2;
    const int nParts = 1;

    /*
    * Try block to detect exceptions raised by any of the calls inside it
    */
    try
    {
        /*
        * Turn off the auto-printing when failure occurs so that we can
        * handle the errors appropriately
        */
        Exception::dontPrint();
        H5File file(FILE_NAME, H5F_ACC_RDONLY );
        DataSet dataset = file.openDataSet( DATASET_NAME );

        H5T_class_t type_class = dataset.getTypeClass();

        auto intype = dataset.getFloatType();
        size_t size = intype.getSize();
        DataSpace dataspace = dataset.getSpace();
        int rank = dataspace.getSimpleExtentNdims();
        /*
        * Get the dimension size of each dimension in the dataspace and
        * display them.
        */
        hsize_t dims_out[2];
        int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
        cout << "rank " << rank << ", dimensions " <<
            (unsigned long)(dims_out[0]) << " x " <<
            (unsigned long)(dims_out[1]) << endl;
        
        numPoints = dims_out[0];
        maxCol = dims_out[1];

        hsize_t      offset[2];   // hyperslab offset in the file
        hsize_t      count[2];    // size of the hyperslab in the file
        offset[0] = 0;
        offset[1] = 0;
        count[0]  = dims_out[0];
        count[1]  = dims_out[1];
        dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );
    
        /*
        * Define the memory dataspace.
        */
        hsize_t dimsm[2]; /* memory space dimensions */
        dimsm[0] = dims_out[0];
        dimsm[1] = dims_out[1];
        DataSpace memspace( RANK_OUT, dimsm );

        int i, j;
        result_doub_arr = new double[dims_out[0] * dims_out[1]]; /* output buffer */
        for (j = 0; j < dims_out[0]; j++)
        {
            for (i = 0; i < dims_out[1]; i++)
            {
                result_doub_arr[j*dims_out[1] + i]= 0;
            }
        }
        dataset.read(result_doub_arr, PredType::NATIVE_DOUBLE, memspace, dataspace);
    }
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printErrorStack();
        return -1;
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printErrorStack();
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printErrorStack();
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printErrorStack();
        return -1;
    }
    return result_doub_arr;
}