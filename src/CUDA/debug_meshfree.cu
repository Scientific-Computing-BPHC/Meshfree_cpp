#include "debug_meshfree.hpp"

void printNormals(int idx_low, int idx_high, Point* globalData)
{   
    cout<<"\n nx followed by ny \n";
    for(int i=idx_low; i<=idx_high; i++)
    {
        cout<<globalData[i].nx<<"  "<<globalData[i].ny<<endl;
    }
}

void printConnectivity(int idx_low, int idx_high, int* xpos_conn, int* xneg_conn, int* ypos_conn, int* yneg_conn)
{   
    cout<<"\n xpos | xneg | ypos | yneg \n";
    for(int i=idx_low; i<=idx_high; i++)
    {
        cout<<"\n ---- Point idx: "<<i<<" --- \n";
        for(int j=0; j<20; j++)
        {
            int access_idx = i*20 + j;
            cout<<xpos_conn[access_idx]<<"  "<<xneg_conn[access_idx]<<ypos_conn[access_idx]<<"  "<<yneg_conn[access_idx]<<endl;
        }
    }
}

__global__ void printFuncDelta(int idx_low, int idx_high, Point* globalData)
{
    printf("\nFuncDelta \n");
    for(int i=idx_low; i<=idx_high; i++)
    {
        printf("\n Point: %d", i);
        printf("\n %lf", globalData[i].nx);
        printf("\n");

    }
}

__global__ void print_qVar_cuda(int idx_low, int idx_high, double* q)
{
    printf("\nqVar \n");
    for(int i=idx_low; i<=idx_high; i++)
    {
        printf("\n Point: %d", i);
        for(int j=0; j<4; j++)
        {
            printf("%lf   ", q[i*4 +j]);
        }
        printf("\n");
    }
}

__global__ void print_qVarDer_cuda(int idx_low, int idx_high, double* dq_1, double* dq_2)
{
    printf("\nqVarDer \n");
    for(int i=idx_low; i<=idx_high; i++)
    {
        printf("\n Point: %d", i);
        for(int j=0; j<4; j++)
        {
            printf("%lf   ", dq_1[i*4 +j]);
        }
        printf("\n");
        for(int j=0; j<4; j++)
        {
            printf("%lf   ", dq_2[i*4 +j]);
        }
        printf("\n");

    }
}

__global__ void print_qVarDer_after_innerloop_cuda(int idx_low, int idx_high, double* dq_1, double* dq_2)
{
    printf("\nqVarDer \n");
    for(int i=idx_low; i<=idx_high; i++)
    {
        printf("\n Point: %d", i);
        for(int j=0; j<4; j++)
        {
            printf("%lf   ", dq_1[i*4 +j]);
        }
        printf("\n");
        for(int j=0; j<4; j++)
        {
            printf("%lf   ", dq_2[i*4 +j]);
        }
        printf("\n");

    }
}
