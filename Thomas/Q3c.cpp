#include <iostream>
#include <cmath>
#include "EulerMethod.hpp"

double* f(double* x);

int main(int argc, char* argv[])
{
    double* initial;
    initial = new double[2];
    initial[0] = 2;
    initial[1] = 1;

    double h = 1.0/12.0;

    int n = 24;

    double* truesol;
    truesol = new double[2];
    truesol[0] = 1;
    truesol[1] = 1;

    std::cout << "n" << "  " << "y_1,n" << "  " << "y_2,n" << "  " 
    << "||y_n -y_*||" << std::endl; 

    for(int i=0; i<n; i++)
    {
        std::cout << i+1 << "  " << ComputeEulerIteration(i+1, initial, f, h)[0] 
        << "  " << ComputeEulerIteration(i+1, initial, f, h)[1] << "  " 
        << Compute2Norm(ComputeEulerIteration(i+1, initial, f, h), truesol, 2) 
        << std::endl;
    }

    return 0;
}

double* f(double* x)
{
    double* y;
    y = new double[2];

    y[0] = -3*pow(x[0], 2)+ 2*x[1] + 1;
    y[1] = -3*pow(x[1], 2)+ 2*x[0] + 1;

    return y;
}

