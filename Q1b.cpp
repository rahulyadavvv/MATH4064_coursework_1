#include <iostream>
#include "TridiagonalSolver.hpp"

int main(int argc, char* argv[])
{
    int n = 4;

    double* lower; 
    double* diagonal;
    double* upper;
    double* rhs;

    lower = new double[n];
    diagonal = new double[n];
    upper = new double[n];
    rhs = new double[n];

    lower[0] = 0;
    lower[1] = 2;
    lower[2] = 1;
    lower[3] = 1;

    diagonal[0] = 6;
    diagonal[1] = 4;
    diagonal[2] = 4;
    diagonal[3] = 6;

    upper[0] = 1;
    upper[1] = 1;
    upper[2] = 2;
    upper[3] = 0;

    rhs[0] = 8;
    rhs[1] = 13;
    rhs[2] = 22;
    rhs[3] = 27;

    double* x = SolveTridiagonalSystem(n, lower, diagonal, upper, rhs);

    for(int i=0; i<n; i++)
    {
        std::cout << x[i] << std::endl;
    }
    std::cout << std::endl;

    return 0;
}