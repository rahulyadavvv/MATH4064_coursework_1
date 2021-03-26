#include <iostream>
#include "TridiagonalMatrix.hpp"

int main(int argc, char* argv[])
{
    int n = 4;

    double* vecL;
    double* vecD;
    double* vecU;
    double* vecRHS;

    vecL = new double[n];
    vecD = new double[n];
    vecU = new double[n];
    vecRHS = new double[n];

    vecL[0] = 0;
    vecL[1] = 2;
    vecL[2] = 1;
    vecL[3] = 1;

    vecD[0] = 6;
    vecD[1] = 4;
    vecD[2] = 4;
    vecD[3] = 6;

    vecU[0] = 1;
    vecU[1] = 1;
    vecU[2] = 2;
    vecU[3] = 0;

    vecRHS[0] = 8;
    vecRHS[1] = 13;
    vecRHS[2] = 22;
    vecRHS[3] = 27;

    double* x = TridiagSolve(n, vecL, vecD, vecU, vecRHS);

    for(int i = 0; i < n; i++)
    {
        std::cout << x[i] << std::endl;
    }

    std::cout << std::endl;

    return 0;
}

