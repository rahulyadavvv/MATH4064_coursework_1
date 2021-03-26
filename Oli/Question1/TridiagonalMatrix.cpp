#include "GeneralFunctions.hpp"

double* TridiagSolve(int n, double* vecL, double* vecD,
                            double* vecU, double* vecRHS)
{

    for (int i = 1; i < n; i++)

    {
        vecD[i] = vecD[i] - vecU[i-1]*vecL[i] / vecD[i-1];
        vecRHS[i] = vecRHS[i] - vecRHS[i-1] * vecL[i] / vecD[i-1];
    }

    double* x = Vector(n);

    for (int i = 1; i < n; i++)

    {
        x[i] = 0;
    }

    x[n-1] = vecRHS[n-1] / vecD[n-1];

    for (int i = n - 2; i > -1; i--)

    {
        x[i] = (vecRHS[i] - vecL[i] * x[i+1]) / vecD[i];
    }

    return x;
}

