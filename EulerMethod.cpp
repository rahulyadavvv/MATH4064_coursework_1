#include <cmath>

double* SolveODESystemEuler(double* initial, int maxIter, double* (*pFunction)(double* x), double h)
{
    double* y_new;
    y_new = new double[2];
    double* y_old;
    y_old = new double[2];
    y_old[0] = initial[0];
    y_old[1] = initial[1];

    for(int i=0; i<maxIter; i++)
    {
        y_new[0] = y_old[0] + h*(*pFunction)(y_old)[0];
        y_new[1] = y_old[1] + h*(*pFunction)(y_old)[1];
        y_old[0] = y_new[0];
        y_old[1] = y_new[1];
    }

    return y_new;
}

double* ComputeEulerIteration(int iter, double* initial, double* (*pFunction)(double* x), double h)
{
    double* y_new;
    y_new = new double[2];
    double* y_old;
    y_old = new double[2];
    y_old[0] = initial[0];
    y_old[1] = initial[1];

    for(int i=0; i<iter; i++)
    {
        y_new[0] = y_old[0] + h*(*pFunction)(y_old)[0];
        y_new[1] = y_old[1] + h*(*pFunction)(y_old)[1];
        y_old[0] = y_new[0];
        y_old[1] = y_new[1];
    }

    return y_new;
}

double Compute2Norm(double* x, double* y, int n)
{
    double norm = 0;

    for(int i=0; i<n; i++)
    {
        norm += pow(x[i]-y[i],2);
    }
    return sqrt(norm);
}