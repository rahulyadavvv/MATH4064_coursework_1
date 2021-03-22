#include "GeneralFunctions.hpp"

double* AssembleLowerDiagonal(int noSubIntervals)
{
    double* l = Vector(noSubIntervals+1);

    l[0] = 0;
    for(int i=1; i<noSubIntervals; i++)
    {
        l[i] = 1;
    }
    l[noSubIntervals] = 2;

    return l;
}

double* AssembleDiagonal(int noSubIntervals)
{
    double* d = Vector(noSubIntervals+1);

    for(int i=0; i<noSubIntervals+1; i++)
    {
        d[i] = 4;
    }

    return d;
}

double* AssembleUpperDiagonal(int noSubIntervals)
{
    double* u = Vector(noSubIntervals+1);

    u[0] = 2;
    for(int i=1; i<noSubIntervals; i++)
    {
        u[i] = 1;
    }
    u[noSubIntervals] = 0;

    return u;
}

double* AssembleRHS(int noSubIntervals, double intervalLength, double (*pFunction)(double x), double (*pDerivativeFunction)(double x))
{
    double* rhs = Vector(noSubIntervals+1);

    double subintervalLength;
    subintervalLength = intervalLength/(double)(noSubIntervals);
    double* grid;
    grid = new double[noSubIntervals+1];

    for(int i=0; i<noSubIntervals+1; i++)
    {
        grid[i] = i*subintervalLength;
    }

    rhs[0] = (*pFunction)(grid[0]) + (subintervalLength/3)*(*pDerivativeFunction)(grid[0]);

    for(int i=1; i<noSubIntervals; i++)
    {
        rhs[i] = (*pFunction)(grid[i]);
    }

    rhs[noSubIntervals] = (*pFunction)(grid[noSubIntervals]) - (subintervalLength/3)*(*pDerivativeFunction)(grid[noSubIntervals]);

    DeleteVector(grid);

    return rhs;
}

// Function to calculate derivative of a scalar real-valued function using a difference approximation
double ApproximateDerivative(double x, double (*pFunction)(double x))
{
    double derivative = 0;

    derivative = ((*pFunction)(x+0.0001)-(*pFunction)(x))/0.0001;

    return derivative;
}

double* AssembleRHSwApprox(int noSubIntervals, double intervalLength, double (*pFunction)(double x))
{
    double* rhs = Vector(noSubIntervals+1);

    double subintervalLength;
    subintervalLength = intervalLength/(double)(noSubIntervals);
    double* grid;
    grid = new double[noSubIntervals+1];

    for(int i=0; i<noSubIntervals+1; i++)
    {
        grid[i] = i*subintervalLength;
    }

    rhs[0] = (*pFunction)(grid[0]) + (subintervalLength/3)*ApproximateDerivative(grid[0], pFunction);

    for(int i=1; i<noSubIntervals; i++)
    {
        rhs[i] = (*pFunction)(grid[i]);
    }

    rhs[noSubIntervals] = (*pFunction)(grid[noSubIntervals]) - (subintervalLength/3)*ApproximateDerivative(grid[noSubIntervals], pFunction);

    DeleteVector(grid);

    return rhs;
}