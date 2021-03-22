#include "TridiagonalSolver.hpp"
#include "TridiagonalAssembler.hpp"
#include "GeneralFunctions.hpp"

double* SolveSystem(int noSubIntervals, double lengthInterval, double (*pFunction)(double x), double (*pDerivativeFunction)(double x))
{
    double* lower = AssembleLowerDiagonal(noSubIntervals);
    double* diagonal = AssembleDiagonal(noSubIntervals);
    double* upper = AssembleUpperDiagonal(noSubIntervals);
    
    double* rhs = AssembleRHS(noSubIntervals, lengthInterval, pFunction, pDerivativeFunction);

    double* partial_coefficients = SolveTridiagonalSystem(noSubIntervals+1, lower, diagonal, upper, rhs);

    double h;
    h = lengthInterval/(double)(noSubIntervals);
    double* grid;
    grid = new double[noSubIntervals+1];

    for(int i=0; i<noSubIntervals+1; i++)
    {
        grid[i] = i*h;

    }

    // Calculate the extra coefficients, c_-1 and c_n+1
    double cleft = partial_coefficients[1] - (h/3)*(*pDerivativeFunction)(grid[0]);

    double cright = partial_coefficients[noSubIntervals-1] + (h/3)*(*pDerivativeFunction)(grid[noSubIntervals]);

    double* coefficients;
    coefficients = new double[noSubIntervals+3];
    coefficients[0] = cleft;
    coefficients[noSubIntervals+2] = cright;
    
    for(int i=0; i<noSubIntervals+1; i++)
    {
        coefficients[i+1] = partial_coefficients[i];
    }

    return coefficients;
}

double* SolveSystem(int noSubIntervals, double lengthInterval, double (*pFunction)(double x))
{
    double* lower = AssembleLowerDiagonal(noSubIntervals);
    double* diagonal = AssembleDiagonal(noSubIntervals);
    double* upper = AssembleUpperDiagonal(noSubIntervals);
    
    double* rhs = AssembleRHSwApprox(noSubIntervals, lengthInterval, pFunction);

    double* partial_coefficients = SolveTridiagonalSystem(noSubIntervals+1, lower, diagonal, upper, rhs);

    double h;
    h = lengthInterval/(double)(noSubIntervals);
    double* grid;
    grid = new double[noSubIntervals+1];

    for(int i=0; i<noSubIntervals+1; i++)
    {
        grid[i] = i*h;

    }

    // Calculate the extra coefficients, c_-1 and c_n+1
    double cleft = partial_coefficients[1] - (h/3)*ApproximateDerivative(grid[0], pFunction);

    double cright = partial_coefficients[noSubIntervals-1] + (h/3)*ApproximateDerivative(grid[noSubIntervals], pFunction);

    double* coefficients;
    coefficients = new double[noSubIntervals+3];
    coefficients[0] = cleft;
    coefficients[noSubIntervals+2] = cright;
    
    for(int i=0; i<noSubIntervals+1; i++)
    {
        coefficients[i+1] = partial_coefficients[i];
    }

    return coefficients;
}