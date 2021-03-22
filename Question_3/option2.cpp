#include <iostream>
#include <cmath>
#include "general.hpp"

//  functions created in this file
double* F(double stepLength, double* y, double* yN);
double** Jacobian(double stepLength, double* px);
double* ComputeNewton(int iterations, double stepLength, double* initialX);
double* BackwardEulerStep(double stepLength, double* y, double* (*f)(int,
double, double*), int kIterations);

double* F(double stepLength, double* y, double* yN)
{
    double* solution = AllocateVector(2);
    solution[0] = y[0] - yN[0] - (-3.0*pow(y[0], 2), 2.0*y[1] + 1);
    solution[1] = y[1] - yN[1] - (-3.0*pow(y[1], 2), 2.0*y[0] + 1);
    return solution;
}

double** Jacobian(double stepLength, double* px)
{
    //create 2x2 matrix
    double** solution = AllocateMatrix(2,2);
    //exact soltion for part a jacobian
    solution[0][0] = 1.0 + 6.0*stepLength*px[0];
    solution[0][1] = -2.0*stepLength;
    solution[1][0] = -2.0*stepLength;
    solution[1][1] = 1.0 + 6.0*stepLength*px[1];

    return solution;
}

//  edit to use previous x thats needed in F
double* ComputeNewton(int iterations, double stepLength, double* initialX)
{
    double* x = AllocateVector(2);
    double* xPrevious = AllocateVector(2);
    x[0] = initialX[0];
    x[1] = initialX[1];
    xPrevious[0] = x[0];
    xPrevious[1] = x[1];
    for (int i = 0; i < iterations; i++)
    {
        //  calculate J(x_i)^(-1)
        double** JxInverse = Inverse2x2(Jacobian(stepLength, x));
        //  calculate value of F at x_i
        double* Fx = F(stepLength, x, xPrevious);
        //  calculate J(x_i)^(-1)*F(x_i)
        double* product = MultiplyMatrixVector(2, 2, JxInverse, Fx);

        //  store previous value for F in next iteration
        xPrevious[0] = x[0];
        xPrevious[1] = x[1];

        //  calculate x_(i+1) = x_i  - J(x_i)^(-1)*F(x_i)
        x[0] -= product[0];
        x[1] -= product[1];

        DeallocateVector(product);
        DeallocateVector(Fx);
        DeallocateMatrix(2, JxInverse);
    }
    DeallocateVector(xPrevious);
    return x;
}

double* BackwardEulerStep(double stepLength, double* y, double* (*f)(int,
double, double*), int kIterations)
{
    double* newY = AllocateVector(2);
    for (int i = 0; i < 2; i++)
    {
        newY[i] = y[i] + stepLength*f(kIterations, stepLength, y)[i];
    }
    return newY;
}

double* ApproximateBackwardEuler(double stepLength, int iterations, int 
kIterations, double* yInitial)
{
    double* solution = AllocateVector(2);
    solution[0] = yInitial[0];
    solution[1] = yInitial[1];

    for (int i = 0; i < iterations; i++)
    {
        solution = BackwardEulerStep(stepLength, solution, ComputeNewton,
        kIterations);
    }
    return solution;
}

int main(int argc, char* argv[]){

    double* yInitial = AllocateVector(2);
    double* approximation = ApproximateBackwardEuler(1.0/2.0, 1, 5, yInitial);
    DeallocateVector(approximation);
    DeallocateVector(yInitial);
}