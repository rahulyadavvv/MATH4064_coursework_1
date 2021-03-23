#include <iostream>
#include <cmath>
#include "general.hpp"

//  functions created in this file
double* F(double stepLength, double* y, double* yN);
double** ExactInverseJacobian(double stepLength, double* px);
double* ComputeNewton(int iterations, double stepLength, double* initialX);
double* BackwardEulerStep(double stepLength, double* y, double* (*f)(int,
double, double*), int kIterations);
double* ApproximateBackwardEuler(double stepLength, int iterations, int 
kIterations, double* yInitial);

double* F(double stepLength, double* y, double* yN)
{
    double* solution = AllocateVector(2);
    solution[0] = y[0] - yN[0] - stepLength*(-3.0*pow(y[0], 2) + 2.0*y[1] + 1);
    solution[1] = y[1] - yN[1] - stepLength*(-3.0*pow(y[1], 2) + 2.0*y[0] + 1);
    return solution;
}

double** ExactInverseJacobian(double stepLength, double* px)
{
    //create 2x2 matrix
    double** solution = AllocateMatrix(2,2);

    //exact soltion for part inverse jacobian
    double det = 1.0/(1 + 6*stepLength*(px[0] + px[1]) +
    pow(stepLength, 2)*(36.0*px[0]*px[1] - 4));
    solution[0][0] = det*(1.0 + 6.0*stepLength*px[1]);
    solution[0][1] = det*(2.0*stepLength);
    solution[1][0] = det*(2.0*stepLength);
    solution[1][1] = det*(1.0 + 6.0*stepLength*px[0]);

    return solution;
}

//  
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
        double** JxInverse = ExactInverseJacobian(stepLength, x);
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
        PrintRowVector(2, x);

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
    double* newY = f(kIterations, stepLength, y);
    //  i don't think this part is right   
    //newY[0] = y[0] + stepLength*newY[0];
    //newY[1] = y[1] + stepLength*newY[1];
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
    yInitial[0] = 2.0;
    yInitial[1] = 1.0;
    double h = 1.0/12.0;
    double* approximation = ApproximateBackwardEuler(h, 1, 15, yInitial);

    PrintColVector(2, approximation);
    DeallocateVector(approximation);
    DeallocateVector(yInitial);
}