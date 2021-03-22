#include <iostream>
#include "JacobianApproximator.hpp"
#include "GeneralQ2.hpp"

double* NewtonMethodwApproximateJacobian(double* initialVector, double* (*pF)(double* x), int n, double epsilon, int noIterations)
{
    double* x;
    x = new double[n];

    for(int i=0; i<n; i++)
    {
        x[i] = initialVector[i];
    }

    double** JxInverse = ComputeApproximateJacobianInverse(x, pF, epsilon, n);
    double* Fx = pF(x);
    double* product = MultiplyMatrixVector(n, n, JxInverse, Fx);

    for(int j=0; j<n; j++)
    {
        x[j] = x[j] - product[j];
    }

    for(int i=1; i<noIterations; i++)
    {
        ComputeApproximateJacobianInverse(JxInverse, x, pF, epsilon, n);
        F(x, Fx);
        MultiplyMatrixVector(n, n, product, JxInverse, Fx);

        for(int j=0; j<n; j++)
        {
            x[j] = x[j] - product[j];
        }
    }
    DeallocateVector(product);
    DeallocateVector(Fx);
    DeallocateMatrix(2, JxInverse);

    return x;
}
