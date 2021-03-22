#include <iostream>
#include <cmath>
#include "JacobianApproximator.hpp"
#include "GeneralFunctions.hpp"

double* F(double* x);
double ComputeMatrix1Norm(double** x, double** y,int n);

int main(int argc, char* argv[])
{
    double* x;
    x = new double[2];

    x[0] = 2;
    x[1] = 2;

    double epsilon = 0.2;

    double** Jacobian = ComputeApproximateJacobian(x, F, 2, epsilon);

    std::cout << Jacobian[0][0] << " " << Jacobian[0][1] << std::endl;
    std::cout << Jacobian[1][0] << " " << Jacobian[1][1] << std::endl;

    std::cout << std::endl;

    DeleteMatrix(Jacobian, 2);

    double** ActualJacobian;
    ActualJacobian = new double*[2];
    for(int i=0; i<2; i++)
    {
        ActualJacobian[i] = new double[2];
    }

    ActualJacobian[0][0] = 4;
    ActualJacobian[0][1] = -1;
    ActualJacobian[1][0] = 4;
    ActualJacobian[1][1] = 4;

    std::cout << "epsilon" << "  " << "error" << std::endl;
    for(int i=0; i<6; i++)
    {
        std::cout << 0.2/pow(2,i) << "  " << ComputeMatrix1Norm(ActualJacobian, 
                ComputeApproximateJacobian(x, F, 2, 0.2/pow(2,i)), 2) << std::endl;
    }
    
    return 0;
}

double* F(double* x)
{
    double* y;
    y = new double[2];

    y[0] = pow(x[0],2)-x[1];
    y[1] = pow(x[0],2)+pow(x[1],2)-2;

    return y;
}

double ComputeMatrix1Norm(double** x, double** y,int n)
{
    double norm = 0;

    double** z;
    z = new double*[n];
    for(int i=0; i<n; i++)
    {
        z[i] = new double[n];
    }
    
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            z[i][j] = fabs(x[i][j]-y[i][j]);
        }
    }
    
    double* sumcols;
    sumcols = new double[n];

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            sumcols[i] += z[j][i];
        }
    }

    norm = sumcols[0];
    for(int i=0; i<n; i++)
    {
        if(sumcols[i]>norm)
        {
            norm = sumcols[i];
        }
    }
    return norm;
}
