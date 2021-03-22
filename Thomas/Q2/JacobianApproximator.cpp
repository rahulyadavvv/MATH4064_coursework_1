#include "GeneralFunctions.hpp"

double** ApproximateJacobian(double* vector, double (*pF1)(double x, double y), 
                    double (*pF2)(double x, double y), double epsilon, int n)
{
    double** Jacobian = Matrix(n);  // Set up jacobian matrix 
    
    Jacobian[0][0] = ((*pF1)(vector[0]+epsilon,vector[1])-(*pF1)(vector[0],vector[1]))/epsilon;
    Jacobian[1][0] = ((*pF2)(vector[0]+epsilon,vector[1])-(*pF2)(vector[0],vector[1]))/epsilon;
    Jacobian[0][1] = ((*pF1)(vector[0],vector[1]+epsilon)-(*pF1)(vector[0],vector[1]))/epsilon;
    Jacobian[1][1] = ((*pF2)(vector[0],vector[1]+epsilon)-(*pF2)(vector[0],vector[1]))/epsilon;

    return Jacobian;
}

double** ApproximateJacobianInverse(double* vector, double (*pF1)(double x, double y), 
                    double (*pF2)(double x, double y), double epsilon, int n)
{
    double** Jacobian = Matrix(n);  // Set up jacobian matrix

    Jacobian = ApproximateJacobian(vector, pF1, pF2, epsilon, n);

    double det;

    det = Jacobian[0][0]*Jacobian[1][1] - Jacobian[1][0]*Jacobian[0][1];

    double** inverseJacobian = Matrix(n);  // Set up inverse jacobian matrix

    inverseJacobian[0][0] = (1/det)*Jacobian[1][1];
    inverseJacobian[1][0] = -(1/det)*Jacobian[1][0];
    inverseJacobian[0][1] = -(1/det)*Jacobian[0][1];
    inverseJacobian[1][1] = (1/det)*Jacobian[0][0];

    DeleteMatrix(Jacobian, n);

    return inverseJacobian;
}

double** ComputeApproximateJacobian(double* x, double* (*pF)(double* x), int n, 
                                    double epsilon)
{
    double** Jacobian = Matrix(n);  // Set up jacobian matrix

    double* save = Vector(n);
    for(int i=0; i<n; i++)
    {
        save[i] = x[i];
    }

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            x[j] = x[j] + epsilon;
            Jacobian[i][j] = ((*pF)(x)[i]-(*pF)(save)[i])/epsilon;
            x[j] = save[j];
        }
    }

    DeleteVector(save);

    return Jacobian;
}

double** ComputeApproximateJacobianInverse(double* x, double* (*pF)(double* x), double epsilon, int n)
{
    double** Jacobian = ComputeApproximateJacobian(x, pF, n, epsilon);

    double det;

    det = Jacobian[0][0]*Jacobian[1][1] - Jacobian[1][0]*Jacobian[0][1];

    double** inverseJacobian = Matrix(n);  // Set up inverse jacobian matrix

    inverseJacobian[0][0] = (1/det)*Jacobian[1][1];
    inverseJacobian[1][0] = -(1/det)*Jacobian[1][0];
    inverseJacobian[0][1] = -(1/det)*Jacobian[0][1];
    inverseJacobian[1][1] = (1/det)*Jacobian[0][0];

    DeleteMatrix(Jacobian, n);

    return inverseJacobian;
}

void ComputeApproximateJacobianInverse(double** matrix, double* x, double* (*pF)(double* x), double epsilon, int n)
{
    double** Jacobian = ComputeApproximateJacobian(x, pF, n, epsilon);

    double det;

    det = Jacobian[0][0]*Jacobian[1][1] - Jacobian[1][0]*Jacobian[0][1];

    matrix[0][0] = (1/det)*Jacobian[1][1];
    matrix[1][0] = -(1/det)*Jacobian[1][0];
    matrix[0][1] = -(1/det)*Jacobian[0][1];
    matrix[1][1] = (1/det)*Jacobian[0][0];

    DeleteMatrix(Jacobian, n);
}
