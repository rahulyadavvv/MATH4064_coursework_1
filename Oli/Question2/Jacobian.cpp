#include "GeneralFunctions.hpp"

double** JacobianApprox(double* vec, double (*pF1)(double x, double y),
                    double (*pF2)(double x, double y), double epsilon, int n)


{
    double** Jacobian = Matrix(n);

    Jacobian[0][0] = ((*pF1)(vec[0] + epsilon,
                        vec[1]) - (*pF1)(vec[0],
                        vec[1])) / epsilon;
    Jacobian[1][0] = ((*pF2)(vec[0] + epsilon,
                        vec[1]) - (*pF2)(vec[0],
                        vec[1])) / epsilon;
    Jacobian[0][1] = ((*pF1)(vec[0],
                        vec[1] + epsilon) - (*pF1)(vec[0],
                        vec[1])) / epsilon;
    Jacobian[1][1] = ((*pF2)(vec[0],
                        vec[1] + epsilon)-(*pF2)(vec[0],
                        vec[1])) / epsilon;

    return Jacobian;
}

double** JacobianInverseApprox(double* vec, double (*pF1)(double x, double y),
                         double (*pF2)(double x, double y), double epsilon, int n)

{

    double** Jacobian = Matrix(n);

    Jacobian = JacobianApprox(vec, pF1, pF2, epsilon, n);

    double Determinant;

    Determinant = Jacobian[0][0]*Jacobian[1][1] - Jacobian[1][0]*Jacobian[0][1];

    double** JacobianInverse = Matrix(n);

    JacobianInverse[0][0] = (1/Determinant) * Jacobian[1][1];
    JacobianInverse[1][0] = -(1/Determinant) * Jacobian[1][0];
    JacobianInverse[0][1] = -(1/Determinant) * Jacobian[0][1];
    JacobianInverse[1][1] = (1/Determinant) * Jacobian[0][0];

    DeleteMatrix(Jacobian, n);

    return JacobianInverse;

}

double** ComputeApprox(double* x, double* (*pF)(double* x),
                              int n, double epsilon)
{

    double** Jacobian = Matrix(n);

    double* vec = Vector(n);

    for(int i = 0; i < n; i++)

    {
        vec[i] = x[i];
    }

    for (int i = 0; i < n; i++)
    {

        for (int j = 0; j < n; j++)
        {
            x[j] = x[j] + epsilon;
            Jacobian[i][j] = ((*pF)(x)[i]-(*pF)(vec)[i]) / epsilon;
            x[j] = vec[j];
        }
    }

    DeleteVector(vec);

    return Jacobian;

}

double** ComputeInverseApprox(double* x, double* (*pF)(double* x), int n, double epsilon)
{
    double** Jacobian = ComputeApprox(x, pF, n, epsilon);

    double Determinant;

    Determinant = Jacobian[0][0]*Jacobian[1][1] - Jacobian[1][0]*Jacobian[0][1];

    double** inverseJacobian = Matrix(n);

    inverseJacobian[0][0] = (1/Determinant) * Jacobian[1][1];
    inverseJacobian[1][0] = -(1/Determinant) * Jacobian[1][0];
    inverseJacobian[0][1] = -(1/Determinant) * Jacobian[0][1];
    inverseJacobian[1][1] = (1/Determinant) * Jacobian[0][0];

    DeleteMatrix(Jacobian, n);

    return inverseJacobian;
}

void ComputeInverseApprox(double** matrix, double* x, double* (*pF)(double* x), double epsilon, int n)
{
    double** Jacobian = ComputeApprox(x, pF, n, epsilon);

    double Determinant;

    Determinant = Jacobian[0][0]*Jacobian[1][1] - Jacobian[1][0]*Jacobian[0][1];

    matrix[0][0] = (1/Determinant) * Jacobian[1][1];
    matrix[1][0] = -(1/Determinant) * Jacobian[1][0];
    matrix[0][1] = -(1/Determinant) * Jacobian[0][1];
    matrix[1][1] = (1/Determinant) * Jacobian[0][0];

    DeleteMatrix(Jacobian, n);
}


