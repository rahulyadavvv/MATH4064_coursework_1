#include <iostream>
#include "general.hpp"

//  calculate exact jacobian of F at value x
double** Jacobian(double* px)
{
    //create 2x2 matrix
    double** solution = AllocateMatrix(2,2);
    //exact soltion for part a jacobian
    solution[0][0] = 2*px[0];
    solution[0][1] = -1.0;
    solution[1][0] = 2*px[0];
    solution[1][1] = 2*px[1];

    return solution;
}


//  inverse of a 2x2 matrix
double** Inverse2x2(double** pMatrix)
{
    //  calculate determininant
    double det = 1.0/(pMatrix[0][0]*pMatrix[1][1] - pMatrix[0][1]*pMatrix[1][0]);
    //  swap a and d then multiply by det
    double a = pMatrix[0][0];
    pMatrix[0][0] = det*pMatrix[1][1];
    pMatrix[1][1] = det*a;
    //  multiply b and c by -det
    pMatrix[0][1] = -1.0*det*pMatrix[0][1];
    pMatrix[1][0] = -1.0*det*pMatrix[1][0];

    return pMatrix;
}

//  compute function at x
double* F(double* px)
{
    double* solution = AllocateVector(2);
    solution[0] = px[0]*px[0] - px[1];
    solution[1] = px[0]*px[0] + px[1]*px[1] - 2;
    return solution;
}

double* MultiplyMatrixVector(int noRows, int noCols, double** pMatrix, double* 
pVector)
{
    double* pResult = AllocateVector(noRows);
    for (int i = 0; i < noRows; i++)
        {
            pResult[i] = 0.0;
            for (int j = 0; j < noCols; j++)
            {
                pResult[i] += pMatrix[i][j]*pVector[j];
            }
        }
    return pResult;
}

int main(int argc, char* argv[]){

    double** x = AllocateMatrix(6, 2);
    double* xReal = AllocateVector(2);

    x[0][0] = 2.0;
    x[0][1] = 2.0;

    for (int i = 0; i < 5; i++) //stop before last row at x_4
    {
        double** JxInverse = Inverse2x2(Jacobian(x[i]));
        double* Fx = F(x[i]);
        double* product = MultiplyMatrixVector(2, 2, JxInverse, Fx);

        x[i + 1][0] = x[i][0] - product[0];
        x[i + 1][1] = x[i][1] - product[1];

        //need to still test

        DeallocateVector(product);
        DeallocateVector(Fx);
        DeallocateMatrix(2, JxInverse);
    }
    
    PrintMatrix(6, 2, x);
    DeallocateVector(xReal);
    DeallocateMatrix(5, x);
}