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

double ComputeOneNorm(int noRows, double* x)
{
    double result = 0.0;
    for (int i = 0; i < noRows; i++)
    {
        result += abs(x[i]);
    }
    return result;
}

//  approximate solution using newton-raphson method (only vector length 2)
double** ComputeNewton(int iterations, double* initialX)
{
    double** x = AllocateMatrix(iterations, 2);
    x[0][0] = initialX[0];
    x[0][1] = initialX[1];
    for (int i = 0; i < iterations - 1; i++) //stop before last row at x_4
    {
        //  calculate J(x_i)^(-1)
        double** JxInverse = Inverse2x2(Jacobian(x[i]));
        //  calculate value of F at x_i
        double* Fx = F(x[i]);
        //  calculate J(x_i)^(-1)*F(x_i)
        double* product = MultiplyMatrixVector(2, 2, JxInverse, Fx);

        //  calculate x_(i+1) = x_i  - J(x_i)^(-1)*F(x_i)
        x[i + 1][0] = x[i][0] - product[0];
        x[i + 1][1] = x[i][1] - product[1];

        DeallocateVector(product);
        DeallocateVector(Fx);
        DeallocateMatrix(2, JxInverse);
    }
    return x;
}

int main(int argc, char* argv[]){

    //  set iterations to 6 used in ComputeNewton() to produce x_{0, ... , 5}
    const int iterations = 6;

    //  initialise initial value of x used in ComputeNewton();
    double* xInitial = AllocateVector(2);
    xInitial[0] = 2.0;
    xInitial[1] = 2.0;

    //  store value for x when F is 0
    double* xReal = AllocateVector(2);
    xReal[0] = 1.0;
    xReal[1] = 1.0;
    
    //  approximate solution using newton-raphson method 
    double** x = ComputeNewton(iterations, xInitial);

    //  output x_k for k = 0, 1, 2, 3, 4, 5
    std::cout << "x_k for k = 0, 1, 2, 3, 4, 5:" << std::endl;
    PrintMatrix(iterations, 2, x);

    double* error = AllocateVector(iterations);
    for (int i = 0; i < iterations; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            x[i][j] = xReal[j] - x[i][j];
        }
        error[i] = ComputeOneNorm(2, x[i]);
    }
    
    //  output error for k = 0, 1, 2, 3, 4, 5
    std::cout << "error for k = 0, 1, 2, 3, 4, 5:" << std::endl;
    PrintColVector(iterations, error);

    //  deallovate all vectors and matrices 
    DeallocateVector(error);
    DeallocateMatrix(iterations, x);
    DeallocateVector(xReal);
    DeallocateVector(xInitial);
}