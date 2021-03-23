#include <iostream>
#include <iomanip>
#include <cmath>
#include "general.hpp"

//  functions created in this file
void F(double stepLength, double* y, double* yN, double* solution);
void ExactInverseJacobian(double stepLength, double* pX, double** solution);
void ComputeNewton(int iterations, double stepLength, double* pX, 
bool produceOutput);
double* ApproximateBackwardEuler(double stepLength, int iterations, int 
kIterations, double* yInitial, bool output);

//  function F, equation 5 in Q3
void F(double stepLength, double* y, double* yN, double* solution)
{
    solution[0] = y[0] - yN[0] - stepLength*(-3.0*pow(y[0], 2) + 2.0*y[1] + 1);
    solution[1] = y[1] - yN[1] - stepLength*(-3.0*pow(y[1], 2) + 2.0*y[0] + 1);
}

//  the inverse of the jacobian, equation 7 in Q3
void ExactInverseJacobian(double stepLength, double* pX, double** solution)
{
    //exact soltion for part inverse jacobian
    double det = 1.0/(1 + 6*stepLength*(pX[0] + pX[1]) +
    pow(stepLength, 2)*(36.0*pX[0]*pX[1] - 4));
    solution[0][0] = det*(1.0 + 6.0*stepLength*pX[1]);
    solution[0][1] = det*(2.0*stepLength);
    solution[1][0] = det*(2.0*stepLength);
    solution[1][1] = det*(1.0 + 6.0*stepLength*pX[0]);
}

/*  Computes the newton-raphoson method for for given n, equation 6 in Q3
    also produces the table output  */
void ComputeNewton(int iterations, double stepLength, double* pX, 
bool produceOutput)
{
    double* xPrevious = AllocateVector(2);
    xPrevious[0] = pX[0];
    xPrevious[1] = pX[1];

    //  table for k iterations output - header
    if (produceOutput)
    {
        std::cout << std::left << std::setw(20) << std::setfill(' ') << "k" <<
        std::left << std::setw(20) << std::setfill(' ') << "y_{1, 1}^(k)" << 
        std::left << std::setw(20) << std::setfill(' ') << "y_{2, 1}^(k)" << 
        std::left << std::setw(20) << std::setfill(' ') << "F_1 (y_{1}^(k))" << 
        std::left << std::setw(20) << std::setfill(' ') << "F_2 (y_{1}^(k))" 
        << std::endl;
    }
    
    double** JxInverse = AllocateMatrix(2, 2);
    double* Fx = AllocateVector(2);

    for (int i = 0; i < iterations; i++)
    {
        //  calculate J(x_i)^(-1)
        ExactInverseJacobian(stepLength, pX, JxInverse);
        //  calculate value of F at x_i
        F(stepLength, pX, xPrevious, Fx);
        //  calculate J(x_i)^(-1)*F(x_i)
        double* product = MultiplyMatrixVector(2, 2, JxInverse, Fx);

        //  store previous value for F in next iteration
        xPrevious[0] = pX[0];
        xPrevious[1] = pX[1];
        
        //  calculate x_(i+1) = x_i  - J(x_i)^(-1)*F(x_i)
        pX[0] -= product[0];
        pX[1] -= product[1];

        //  table for k iterations output - rows
        if (produceOutput)
        {
            std::cout << std::left << std::setw(20) << std::setfill(' ') << i <<
            std::left << std::setw(20) << std::setfill(' ') << pX[0] << 
            std::left << std::setw(20) << std::setfill(' ') << pX[1] << 
            std::left << std::setw(20) << std::setfill(' ') << Fx[0] << 
            std::left << std::setw(20) << std::setfill(' ') << Fx[1] 
            << std::endl;
        }

        DeallocateVector(product);
    }
    DeallocateVector(Fx);
    DeallocateMatrix(2, JxInverse);
    DeallocateVector(xPrevious);
}

//  the Backward Euler method
double* ApproximateBackwardEuler(double stepLength, int iterations, int 
kIterations, double* yInitial, bool output)
{
    double* solution = AllocateVector(2);
    solution[0] = yInitial[0];
    solution[1] = yInitial[1];

    if (!output)
    {
        std::cout << std::left << std::setw(20) << std::setfill(' ') << "n" <<
        std::left << std::setw(20) << std::setfill(' ') << "y_{1, n}" << 
        std::left << std::setw(20) << std::setfill(' ') << "y_{2, n}" << 
        std::left << std::setw(20) << std::setfill(' ') << "||y_n - y*||_2" 
        << std::endl;
    }
    
    for (int i = 0; i < iterations; i++)
    {
        ComputeNewton(kIterations, stepLength, solution, output);
    }
    return solution;
}

int main(int argc, char* argv[]){

    double* yInitial = AllocateVector(2);
    yInitial[0] = 2.0;
    yInitial[1] = 1.0;
    double h = 1.0/12.0;
    //  output == true -->  output newton for each k
    bool output = true;
    double* approximation = ApproximateBackwardEuler(h, 1, 6, yInitial, output);
    /*   code for table output contained in Computenewton   */
    DeallocateVector(approximation);

    //  output == false --> output error for n
    //approximation = ApproximateBackwardEuler(h, 25, 5, yInitial, false);
    //DeallocateVector(approximation);

    DeallocateVector(yInitial);
}