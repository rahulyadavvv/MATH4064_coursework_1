#include <iostream>
#include "NewtonApproximateJacobian.hpp"
#include "GeneralQ2.hpp"

int main(int argc, char* argv[])
{
    int n = 2;
    int noIterations = 10;

    double* initialGuess = AllocateVector(n);
    initialGuess[0] = 2;
    initialGuess[1] = 2;

    double* actualVector = AllocateVector(n);
    actualVector[0] = 1;
    actualVector[1] = 1;

    double epsilon = 0.2;

    std::cout << "k" << "  " << "error in one norm" << std::endl;
    double error;
    double* errorVector = AllocateVector(n);
    for(int i=0; i<n; i++)
    {
        errorVector[i] = actualVector[i] - initialGuess[i];
    }
    error = ComputeOneNorm(n, errorVector);
    
    std::cout << "0"  << "  " << error << std::endl; 

    double* x = NewtonMethodwApproximateJacobian(initialGuess, F, n, epsilon, 1);

    int k = 1;

    while(error>10e-12)
    {
        std::cout << k << "  ";
        for(int i=0; i<n; i++)
        {
            errorVector[i] = actualVector[i] - x[i];
        }
        error = ComputeOneNorm(n, errorVector);
        std::cout << error << std::endl;
        x = NewtonMethodwApproximateJacobian(x, F, n, epsilon, 1);
        k++;
    }
    DeallocateVector(errorVector);
    DeallocateVector(x);
    DeallocateVector(initialGuess);
    DeallocateVector(actualVector);

    return 0;
}


