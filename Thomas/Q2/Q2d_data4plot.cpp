#include <iostream>
#include <fstream>
#include <cassert>
#include "NewtonApproximateJacobian.hpp"
#include "GeneralQ2.hpp"

void WriteToDatFile(int n, double* vector1, double* vector2, std::string fileName);

int main(int argc, char* argv[])
{
    int n = 2;
    int noIterations = 15;

    double* initialGuess = AllocateVector(n);
    initialGuess[0] = 2;
    initialGuess[1] = 2;

    double* actualVector = AllocateVector(n);
    actualVector[0] = 1;
    actualVector[1] = 1;

    double epsilon = 0.2;

    double error;
    double* errorVector = AllocateVector(n);
    for(int i=0; i<n; i++)
    {
        errorVector[i] = actualVector[i] - initialGuess[i];
    }
    error = ComputeOneNorm(n, errorVector);

    double* k = AllocateVector(noIterations);
    double* errorValues = AllocateVector(noIterations);

    errorValues[0] = error;

    for(int i=0; i<noIterations; i++)
    {
        k[i] = i;
    }

    double* x = NewtonMethodwApproximateJacobian(initialGuess, F, n, epsilon, 1);
    for(int j=1; j<noIterations; j++)
    {
        for(int i=0; i<n; i++)
        {
            errorVector[i] = actualVector[i] - x[i];
        }
        error = ComputeOneNorm(n, errorVector);
        errorValues[j] = error;
        x = NewtonMethodwApproximateJacobian(x, F, n, epsilon, 1);
    }
    
    DeallocateVector(errorVector);
    DeallocateVector(x);
    DeallocateVector(initialGuess);
    DeallocateVector(actualVector);

    WriteToDatFile(noIterations, k, errorValues, "errorvsk.dat");

    return 0;
}

void WriteToDatFile(int n, double* vector1, double* vector2, std::string fileName)
{    
    std::ofstream writeFile; // Define output stream 
    writeFile.open(fileName); // Open file 
    assert(writeFile.is_open()); // Check file is open 
    for (int i=0; i<n; i++) 
    {
        writeFile << vector1[i] << " " << vector2[i] << "\n"; 
    }
    writeFile.close(); // Close file
}
