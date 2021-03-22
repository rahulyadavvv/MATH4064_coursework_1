#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert> 
#include "EulerMethod.hpp"
#include "GeneralFunctions.hpp"

double* f(double* x);
void WriteToDatFile(int n, double* vector1, double* vector2, std::string fileName);


int main(int argc, char* argv[])
{
    double* initial;
    initial = new double[2];
    initial[0] = 2;
    initial[1] = 1;

    double h = 1.0/12.0;

    int n_max = 24;

    double* truesol;
    truesol = new double[2];
    truesol[0] = 1;
    truesol[1] = 1;

    double* n = Vector(n_max+1);
    double* error = Vector(n_max+1);

    std::cout << "n" << "  " << "y_1,n" << "  " << "y_2,n" << "  " 
    << "||y_n -y_*||" << std::endl; 

    for(int i=0; i<=n_max; i++)
    {
        n[i] = i;
        std::cout << i << "  " << ComputeEulerIteration(i+1, initial, f, h)[0] 
        << "  " << ComputeEulerIteration(i+1, initial, f, h)[1] << "  " 
        << Compute2Norm(ComputeEulerIteration(i+1, initial, f, h), truesol, 2) 
        << std::endl;
        error[i] = Compute2Norm(ComputeEulerIteration(i+1, initial, f, h), truesol, 2);
    }

    WriteToDatFile(n_max+1, n, error, "errorvsn.dat");

    return 0;
}

double* f(double* x)
{
    double* y;
    y = new double[2];

    y[0] = -3*pow(x[0], 2)+ 2*x[1] + 1;
    y[1] = -3*pow(x[1], 2)+ 2*x[0] + 1;

    return y;
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
