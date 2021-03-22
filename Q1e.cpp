#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cassert>
#include "SystemSolver.hpp"
#include "GeneralFunctions.hpp"

double Q1d(double x);
double Q1dDerivative(double x);
double B(double x);
double B_local(double x, int i, double h, double* grid);
double spline(double x, double* c, int n, double h);
void WriteToDatFile(int n, double* mesh, double* solution_vector, std::string fileName);

int main(int argc, char* argv[])
{
    double l = 3;
    double h = 0.01;
    int n = l/h;

    double* c = SolveSystem(n, l, Q1d, Q1dDerivative);

    double x = 1.0/3.0;

    std::cout << "h" << " " << "q_h(1/3)" << "  " << "|f(1/3)-q_h(1/3)|" << std::endl;
    for(int i=1; i<20; i++)
    {
        double h1 = 1.0/pow(2,i);
        double n1 = l/h1;
        c = SolveSystem(n1, l, Q1d, Q1dDerivative);
        double q = spline(x, c, n1, h1);
        std::cout << h1 << " " << q << " " << fabs(q-Q1d(1.0/3.0)) << std::endl;
    }

    int no_values = 40;

    double* hvalues = Vector(no_values);
    double* errorvalues = Vector(no_values); 

    for(int i=0; i<no_values; i++)
    {
        hvalues[i] = 1/pow(2,i+1);
        double n1 = l/hvalues[i];
        c = SolveSystem(n1, l, Q1d, Q1dDerivative);
        double q = spline(x, c, n1, hvalues[i]);
        errorvalues[i] = fabs(q-Q1d(1.0/3.0));
    }

    WriteToDatFile(no_values, hvalues, errorvalues, "hvserror.dat");

    return 0;
}

double Q1d(double x)
{
    return exp(x/2)*cos((1.0/2.0)*M_PI*x);
}

double Q1dDerivative(double x)
{
    return (1.0/2.0)*exp(x/2.0)*cos((1.0/2.0)*M_PI*x) - (M_PI/2.0)*exp(x/2.0)*sin((1.0/2.0)*M_PI*x);
}


double B(double x)
{
    if (x<=-2)
    {
        return 0;
    }
    else if (-2<x && x<=-1)
    {
        return pow(x+2,3);
    }
    else if (-1<x && x<=0)
    {
        return 1+3*(x+1)+3*pow(x+1,2)-3*pow(x+1,3);
    }
    else if (0<x && x<=1)
    {
        return 1+3*(1-x)+3*pow(1-x,2)-3*pow(1-x,3);
    }
    else if (1<x && x<=2)
    {
        return pow(2-x,3);
    }
    else
    {
        return 0;
    }
    
}

double B_local(double x, int i, double h, double* grid)
{
    return B((x-grid[i])/h);
}

double spline(double x, double* c, int n, double h)
{
    double* grid;
    grid = new double[n+3];

    grid[0] = -h;
    for(int i=1; i<n+3; i++)
    {
        grid[i] = i*h;
    }

    int k = floor(x/h) + 1;

    double s = 0;
    for(int i=k-2; i<k+2; i++)
    {
        s += c[i]*B_local(x, i, h, grid);
    }

    return s;
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