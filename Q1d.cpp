#include <iostream>
#include <cmath>
#include "SystemSolver.hpp"

double Q1d(double x);
double Q1dDerivative(double x);

int main(int argc, char* argv[])
{
    double l = 3;
    double h = 0.5;
    int n = 6;

    double* c = SolveSystem(n, l, Q1d, Q1dDerivative);

    for(int i=0; i<n+3; i++)
    {
        std::cout << c[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << std::endl;

    double* c_approx = SolveSystem(n, l, Q1d);

    for(int i=0; i<n+3; i++)
    {
        std::cout << c_approx[i] << std::endl;
    }
    std::cout << std::endl;

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