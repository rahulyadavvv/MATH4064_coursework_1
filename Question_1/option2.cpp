#include <iostream>
#include "general.hpp"

//  functions created in this file
double* InitialiseMainDiagonal(int length);
double* InitialiseLowerDiagonal(int length);
double* InitialiseUpperDiagonal(int length);
double* InitialiseRHS(int noSteps, double intervalLength, 
double (*pFunction)(double), double (*pFunctionDerivative)(double));


double* InitialiseMainDiagonal(int length)
{
    length++;
    double* result = AllocateVector(length);
    for (int i = 0; i < length; i++)
    {
        result[i] = 4.0;
    }
    return result;
}

double* InitialiseLowerDiagonal(int length)
{
    length++;
    double* result = AllocateVector(length);
    result[0] = 0.0;
    for (int i = 1; i < length - 1; i++)
    {
        
        result[i] = 1.0;
    }
    result[length - 1] = 2.0;
    return result;
}

double* InitialiseUpperDiagonal(int length)
{
    length++;
    double* result = AllocateVector(length);
    result[0] = 2.0;
    for (int i = 1; i < length - 1; i++)
    {
        
        result[i] = 1.0;
    }
    result[length - 1] = 0.0;
    return result;
}

double* InitialiseRHS(int noSteps, double intervalLength, 
double (*pFunction)(double), double (*pFunctionDerivative)(double))
{
    double h = intervalLength/((double)noSteps);

    double* result = AllocateVector(noSteps + 1);
    result[0] = pFunction(0.0) + (h/3.0)*pFunctionDerivative(0.0);
    for (int i = 1; i < noSteps; i++)
    {
        result[i] = pFunction((double)(i)*h);
    }
    result[noSteps] = pFunction((double)(noSteps)*h) - 
    (h/3.0)*pFunctionDerivative((double)(noSteps)*h);

    return result;
}

//  used for testing InitialiseRHS()
double F(double x){return x;}
double Fprime(double x){return 1.0;}

int main(int argc, char* argv[]){
    //  get user to input length and n points
    int n = 6;
    double length = 3.0;

    /* allocate memory for diagonals and assign d = main diagonal, u = upper
    diagonal, l = lower diagonal     */
    double* d = InitialiseMainDiagonal(n);
    double* u = InitialiseUpperDiagonal(n);
    double* l = InitialiseLowerDiagonal(n);
    //  using test function for InitialiseRHS()
    double* rhs = InitialiseRHS(n, length, F, Fprime);
    
    //  print vectors for upper diagonal, main diagonal and lower diagonal
    n++;
    std::cout << "upper diagonal: ";
    PrintRowVector(n, u);
    std::cout << "main diagonal:  ";
    PrintRowVector(n, d);
    std::cout << "lower diagonal: ";
    PrintRowVector(n, l);
    std::cout << "right-hand side:";
    PrintRowVector(n, rhs);
    std::cout << "note: not RHS from Q1, it's the RHS generated using f(x) = x"
    << std::endl;

    //  deallocate memory for upper diagonal, main diagonal and lower diagonal
    DeallocateVector(rhs);
    DeallocateVector(l);
    DeallocateVector(u);
    DeallocateVector(d);
}