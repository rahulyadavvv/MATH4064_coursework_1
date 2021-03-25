#include <iostream>
#include "general.hpp"

double* InitialiseMainDiagonal(int length)
{
    double* result = AllocateVector(length);
    for (int i = 0; i < length; i++)
    {
        result[i] = 4.0;
    }
    return result;
}

double* InitialiseLowerDiagonal(int length)
{
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
    double* result = AllocateVector(length);
    result[0] = 2.0;
    for (int i = 1; i < length - 1; i++)
    {
        
        result[i] = 1.0;
    }
    result[length - 1] = 0.0;
    return result;
}

int main(int argc, char* argv[]){
    //  get user to input length and n points
    int n = 6;
    double length = 3.0;

    //  calculate and store step size in h 
    double h = length/((double)n);

    /* allocate memory for diagonals and assign d = main diagonal, u = upper
    diagonal, l = lower diagonal     */
    double* d = InitialiseMainDiagonal(n);
    double* u = InitialiseUpperDiagonal(n);
    double* l = InitialiseLowerDiagonal(n);
    
    //  print vectors for upper diagonal, main diagonal and lower diagonal
    std::cout << "upper diagonal: ";
    PrintRowVector(n, u);
    std::cout << "main diagonal:  ";
    PrintRowVector(n, d);
    std::cout << "lower diagonal: ";
    PrintRowVector(n, l);

    //  deallocate memory for upper diagonal, main diagonal and lower diagonal
    DeallocateVector(l);
    DeallocateVector(u);
    DeallocateVector(d);
}