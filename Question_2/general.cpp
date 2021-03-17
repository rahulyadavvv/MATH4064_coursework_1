#include <iostream>
#include <cmath>

/*  AllocateVector allocates memory for vector of dimensions (n x 1), where
    (integer) noRows = n. */
double* AllocateVector(int noRows)
{
    double* vector;
    vector = new double [noRows];
    return vector;
}

/*  DeallocateVector deallocates memory allocated for vector. */
void DeallocateVector(double* vector)
{
    delete[] vector;
}

/*  AllocateMatrix allocates memory for matrix of dimensions (n x m), where
    (integer) noRows = n and (integer) noCols = m. */
double** AllocateMatrix(int noRows, int noCols)
{
    double** matrix;
    matrix = new double* [noRows];
    for (int i = 0; i < noRows; i++)
    {
        matrix[i] = new double [noCols];
    }
    return matrix;
}

/*  DeallocateMatrix deallocates memory allocated for matrix of dimensions
    (n x m), where (integer) noRows = n. */
void DeallocateMatrix(int noRows, double** matrix)
{
    for (int i = 0; i < noRows; i++)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}

/*  Print matrix outputs matrix of dimensions (n x m), where (integer) noRows = n 
    and (integer) noCols = m. */
void PrintMatrix(int noRows, int noCols, double** matrix)
{
    for (int i = 0; i < noRows; i++)
    {
        for (int j = 0; j < noCols; j++)
        {
             std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/*  PrintRowVector outputs vector of dimensions (1 x m), where (integer) 
    noCols = m. */
void PrintRowVector(int noCols, double* vector)
{
    std::cout << "[ ";
    for (int i = 0; i < noCols; i++)
    {
        std::cout << vector[i] << " ";
    }
    std::cout << "]" << std::endl;
}

/*  PrintColVector outputs vector of dimensions (n x 1), where (integer) 
    noRows = n. */
void PrintColVector(int noRows, double* vector)
{
    for (int i = 0; i < noRows; i++)
    {
        std::cout << vector[i] << std::endl;
    }
}

/* AddVectors performs A + B, where vectors A and B have dimensions (n x 1) or 
   (1 x n), where (integer) vectorSize = n. */
double* AddVectors(int vectorSize, double* vectorA, double* vectorB)
{
    double* result = AllocateVector(vectorSize);
    for (int i = 0; i < vectorSize; i++)
    {
        result[i] = vectorA[i] + vectorB[i];
    }
    return result;
}

/* AddVectors performs A - B, where vectors A and B have dimensions (n x 1) or 
   (1 x n), where (integer) vectorSize = n. */
double* SubtractVectors(int vectorSize, double* vectorA, double* vectorB)
{
    double* result = AllocateVector(vectorSize);
    for (int i = 0; i < vectorSize; i++)
    {
        result[i] = vectorA[i] - vectorB[i];
    }
    return result;
}

/*  DotProduct performs A.B, where vectors A and B have dimensions (n x 1), 
    where (integer) noRows = n. */
double DotProduct(int noRows, double* vectorA, double* vectorB)
{
    double result = 0.0;
    for (int i = 0; i < noRows; i++)
    {
        result += vectorA[i]*vectorB[i];
    }
    return result;   
}

/*  ScaleVector multiplies every element in (double*) vector by 
    (double) scaleFactor. */
double* ScaleVector(int noRows, double scaleFactor, double* vector)
{
    double* result = AllocateVector(noRows);
    for (int i = 0; i < noRows; i++)
    {
        result[i] = scaleFactor*vector[i];
    }
    return result;
}

/*  ComputeTwoNorm computes the 2-norm of (double*) vector, 
    ||(double*) vector||_2, where (double*) vector is of dimensions (n x 1) and 
    (integer) noRows = n. */
double ComputeTwoNorm(int noRows, double* vector)
{
    double result = 0.0;
    for (int i = 0; i < noRows; i++)
    {
        result += pow(vector[i], 2);
    }
    result = sqrt(result);
    return result;
}