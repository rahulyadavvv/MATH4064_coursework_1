double* Vector(int length)
{
    double* vector;
    vector = new double[length];
    return vector;
}

void DeleteVector(double* vector)
{
    delete[] vector;
}

double** Matrix(int n)
{
    double** matrix;
    matrix = new double*[n];
    for(int i=0; i<n; i++)
    {
        matrix[i] = new double[n];
    }
    return matrix;
}

void DeleteMatrix(double** matrix, int n)
{
    for(int i=0; i<n; i++)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}