double* Vector(int length)
{
    double* vector;
    vector = new double[length];

    return vector;
}

double* SolveTridiagonalSystem(int n, double* lower, double* diagonal, double* upper, double* rhs)
{
    for(int i=1; i<n; i++)
    {
        diagonal[i] = diagonal[i] - upper[i-1]*lower[i]/diagonal[i-1];
        rhs[i] = rhs[i] - rhs[i-1]*lower[i]/diagonal[i-1];
    }

    double* x = Vector(n);
    for(int i=0; i<n; i++)
    {
        x[i] = 0;
    }

    x[n-1] = rhs[n-1]/diagonal[n-1];

    for(int i=n-2; i>-1; i--)
    {
        x[i] = (rhs[i] - upper[i]*x[i+1])/diagonal[i];
    }

    return x;
}

