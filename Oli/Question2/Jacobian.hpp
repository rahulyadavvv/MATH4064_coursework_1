double** ComputeApprox(double* x, double* (*pF)(double* x), int n, double epsilon);

double** ComputeInverseApprox(double * x, double* (*pF)(double* x), int n, double epsilon);

void ComputeInverseApprox(double** matrix, double* x, double* (*pF)(double* x), double epsilon, int n);
