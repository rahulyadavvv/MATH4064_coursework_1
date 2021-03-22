double* AssembleLowerDiagonal(int noSubIntervals);

double* AssembleDiagonal(int noSubIntervals);

double* AssembleUpperDiagonal(int noSubIntervals);

double* AssembleRHS(int noSubIntervals, double l, double (*pFunction)(double x), double (*pDerivativeFunction)(double x));

double ApproximateDerivative(double x, double (*pFunction)(double x));

double* AssembleRHSwApprox(int noSubIntervals, double intervalLength, double (*pFunction)(double x));

