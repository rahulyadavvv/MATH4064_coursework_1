double* SolveODESystemEuler(double* initial, int max_iter, double* (*pFunction)(double* x), double h);

double* ComputeEulerIteration(int iter, double* initial, double* (*pFunction)(double* x), double h);

double Compute2Norm(double* x, double* y, int n);

