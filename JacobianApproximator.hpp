double** ApproximateJacobian( double* vector, double (*pF1)(double x, double y), 
                    double (*pF2)(double x, double y), double epsilon, int n);

double** ApproximateJacobianInverse(double* vector, double (*pF1)(double x, double y), 
                    double (*pF2)(double x, double y), double epsilon, int n);

double** ComputeApproximateJacobian(double* x, double* (*pF)(double* x), int n, 
                                    double epsilon);

