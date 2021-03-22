double** ComputeApproximateJacobian(double* x, double* (*pF)(double* x), int n, 
                                    double epsilon);

double** ComputeApproximateJacobianInverse(double* x, double* (*pF)(double* x), 
                                            double epsilon, int n);

void ComputeApproximateJacobianInverse(double** matrix, double* x, 
                            double* (*pF)(double* x), double epsilon, int n);
