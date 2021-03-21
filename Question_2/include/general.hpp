double* AllocateVector(int noRows);
void DeallocateVector(double* vector);
double** AllocateMatrix(int noRows, int noCols);
void DeallocateMatrix(int noRows, double** matrix);
void PrintMatrix(int noRows, int noCols, double** matrix);
void PrintRowVector(int noCols, double* vector);
void PrintColVector(int noRows, double* vector);
double* AddVectors(int vectorSize, double* vectorA, double* vectorB);
double* SubtractVectors(int vectorSize, double* vectorA, double* vectorB);
double DotProduct(int noRows, double* vectorA, double* vectorB);
double* ScaleVector(int noRows, double scaleFactor, double* vector);
double ComputeTwoNorm(int noRows, double* vector);
double* MultiplyMatrixVector(int noRows, int noCols, double** pMatrix, double* 
pVector);
double ComputeOneNorm(int noRows, double* x);
double** Inverse2x2(double** pMatrix);

