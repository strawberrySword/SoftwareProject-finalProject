double **parseFile(char *fname, int *n, int *d);
void findArrayDimentions(FILE *fp, int *n, int *d);

double **calcSymilarityMatrix(double **dataPoints, int n, int d);
double *calcDiagonalDegreeMatrix(double **A, int n);
double **calcNormalizedSymilarityMatrix(double *D, double **A, int n);

double calcEuclideanDistanceSquared(double *x, double *y, int d);

void printMatrix(double **M, int n);
void printDiagMatrix(double *D, int n);