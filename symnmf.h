#define BETA 0.5
#define ERR "An Error Has Occurred"
#define EPSILON 0.0001
#define MAX_ITER 300

double **parseFile(char *fname, int *n, int *d);
void findArrayDimentions(FILE *fp, int *n, int *d);

double **calcSymilarityMatrix(double **dataPoints, int n, int d);
double *calcDiagonalDegreeMatrix(double **A, int n);
double **calcNormalizedSymilarityMatrix(double *D, double **A, int n);
double **calcOptimalDecompMatrix(double **initialH, double **W, int n, int k);

void updateDecompMatrix(double **initialH, double **W, double **next, double **HHT, double **HHTH, double **WH, int n, int k);

double calcEuclideanDistanceSquared(double *x, double *y, int d);
double calcFrobeniusNorm(double **A, double **B, int n, int k);
int calcMatrixMult(double **A, double **B, double **C, int n, int k, int m, int transpose);

void printMatrix(double **M, int n);
void printDiagMatrix(double *D, int n);