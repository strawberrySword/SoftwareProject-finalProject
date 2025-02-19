#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

int main(int argc, char *argv[])
{
    char *fname;
    double **dataPoints, **A, *D, **W;
    char *goal;
    int n, d;
    n = 0;
    d = 0;
    if (argc != 3)
    {
        perror("Invalid number of arguments.");
        return 0;
    }

    fname = argv[2];
    goal = argv[1];
    dataPoints = parseFile(fname, &n, &d);
    if (strcmp(goal, "sym") == 0)
    {
        A = calcSymilarityMatrix(dataPoints, n, d);
        printMatrix(A, n);
        freeMatrix(A, n);
        freeMatrix(dataPoints, n);
        return 1;
    }

    if (strcmp(goal, "ddg") == 0)
    {
        A = calcSymilarityMatrix(dataPoints, n, d);
        D = calcDiagonalDegreeMatrix(A, n);
        printDiagMatrix(D, n);
        freeMatrix(A, n);
        free(D);
        freeMatrix(dataPoints, n);
        return 1;
    }

    if (strcmp(goal, "norm") == 0)
    {
        A = calcSymilarityMatrix(dataPoints, n, d);
        D = calcDiagonalDegreeMatrix(A, n);
        W = calcNormalizedSymilarityMatrix(D, A, n);
        printMatrix(W, n);
        freeMatrix(A, n);
        free(D);
        freeMatrix(W, n);
        freeMatrix(dataPoints, n);
        return 1;
    }

    return 1;
}

double **parseFile(char *fname, int *n, int *d)
{
    int i, j;

    FILE *fp;
    double **dataPoints;
    double current;

    fp = fopen(fname, "r");
    if (fp == NULL)
    {
        perror("An Error Has Occured\n");
        exit(1);
    }

    findArrayDimentions(fp, n, d);
    rewind(fp);

    dataPoints = (double **)calloc(*n, sizeof(double *));
    if (dataPoints == NULL)
    {
        perror("An Error Has Occured\n");
        exit(1);
    }
    for (i = 0; i < *n; i++)
    {
        dataPoints[i] = (double *)calloc(*d, sizeof(double));
        if (dataPoints[i] == NULL)
        {
            printf("An Error Has Occurred in malloc \n");
            freeMatrix(dataPoints, i);
            exit(1);
        }
        for (j = 0; (j < *d) && (fscanf(fp, "%lf", &current) != EOF); j++)
        {
            fgetc(fp);
            dataPoints[i][j] = current;
        }
        if (ferror(fp))
        {
            perror("Error reading from file");
            fclose(fp);
            exit(1);
        }
    }

    return dataPoints;
}

void freeMatrix(double **matrix, int length)
{
    int i;
    for (i = 0; i < length; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
    return;
}

void findArrayDimentions(FILE *fp, int *n, int *d)
{
    /* this section finds d and n by counting \n and ','. this relies HEAVILY on the format of the input */
    char c;
    double temp;
    while (fscanf(fp, "%lf", &temp) != EOF)
    {
        c = fgetc(fp);
        if (*n == 0 && c == ',')
        {
            *d = *d + 1;
        }
        if (c == '\n' || c == EOF)
        {
            *n = *n + 1;
        }
    }
    if (ferror(fp))
    {
        perror("An Error Has Occured");
        fclose(fp);
        exit(1);
    }
    /* *n = *n + 1;*/
    *d = *d + 1;
    return;
}

double **calcSymilarityMatrix(double **dataPoints, int n, int d)
{
    int i, j;
    double **A, distance;

    A = (double **)calloc(n, sizeof(double *));
    if (A == NULL)
    {
        perror("An Error Has Occurred in malloc\n");
        exit(1);
    }
    for (i = 0; i < n; i++)
    {
        A[i] = (double *)calloc(n, sizeof(double));
        if (A[i] == NULL)
        {
            printf("An Error Has Occurred in malloc \n");
            freeMatrix(dataPoints, i);
            exit(1);
        }
        for (j = 0; j < n; j++)
        {
            distance = calcEuclideanDistanceSquared(dataPoints[i], dataPoints[j], d);
            A[i][j] = exp(-distance / 2);
        }
        A[i][i] = 0;
    }

    return A;
}

double *calcDiagonalDegreeMatrix(double **A, int n)
{
    int i, j;
    double sum, *D;

    D = (double *)calloc(n, sizeof(double));
    if (D == NULL)
    {
        perror("An Error Has Occurred in malloc\n");
        exit(1);
    }
    for (i = 0; i < n; i++)
    {
        sum = 0;
        for (j = 0; j < n; j++)
        {
            sum += A[i][j];
        }
        D[i] = sum;
    }

    return D;
}

double **calcNormalizedSymilarityMatrix(double *D, double **A, int n)
{
    int i, j;
    double *Diag, **W;

    Diag = (double *)calloc(n, sizeof(double));
    if (Diag == NULL)
    {
        perror("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++)
    {
        Diag[i] = 1 / sqrt(D[i]);
    }

    /* We use the fact that (DAD) = d_i * a_ij * d_j*/
    W = (double **)calloc(n, sizeof(double *));
    if (W == NULL)
    {
        perror("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++)
    {
        W[i] = (double *)calloc(n, sizeof(double));
        if (W[i] == NULL)
        {
            printf("An Error Has Occurred\n");
            freeMatrix(W, i);
            exit(1);
        }
        for (j = 0; j < n; j++)
        {
            W[i][j] = A[i][j] * Diag[i] * Diag[j];
        }
    }

    return W;
}

double calcEuclideanDistanceSquared(double *x, double *y, int d)
{
    int i;
    double result = 0;
    for (i = 0; i < d; i++)
    {
        result += pow(x[i] - y[i], 2);
    }
    return result;
}

double **calcOptimalDecompMatrix(double **initialH, double **W, int n, int k)
{
    int i;
    double frobeniusNorm, **temp, **next, **HHT, **HHTH, **WH;
    /* Initialize all matrices*/
    next = (double **)calloc(n, sizeof(double *));

    HHT = (double **)calloc(n, sizeof(double *));

    HHTH = (double **)calloc(n, sizeof(double *));

    WH = (double **)calloc(n, sizeof(double *));

    if (next == NULL || HHT == NULL || HHTH == NULL || WH == NULL)
    {
        printf("An Error Has Occurred\n");
        free(next);
        free(HHT);
        free(HHTH);
        free(WH);
        exit(1);
    }

    for (i = 0; i < n; i++)
    {
        next[i] = (double *)calloc(k, sizeof(double));
        HHTH[i] = (double *)calloc(k, sizeof(double));
        WH[i] = (double *)calloc(k, sizeof(double));
        HHT[i] = (double *)calloc(n, sizeof(double));

        if (WH[i] == NULL || next[i] == NULL || HHTH[i] == NULL || HHT[i] == NULL)
        {
            perror("An Error Has Occurred\n");
            freeMatrix(next, i);
            freeMatrix(HHTH, i);
            freeMatrix(WH, i);
            freeMatrix(HHT, i);
            exit(1);
        }
    }

    /* main loop*/
    for (i = 0; i < MAX_ITER; i++)
    {
        updateDecompMatrix(initialH, W, next, HHT, HHTH, WH, n, k);

        frobeniusNorm = calcFrobeniusNorm(next, initialH, n, k);

        temp = initialH;
        initialH = next;
        next = temp;

        if (frobeniusNorm < EPSILON)
            return initialH;
    }
    return initialH;
}

void updateDecompMatrix(double **initialH, double **W, double **next, double **HHT, double **HHTH, double **WH, int n, int k)
{
    int i, j;

    calcMatrixMultTranspose(initialH, initialH, HHT, n, k, n);

    calcMatrixMult(HHT, initialH, HHTH, n, n, k);

    calcMatrixMult(W, initialH, WH, n, n, k);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < k; j++)
        {
            next[i][j] = initialH[i][j] * (1 - BETA + BETA * (WH[i][j] / HHTH[i][j]));
        }
    }
}

double calcFrobeniusNorm(double **A, double **B, int n, int k)
{
    double ret;
    int i, j;
    ret = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < k; j++)
        {
            ret += (A[i][j] - B[i][j]) * (A[i][j] - B[i][j]);
        }
    }
    return ret;
}

int calcMatrixMult(double **A, double **B, double **C, int n, int k, int m)
{
    int i, j, l;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            C[i][j] = 0;
            for (l = 0; l < k; l++)
            {
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
    return 0;
}

int calcMatrixMultTranspose(double **A, double **B, double **C, int n, int k, int m)
{
    int i, j, l;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            C[i][j] = 0;
            for (l = 0; l < k; l++)
            {
                C[i][j] += A[i][l] * B[j][l];
            }
        }
    }
    return 0;
}

void printMatrix(double **M, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        printf("%.4f", M[i][0]);
        for (j = 1; j < n; j++)
        {
            printf(",%.4f", M[i][j]);
        }
        printf("\n");
    }
}

void printDiagMatrix(double *D, int n)
{
    int i, j;

    for (i = 0; i < n; i++)
    {
        if (i == 0)
        {
            printf("%.4f", D[0]);
        }
        else
        {
            printf("0.0000");
        }
        for (j = 1; j < n; j++)
        {
            if (i == j)
            {
                printf(",%.4f", D[i]);
            }
            else
            {
                printf(",0.0000");
            }
        }
        printf("\n");
    }
}
