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
        return 1;
    }

    if (strcmp(goal, "ddg") == 0)
    {
        A = calcSymilarityMatrix(dataPoints, n, d);
        D = calcDiagonalDegreeMatrix(A, n);
        printDiagMatrix(D, n);
        return 1;
    }

    if (strcmp(goal, "norm") == 0)
    {
        A = calcSymilarityMatrix(dataPoints, n, d);
        D = calcDiagonalDegreeMatrix(A, n);
        W = calcNormalizedSymilarityMatrix(D, A, n);
        printMatrix(W, n);
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
    for (i = 0; i < *n; i++)
    {
        dataPoints[i] = (double *)calloc(*d, sizeof(double));
        if (dataPoints[i] == NULL)
        {
            printf("An Error Has Occurred in malloc \n");
            /* TODO: set return value to failure and free datapoints*/
        }
        for (j = 0; (j < *d) && (fscanf(fp, "%lf", &current) != EOF); j++)
        {
            fgetc(fp);
            dataPoints[i][j] = current;
        }
    }

    return dataPoints;
}

void findArrayDimentions(FILE *fp, int *n, int *d)
{
    /* this section finds d and n by counting \n and ','. this relies HEAVILY on the format of the input */
    char c;
    while (!feof(fp))
    {
        c = fgetc(fp);
        if (c == ',')
        {
            *d = *d + 1;
        }
        if (c == '\n')
        {
            *n = *n + 1;
        }
    }
    printf("%i,%i\n", *n, *d);
    *n = *n + 1;
    *d = ((*d) / (*n)) + 1;
    return;
}

double **calcSymilarityMatrix(double **dataPoints, int n, int d)
{
    int i, j;
    double **A, distance;

    A = (double **)calloc(n, sizeof(double *));
    for (i = 0; i < n; i++)
    {
        A[i] = (double *)calloc(n, sizeof(double));
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
    for (i = 0; i < n; i++)
    {
        Diag[i] = 1 / sqrt(D[i]);
    }

    /* We use the fact that (DAD) = d_i * a_ij * d_j*/
    W = (double **)calloc(n, sizeof(double *));
    for (i = 0; i < n; i++)
    {
        W[i] = (double *)calloc(n, sizeof(double));
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
    double result;
    for (i = 0; i < d; i++)
    {
        result += pow(x[i] - y[i], 2);
    }
    return result;
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