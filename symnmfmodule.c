#define PY_SSIZE_T_CLEAN
#include "Python.h"

#include "symnmf.h"
#include <stdlib.h>
#include <stdio.h>

static PyObject *symnmf(PyObject *self, PyObject *args)
{
    PyObject *initialHList, *WList, *HList;
    PyObject *row;
    PyObject *item;
    int n, k, i, j;
    int returnValue;

    double **initialH, **W, **H;
    returnValue = 0;

    if (!PyArg_ParseTuple(args, "OO", &initialHList, &WList))
    {
        return NULL;
    }

    n = (int)PyObject_Length(initialHList);
    k = (int)PyObject_Length(PyList_GetItem(initialHList, 0));

    initialH = (double **)calloc(n, sizeof(double *));
    W = (double **)calloc(n, sizeof(double *));

    if (initialH == NULL || W == NULL)
    {
        printf(ERR);
        returnValue = 1;
        goto FREE_X;
    }

    for (i = 0; i < n; i++)
    {
        row = PyList_GetItem(initialHList, i);
        initialH[i] = (double *)calloc(k, sizeof(double));
        if (initialH[i] == NULL)
        {
            printf(ERR);
            returnValue = 1;
            goto FREE_X;
        }
        for (j = 0; j < k; j++)
        {
            item = PyList_GetItem(row, j);
            initialH[i][j] = PyFloat_AsDouble(item);
        }
    }

    for (i = 0; i < n; i++)
    {
        row = PyList_GetItem(WList, i);
        W[i] = (double *)calloc(n, sizeof(double));
        if (W[i] == NULL)
        {
            printf(ERR);
            returnValue = 1;
            goto FREE_X;
        }
        for (j = 0; j < n; j++)
        {
            item = PyList_GetItem(row, j);
            W[i][j] = PyFloat_AsDouble(item);
        }
    }

    H = calcOptimalDecompMatrix(initialH, W, n, k);

    HList = PyList_New(n);

    for (int i = 0; i < n; i++)
    {
        row = PyList_New(k);
        for (int j = 0; j < k; j++)
        {
            item = Py_BuildValue("d", H[i][j]);
            PyList_SetItem(row, j, item);
        }
        PyList_SetItem(HList, i, row);
    }

FREE_X:
    for (int i = 0; i < n; i++)
    {
        free(initialH[i]);
        free(W[i]);
    }
    free(initialH);
    free(W);

    if (returnValue == 1)
    {
        return NULL;
    }
    return HList;
}

static PyObject *sym(PyObject *self, PyObject *args)
{
    PyObject *XList, *AList;
    PyObject *row;
    PyObject *item;
    int n, d, i, j;
    int returnValue;

    double **X, **A;
    returnValue = 0;

    if (!PyArg_ParseTuple(args, "O", &XList))
    {
        return NULL;
    }

    n = (int)PyObject_Length(XList);
    d = (int)PyObject_Length(PyList_GetItem(XList, 0));

    X = (double **)calloc(n, sizeof(double *));

    if (X == NULL)
    {
        printf(ERR);
        returnValue = 1;
        goto FREE_X;
    }

    for (i = 0; i < n; i++)
    {
        row = PyList_GetItem(XList, i);
        X[i] = (double *)calloc(d, sizeof(double));
        if (X[i] == NULL)
        {
            printf(ERR);
            returnValue = 1;
            goto FREE_X;
        }
        for (j = 0; j < d; j++)
        {
            item = PyList_GetItem(row, j);
            X[i][j] = PyFloat_AsDouble(item);
        }
    }

    A = calcSymilarityMatrix(X, n, d);

    AList = PyList_New(n);

    for (int i = 0; i < n; i++)
    {
        row = PyList_New(n);
        for (int j = 0; j < n; j++)
        {
            item = Py_BuildValue("d", A[i][j]);
            PyList_SetItem(row, j, item);
        }
        PyList_SetItem(AList, i, row);
    }

    for (int i = 0; i < n; i++)
    {
        free(A[i]);
    }
    free(A);

FREE_X:
    for (int i = 0; i < n; i++)
    {
        free(X[i]);
    }
    free(X);

    if (returnValue == 1)
    {
        return NULL;
    }
    return AList;
}

static PyObject *ddg(PyObject *self, PyObject *args)
{
    PyObject *XList, *DList;
    PyObject *row;
    PyObject *item;
    int n, d, i, j;
    int returnValue;

    double **X, **A, *D;
    returnValue = 0;

    if (!PyArg_ParseTuple(args, "O", &XList))
    {
        return NULL;
    }

    n = (int)PyObject_Length(XList);
    d = (int)PyObject_Length(PyList_GetItem(XList, 0));

    X = (double **)calloc(n, sizeof(double *));

    if (X == NULL)
    {
        printf(ERR);
        returnValue = 1;
        goto FREE_X;
    }

    for (i = 0; i < n; i++)
    {
        row = PyList_GetItem(XList, i);
        X[i] = (double *)calloc(d, sizeof(double));
        if (X[i] == NULL)
        {
            printf(ERR);
            returnValue = 1;
            goto FREE_X;
        }
        for (j = 0; j < d; j++)
        {
            item = PyList_GetItem(row, j);
            X[i][j] = PyFloat_AsDouble(item);
        }
    }

    A = calcSymilarityMatrix(X, n, d);
    D = calcDiagonalDegreeMatrix(A, n);

    DList = PyList_New(n);

    for (int i = 0; i < n; i++)
    {
        item = Py_BuildValue("d", D[i]);
        PyList_SetItem(DList, i, item);
    }

    for (int i = 0; i < n; i++)
    {
        free(A[i]);
    }
    free(A);
    free(D);
FREE_X:
    for (int i = 0; i < n; i++)
    {
        free(X[i]);
    }
    free(X);

    if (returnValue == 1)
    {
        return NULL;
    }
    return DList;
}

static PyObject *norm(PyObject *self, PyObject *args)
{
    PyObject *XList, *WList;
    PyObject *row;
    PyObject *item;
    int n, d, i, j;
    int returnValue;

    double **X, **A, **W, *D;
    returnValue = 0;

    if (!PyArg_ParseTuple(args, "O", &XList))
    {
        return NULL;
    }

    n = (int)PyObject_Length(XList);
    d = (int)PyObject_Length(PyList_GetItem(XList, 0));

    X = (double **)calloc(n, sizeof(double *));

    if (X == NULL)
    {
        printf(ERR);
        returnValue = 1;
        goto FREE_X;
    }

    for (i = 0; i < n; i++)
    {
        row = PyList_GetItem(XList, i);
        X[i] = (double *)calloc(d, sizeof(double));
        if (X[i] == NULL)
        {
            printf(ERR);
            returnValue = 1;
            goto FREE_X;
        }
        for (j = 0; j < d; j++)
        {
            item = PyList_GetItem(row, j);
            X[i][j] = PyFloat_AsDouble(item);
        }
    }

    A = calcSymilarityMatrix(X, n, d);
    D = calcDiagonalDegreeMatrix(A, n);
    W = calcNormalizedSymilarityMatrix(D, A, n);

    WList = PyList_New(n);

    for (int i = 0; i < n; i++)
    {
        row = PyList_New(n);
        for (int j = 0; j < n; j++)
        {
            item = Py_BuildValue("d", W[i][j]);
            PyList_SetItem(row, j, item);
        }
        PyList_SetItem(WList, i, row);
    }

    for (int i = 0; i < n; i++)
    {
        free(W[i]);
        free(A[i]);
    }
    free(W);
    free(A);
    free(D);

FREE_X:
    for (int i = 0; i < n; i++)
    {
        free(X[i]);
    }
    free(X);

    if (returnValue == 1)
    {
        return NULL;
    }
    return WList;
}

static PyMethodDef symnmfMethods[] = {
    {"symnmf",
     (PyCFunction)symnmf,
     METH_VARARGS,
     PyDoc_STR("find k means for set of datapoints. pass- centroids: array of initialized centroids, dataPoints: array of datapoints iter: maximum iterations, e: accaptable accuracy for the means.")},
    {"sym",
     (PyCFunction)sym,
     METH_VARARGS,
     PyDoc_STR("find k means for set of datapoints. pass- centroids: array of initialized centroids, dataPoints: array of datapoints iter: maximum iterations, e: accaptable accuracy for the means.")},
    {"ddg",
     (PyCFunction)ddg,
     METH_VARARGS,
     PyDoc_STR("find k means for set of datapoints. pass- centroids: array of initialized centroids, dataPoints: array of datapoints iter: maximum iterations, e: accaptable accuracy for the means.")},
    {"norm",
     (PyCFunction)norm,
     METH_VARARGS,
     PyDoc_STR("find k means for set of datapoints. pass- centroids: array of initialized centroids, dataPoints: array of datapoints iter: maximum iterations, e: accaptable accuracy for the means.")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef symnmfModule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    NULL,
    -1,
    symnmfMethods};

PyMODINIT_FUNC PyInit_symnmfmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfModule);
    if (!m)
    {
        return NULL;
    }
    return m;
}