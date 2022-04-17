#include "gmresutils.h"

//Everything is 1-indexed
//calculate b = Ax
void Mult(MatrixPtr A, double *x, double *b) {
    spMultiply(A, b, x, NULL, NULL);
}

double Dot(double *x, double *y, int n) {
    double result = 0.0;
    for (int i = 1; i <= n; i++)
        result += x[i] * y[i];
    return result;
}

double Norm(double *x, int n) {
    double result = 0.0;
    for (int i = 1; i <= n; i++)
        result += x[i] * x[i];
    return sqrt(result);
}

//calculate y = cx
void VectorConstMult(double *x, double c, double *y, int n) {
    for (int i = 1; i <= n; i++)
        y[i] = c * x[i];
}

//calculate x = x + b*y
void VectorAdd(double *x, double *y, double b, int n) {
    for (int i = 1; i <= n; i++)
        x[i] += b * y[i];
}

void PrintVector(double *x, int n) {
    for (int i = 1; i <= n; i++)
        printf("%f ", x[i]);
    printf("\n");
    fflush(stdout);
}