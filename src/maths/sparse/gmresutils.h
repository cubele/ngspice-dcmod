#ifndef GMRESUTILS_H
#define GMRESUTILS_H

#include "ngspice/smpdefs.h"
#include "ngspice/spmatrix.h"
#include "ngspice/cktdefs.h"
#include "spdefs.h"
#include "spconfig.h"

void Mult(MatrixPtr A, double *x, double *b);
double Dot(double *x, double *y, int n);
double Norm(double *x, int n);
void VectorConstMult(double *x, double c, double *y, int n);
void VectorAdd(double *x, double *y, double b, int n);
void PrintVector(double *x, int n);
void copyMatrix(MatrixPtr Matrix, MatrixPtr dest);
void LoadGmin(MatrixPtr Matrix, double Gmin);

#endif