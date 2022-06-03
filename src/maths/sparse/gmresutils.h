#ifndef GMRESUTILS_H
#define GMRESUTILS_H

#include "ngspice/smpdefs.h"
#include "ngspice/spmatrix.h"
#include "ngspice/cktdefs.h"
#include "spdefs.h"
#include "spconfig.h"
#include "feGRASS.hpp"

#define GMRESmaxiter (300)
#define ratiodiff (0.03)
#define initratio (0.12)
#define GMRESeps (1e-12)
#define REALeps (1e-6)
#define GMRESreboots (6)
struct GMRESarr{
    double h[GMRESmaxiter + 3][GMRESmaxiter + 3];
    double c[GMRESmaxiter + 3], s[GMRESmaxiter + 3], y[GMRESmaxiter + 3];
    double *x0;
    double *v[GMRESmaxiter + 3];
    double *r0, *w, *q;
    double *sol;
    double eps;
    int n;

    MatrixPtr Prec, Orig;
    double *Precarr; //contiguous array for preconditioner to acclerate solving
    int *idx;
    int LUsize;
    int *rowind, *colind; //next row/col starts at colind[n]
    int Lsize, Usize;
    double GMREStime, Prectime;
    int origiters, totaliters, extraiters, totalrounds;
    int hadPrec, PrecNeedReset, iterno, PrecNeedUpdate;

    matGraph *G;
    double ratio;
};

void Mult(MatrixPtr A, double *x, double *b);
double Dot(double *x, double *y, int n);
double Norm(double *x, int n);
void VectorConstMult(double *x, double c, double *y, int n);
void VectorAdd(double *x, double *y, double b, int n);
void PrintVector(double *x, int n);
void copyMatrix(MatrixPtr Matrix, MatrixPtr dest);
void LoadGmin(MatrixPtr Matrix, double Gmin);

#endif