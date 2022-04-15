#include "ngspice/spmatrix.h"

#ifndef GMRES_H
#define GMRES_H

#define GMRESmaxiter (50)
typedef struct gmreselements {
    double h[GMRESmaxiter + 3][GMRESmaxiter + 3];
    double c[GMRESmaxiter + 3], s[GMRESmaxiter + 3], y[GMRESmaxiter + 3];
    double *x0;
    double *v[GMRESmaxiter + 3];
    double *r0, *w, *q;
    double *sol;
    int n;

    MatrixPtr Prec;
    double *Precarr; //contiguous array for preconditioner to acclerate solving
    int *idx;
    int *rowind, *colind; //next row/col starts at colind[n]
    int Lsize, Usize;
    double GMREStime, Prectime;
    int origiters, totaliters;
    int hadPrec, PrecNeedReset;
}GMRESarr;

int gmresSolvePreconditoned(GMRESarr *arr, MatrixPtr Matrix, double *RHS, double *Solution);
void initPreconditoner(MatrixPtr Matrix, GMRESarr *arr);
void constructGMRES(GMRESarr *arr);
void initGMRES(GMRESarr *arr, int n);
void freeGMRES(GMRESarr *arr);
void continuify(GMRESarr *arr);

#endif