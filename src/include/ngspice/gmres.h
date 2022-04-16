#include "ngspice/spmatrix.h"
#include "ngspice/cktdefs.h"

#ifndef GMRES_H
#define GMRES_H

#define GMRESmaxiter (50)
typedef struct {
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

#ifdef __cplusplus
extern "C" {
#endif
int gmresSolvePreconditoned(GMRESarr *, MatrixPtr, double *, double *);
void initPreconditoner(MatrixPtr, GMRESarr *);
void constructGMRES(GMRESarr *);
void initGMRES(GMRESarr *, int);
void freeGMRES(GMRESarr *);
void continuify(GMRESarr *);
int NIiter_fast(CKTcircuit *, GMRESarr *, int);
int CKTloadPreconditioner(CKTcircuit *, GMRESarr *);
#ifdef __cplusplus
}
#endif

#endif