#include "ngspice/spmatrix.h"
#include "ngspice/cktdefs.h"

#ifndef GMRES_H
#define GMRES_H

typedef struct GMRESarr GMRESarr;

int gmresSolvePreconditoned(GMRESarr *, CKTcircuit *, MatrixPtr, double, double *, double *);
void initPreconditoner(MatrixPtr, double, double, GMRESarr *);
void constructGMRES(GMRESarr **);
void initGMRES(GMRESarr *, int);
void freeGMRES(GMRESarr *);
void continuify(GMRESarr *);
int NIiter_fast(CKTcircuit *, GMRESarr *, int);
int CKTloadPreconditioner(CKTcircuit *, GMRESarr *);
void resetPrec(GMRESarr *);

#endif