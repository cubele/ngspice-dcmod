#include "ngspice/spmatrix.h"
#include "ngspice/cktdefs.h"

#ifndef GMRES_H
#define GMRES_H

#define GMRESmaxiter (50)
typedef struct GMRESarr GMRESarr;

#ifdef __cplusplus
extern "C" {
#endif
int gmresSolvePreconditoned(GMRESarr *, MatrixPtr, double *, double *);
void initPreconditoner(MatrixPtr, GMRESarr *);
void constructGMRES(GMRESarr **);
void initGMRES(GMRESarr *, int);
void freeGMRES(GMRESarr *);
void continuify(GMRESarr *);
int NIiter_fast(CKTcircuit *, GMRESarr *, int);
int CKTloadPreconditioner(CKTcircuit *, GMRESarr *);
#ifdef __cplusplus
}
#endif

#endif