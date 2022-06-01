#include "ngspice/gmres.h"
#include "ngspice/smpdefs.h"
#include "ngspice/ngspice.h"
#include "ngspice/devdefs.h"
#include "ngspice/sperror.h"
#include "spdefs.h"
#include "spconfig.h"
#ifdef XSPICE
#include "ngspice/enh.h"
#include "ngspice/mif.h"
#endif

#include "feGRASS.hpp"
#include "gmresutils.h"
#include <math.h>
#include <assert.h>
#include <time.h>
#define GMRESmaxiter (1000)
#define ratiodiff (0.03)
#define initratio (0.15)

struct GMRESarr{
    double h[GMRESmaxiter + 3][GMRESmaxiter + 3];
    double c[GMRESmaxiter + 3], s[GMRESmaxiter + 3], y[GMRESmaxiter + 3];
    double *x0;
    double *v[GMRESmaxiter + 3];
    double *r0, *w, *q;
    double *sol;
    int n;

    MatrixPtr Prec, Orig;
    double *Precarr; //contiguous array for preconditioner to acclerate solving
    int *idx;
    int LUsize;
    int *rowind, *colind; //next row/col starts at colind[n]
    int Lsize, Usize;
    double GMREStime, Prectime;
    int origiters, totaliters, extraiters, totalrounds;
    int hadPrec, PrecNeedReset;

    matGraph *G;
    double ratio;
};

void resetPrec(GMRESarr *arr) {
    arr->PrecNeedReset = 1;
}

void constructGMRES(GMRESarr **arr) {
    *arr = SP_MALLOC(GMRESarr, 1);
    (*arr)->n = 0;
    (*arr)->hadPrec = 0;
    (*arr)->PrecNeedReset = 1;
}

void continuify(GMRESarr *arr) {
    ElementPtr pElement;
    RealNumber Temp;
    int  I, Size;
    ElementPtr pPivot;
    MatrixPtr Matrix = arr->Prec;
    int cnt = 0;

    Size = Matrix->Size;

    /* Forward elimination. Solves Lc = b.*/
    for (I = 1; I <= Size; I++) {
        pPivot = Matrix->Diag[I];
        arr->Precarr[cnt++] = pPivot->Real;
        pElement = pPivot->NextInCol;
        while (pElement != NULL) {
            arr->idx[cnt] = pElement->Row;
            arr->Precarr[cnt++] = pElement->Real;
            pElement = pElement->NextInCol;
        }
        arr->colind[I] = cnt;
    }

    /* Backward Substitution. Solves Ux = c.*/
    for (I = Size; I > 0; I--) {
        pElement = Matrix->Diag[I]->NextInRow;
        while (pElement != NULL)
        {
            arr->idx[cnt] = pElement->Col;
            arr->Precarr[cnt++] = pElement->Real;
            pElement = pElement->NextInRow;
        }
        arr->rowind[I] = cnt;
    }
    assert(arr->Prec->Elements == cnt);
}

void initPreconditoner(MatrixPtr Matrix, double relthres, double absthres, GMRESarr *arr) {
    int error;
    clock_t start = clock();
    error = SMPpreOrder(arr->Prec);
    spSetReal(Matrix);
    error = spOrderAndFactor(arr->Prec, NULL, relthres, absthres, YES);
    clock_t endfactor = clock();
    printf("preconditioner LU factor time: %f\n", (double)(endfactor - start) / CLOCKS_PER_SEC);
    if (error) {
        printf("spOrderAndFactor: %d\n", error);
        fflush(stdout);
        exit(1);
    }

    if (!arr->hadPrec) {
        printf("New Preconditioner elements: %d\n", arr->Prec->Elements);
        arr->LUsize = arr->Prec->Elements;
        arr->Precarr = SP_MALLOC(double, arr->Prec->Elements);
        arr->idx = SP_MALLOC(int, arr->Prec->Elements);
    } else if (arr->Prec->Elements > arr->LUsize) {
        printf("Preconditioner elements increased from %d to %d\n", arr->LUsize, arr->Prec->Elements);
        arr->LUsize = arr->Prec->Elements;
        arr->Precarr = SP_REALLOC(arr->Precarr, double, arr->Prec->Elements);
        arr->idx = SP_REALLOC(arr->idx, int, arr->Prec->Elements);
    }
    continuify(arr);
    clock_t end = clock();
    arr->Prectime += (double) (end - start) / CLOCKS_PER_SEC;
    arr->hadPrec = 1;
    arr->PrecNeedReset = 0;
    printf("Total Preconditioner construction time: %f\n", arr->Prectime);
}

void initGMRES(GMRESarr *arr, int n) {
    if (arr->n != 0)
        return;
    printf("initGMRES\n");
    printf("n = %d\n", n);
    arr->n = n;
    arr->x0 = SP_MALLOC(double, n + 1);
    for (int i = 1; i <= n; i++)
        arr->x0[i] = 0.0;
    arr->r0 = SP_MALLOC(double, n + 1);
    arr->w = SP_MALLOC(double, n + 1);
    arr->q = SP_MALLOC(double, n + 1);
    arr->colind = SP_MALLOC(int, n + 1);
    arr->rowind = SP_MALLOC(int, n + 1);
    for (int i = 0; i < GMRESmaxiter + 3; i++)
        arr->v[i] = SP_MALLOC(double, n + 1);
    for (int i = 1; i <= GMRESmaxiter + 1; ++i) {
        for (int j = 1; j <= GMRESmaxiter + 1; ++j) {
            arr->h[i][j] = 0;
        }
    }

    int error;
    arr->Orig = spCreate(n, 0, &error);
    arr->ratio = initratio;
    initGraph(&arr->G, n);
}

void freeGMRES(GMRESarr *arr) {
    if (arr->n == 0)
        return;
    SP_FREE(arr->x0);
    SP_FREE(arr->r0);
    SP_FREE(arr->w);
    SP_FREE(arr->q);
    SP_FREE(arr->colind);
    SP_FREE(arr->rowind);
    for (int i = 0; i < GMRESmaxiter + 3; i++)
        SP_FREE(arr->v[i]);
    SP_FREE(arr->Precarr);
    SP_FREE(arr->idx);
    SP_FREE(arr);
}

void fastSolve(GMRESarr *arr, double * RHS, double * Solution) {
    double *Intermediate;
    RealNumber Temp;
    int I, *pExtOrder, Size;
    MatrixPtr Matrix = arr->Prec;
    int now = 0;

    Intermediate = Matrix->Intermediate;
    Size = Matrix->Size;

    /* Initialize Intermediate vector. */
    pExtOrder = &Matrix->IntToExtRowMap[Size];
    for (I = Size; I > 0; I--)
        Intermediate[I] = RHS[*(pExtOrder--)];

    /* Forward elimination. Solves Lc = b.*/
    for (I = 1; I <= Size; I++) {
	/* This step of the elimination is skipped if Temp equals zero. */
        if ((Temp = Intermediate[I]) != 0.0) {
            Intermediate[I] = (Temp *= arr->Precarr[now++]);
            while (now < arr->colind[I]) {
		        Intermediate[arr->idx[now]] -= Temp * arr->Precarr[now];
                ++now;
            }
        }
        now = arr->colind[I];
    }

    /* Backward Substitution. Solves Ux = c.*/
    for (I = Size; I > 0; I--) {
	    Temp = Intermediate[I];
        while (now < arr->rowind[I]) {
	        Temp -= arr->Precarr[now] * Intermediate[arr->idx[now]];
            ++now;
        }
        Intermediate[I] = Temp;
        now = arr->rowind[I];
    }

    /* Unscramble Intermediate vector while placing data in to Solution vector. */
    pExtOrder = &Matrix->IntToExtColMap[Size];
    for (I = Size; I > 0; I--)
        Solution[*(pExtOrder--)] = Intermediate[I];

    return;
}

int gmresSolvePreconditoned(GMRESarr *arr, MatrixPtr origMatrix, double Gmin, double *RHS, double *Solution)
{
    clock_t start = clock();
    //original matrix may contain spaces reserved for LUfac
    //copyMatrix(origMatrix, arr->Orig);
    LoadGmin(origMatrix, Gmin);
    MatrixPtr Matrix = origMatrix;
    int n = Matrix->Size, iters = 0;
    double eps = 1e-10;
    double *x0 = arr->x0;
    double *tmp = SP_MALLOC(double, n + 1);
    int maxiter = GMRESmaxiter;
    for (int reboot = 0; reboot < 1; reboot++) {
        double **v = arr->v;
        double *h[GMRESmaxiter + 3];
        for (int i = 0; i < GMRESmaxiter + 3; i++)
            h[i] = arr->h[i];
        double *c = arr->c, *s = arr->s;
        double *r0 = arr->r0;
        double *w = arr->w;
        double *q = arr->q;
        
        for (int i = 0; i <= n; i++)
            r0[i] = 1.0, w[i] = 0.0, q[i] = 0.0;
        for (int i = 1; i <= n; i++)
            v[1][i] = r0[i];
        
        Mult(Matrix, x0, r0); //r0 = Ax0
        for (int i = 1; i <= n; i++)
            r0[i] = RHS[i] - r0[i]; //r0 = b - Ax0
        //fastSolve(arr, r0, r0); //r0 = (Prec)^{-1} r0
        double beta = Norm(r0, n);

        VectorConstMult(r0, 1.0 / beta, v[1], n); //v[1] = r0/beta
        q[1] = beta;

        double relres = 0;
        int m = 0;
        for (int j = 1; j <= maxiter; ++j) {
            //Mult(Matrix, v[j], w); //w[j] = Av[j]
            //fastSolve(arr, w, w); //w[j] = (Prec)^{-1} w[j]
            fastSolve(arr, v[j], w); //w[j] = (Prec)^{-1} v[j]
            Mult(Matrix, w, w); //w[j] = A(Prec)^{-1}v[j]
            for (int i = 1; i <= j; ++i) {
                h[i][j] = Dot(v[i], w, n);
                VectorAdd(w, v[i], -h[i][j], n);//w[j] = w[j] - h[i][j]v[i]
            }
            h[j + 1][j] = Norm(w, n);
            for (int i = 1; i <= j - 1; ++i) {
                double hij = h[i][j], hij1 = h[i + 1][j];
                h[i][j] = c[i] * hij + s[i] * hij1;
                h[i + 1][j] = -s[i] * hij + c[i] * hij1;
            }
            if (ABS(h[j + 1][j]) < 1e-16) {
                m = j;
                break;
            }
            VectorConstMult(w, 1.0 / h[j + 1][j], v[j + 1], n); //v[j+1] = w/h[j+1][j]
            if (ABS(h[j][j]) > ABS(h[j + 1][j])) {
                double et = h[j + 1][j] / h[j][j];
                c[j] = 1.0 / sqrt(1.0 + et * et);
                s[j] = c[j] * et;
            } else {
                double et = h[j][j] / h[j + 1][j];
                s[j] = 1.0 / sqrt(1.0 + et * et);
                c[j] = s[j] * et;
            }

            h[j][j] = c[j] * h[j][j] + s[j] * h[j + 1][j];
            h[j + 1][j] = 0;

            q[j + 1] = -s[j] * q[j];
            q[j] = c[j] * q[j];
            relres = ABS(q[j + 1]) / beta;
            if (relres < eps) {
                m = j;
                break;
            }
            if (j % 100 == 0) {
                printf("relres: %e iter: %d\n", relres, j);
            }
            m = j;
        }
        printf("final relres: %e\n", relres);
        iters += m;
        double *y = arr->y;

        //solve h*y = q
        for (int i = m; i >= 1; --i) {
            y[i] = q[i];
            for (int j = i + 1; j <= m; ++j)
                y[i] -= h[i][j] * y[j];
            y[i] /= h[i][i];
        }

        //x_m = x_0 + Prec^{-1}V_my_m
        for (int i = 1; i <= m; i++) {
            VectorConstMult(v[i], y[i], tmp, n);
            fastSolve(arr, tmp, tmp);
            VectorAdd(x0, tmp, 1, n);
            //VectorAdd(x0, v[i], y[i], n);
        }
        if (relres < eps)
            break;
    }
    for (int I = n; I > 0; I--)
        Solution[I] = x0[I];
    Mult(Matrix, Solution, tmp);

    double diff = 0;
    for (int i = 1; i <= n; i++)
        diff += (RHS[i] - tmp[i]) * (RHS[i] - tmp[i]);
    diff = sqrt(diff);
    printf("GMRES end: %e residual\n", diff);

    double absdiff = 0, reldiff = 0;
    for (int i = 1; i <= n; i++) {
        absdiff = MAX(absdiff, ABS(RHS[i] - tmp[i]));
        double rdiff = ABS(RHS[i] - tmp[i]) / MAX(ABS(RHS[i]), ABS(tmp[i]));
        reldiff = MAX(reldiff, rdiff);
    }
    printf("max absdiff: %e max reldiff: %e\n", absdiff, reldiff);

    //SMPclear(arr->Orig);
    clock_t end = clock();
    arr->GMREStime += (double)(end - start) / CLOCKS_PER_SEC;
    printf("GMRES time: %g iters: %d\n", (double)(end - start) / CLOCKS_PER_SEC, iters);
    SP_FREE(tmp);
    return iters;
}

static int ZeroNoncurRow(SMPmatrix *matrix, CKTnode *nodes, int rownum)
{
    CKTnode     *n;
    double      *x;
    int         currents;
    currents = 0;
    for (n = nodes; n; n = n->next) {
        x = (double *) SMPfindElt(matrix, rownum, n->number, 0);
        if (x) {
            if (n->type == SP_CURRENT)
                currents = 1;
            else
                *x = 0.0;
        }
    }
    return currents;
}

int isLinear(char *name) {
    if (strcmp(name, "Resistor") == 0) {
        return 1;
    } else {
        return 0;
    }
}

#define FREE(x) do { if(x) { txfree(x); (x) = NULL; } } while(0)
//extract the linear components of the circuit while loading matrix
int CKTloadPreconditioner(CKTcircuit *ckt, GMRESarr *arr) {
    int i;
    int size;
    double startTime;
    CKTnode *node;
    int error;
#ifdef XSPICE
    double gshunt;
    int num_nodes;
#endif
    startTime = SPfrontEnd->IFseconds();
    size = SMPmatSize(ckt->CKTmatrix);
    for (i = 0; i <= size; i++) {
        ckt->CKTrhs[i] = 0;
    }
    SMPclear(ckt->CKTmatrix);
    for (i = 0; i < DEVmaxnum; i++) {
        if (DEVices[i] && DEVices[i]->DEVload && ckt->CKThead[i] && isLinear(DEVices[i]->DEVpublic.name)) {
            error = DEVices[i]->DEVload (ckt->CKThead[i], ckt);
            if (ckt->CKTnoncon)
                ckt->CKTtroubleNode = 0;
            if (error) return(error);
        }
    }

    MatrixPtr Matrix = ckt->CKTmatrix;
    int n = Matrix->Size;
    //load the linear part into graph
    for (int I = 1; I <= n; I++) {
        ElementPtr pElement = Matrix->FirstInCol[I];
        while (pElement != NULL)
        {
            int Row = Matrix->IntToExtRowMap[pElement->Row];
            int Col = Matrix->IntToExtColMap[I];
            if (pElement->Real != 0 && Row > Col) {
                addEdge(arr->G, Row, Col, pElement->Real);
            }
            pElement = pElement->NextInCol;
        }
    }

    for (i = 0; i < DEVmaxnum; i++) {
        if (DEVices[i] && DEVices[i]->DEVload && ckt->CKThead[i] && !isLinear(DEVices[i]->DEVpublic.name)) {
            error = DEVices[i]->DEVload (ckt->CKThead[i], ckt);
            if (ckt->CKTnoncon)
                ckt->CKTtroubleNode = 0;
            if (error) return(error);
        }
    }

#ifdef XSPICE
    g_mif_info.circuit.init = MIF_FALSE;
    g_mif_info.circuit.anal_init = MIF_FALSE;
    if (ckt->enh->rshunt_data.enabled) {
        gshunt = ckt->enh->rshunt_data.gshunt;
        num_nodes = ckt->enh->rshunt_data.num_nodes;
        for (i = 0; i < num_nodes; i++) {
            *(ckt->enh->rshunt_data.diag[i]) += gshunt;
        }
    }
#endif
    if (ckt->CKTmode & MODEDC) {
        /* consider doing nodeset & ic assignments */
        if (ckt->CKTmode & (MODEINITJCT | MODEINITFIX)) {
            /* do nodesets */
            for (node = ckt->CKTnodes; node; node = node->next) {
                if (node->nsGiven) {
                    printf("node %s\n", node->name);
                    if (ZeroNoncurRow(ckt->CKTmatrix, ckt->CKTnodes,
                                      node->number)) {
                        ckt->CKTrhs[node->number] = 1.0e10 * node->nodeset *
                                                      ckt->CKTsrcFact;
                        *(node->ptr) = 1e10;
                    } else {
                        ckt->CKTrhs[node->number] = node->nodeset *
                                                      ckt->CKTsrcFact;
                        *(node->ptr) = 1;
                    }
                }
            }
        }
        if ((ckt->CKTmode & MODETRANOP) && (!(ckt->CKTmode & MODEUIC))) {
            for (node = ckt->CKTnodes; node; node = node->next) {
                if (node->icGiven) {
                    printf("node %s\n", node->name);
                    if (ZeroNoncurRow(ckt->CKTmatrix, ckt->CKTnodes,
                                      node->number)) {
                        ckt->CKTrhs[node->number] = 1.0e10 * node->ic *
                                                      ckt->CKTsrcFact;
                        *(node->ptr) += 1.0e10;
                    } else {
                        ckt->CKTrhs[node->number] = node->ic*ckt->CKTsrcFact; /* AlansFixes */
                        *(node->ptr) = 1;
                    }
                }
            }
        }
    }

    if (!arr->hadPrec) {
        arr->Prec = spCreate(n, 0, &error);
    } else {
        SMPdestroy(arr->Prec);
        arr->Prec = spCreate(n, 0, &error);
    }

    clock_t start = clock();
    printf("ratio = %f\n", arr->ratio);
    sparsify(arr->G, arr->ratio);
    int nnz = 0, orignnz = 0;
    for (int I = 1; I <= n; I++) {
        ElementPtr pElement = Matrix->FirstInCol[I];
        while (pElement != NULL) {
            int Row = Matrix->IntToExtRowMap[pElement->Row];
            int Col = Matrix->IntToExtColMap[I];
            orignnz++;
            double nv = checkEdge(arr->G, Row, Col, pElement->Real);
            if (nv != 0 || Row == Col) {
                ++nnz;
                if(Row == Col) {
                    nv += ckt->CKTdiagGmin;
                }
                SMPaddElt(arr->Prec, Row, Col, nv);
            }
            pElement = pElement->NextInCol;
        }
    }
    printf("nnz in precondtioner: %d\n", nnz);
    printf("nnz in original matrix: %d\n", orignnz);
    clearGraph(arr->G);
    clock_t end = clock();
    arr->Prectime = (double)(end - start) / CLOCKS_PER_SEC;
    printf("preconditioner creation time: %f\n", arr->Prectime);
    fflush(stdout);

    ckt->CKTstat->STATloadTime += SPfrontEnd->IFseconds()-startTime;
    return(OK);
}

int NIiter_fast(CKTcircuit *ckt, GMRESarr *arr, int maxIter)
{
    double startTime, *OldCKTstate0 = NULL;
    int error, i, j;

    int iterno = 0;
    int ipass = 0;

    /* some convergence issues that get resolved by increasing max iter */
    if (maxIter < 100)
        maxIter = 100;

    if ((ckt->CKTmode & MODETRANOP) && (ckt->CKTmode & MODEUIC)) {
        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);
        error = CKTload(ckt);
        if (error)
            return(error);
        return(OK);
    }

#ifdef WANT_SENSE2
    if (ckt->CKTsenInfo) {
        error = NIsenReinit(ckt);
        if (error)
            return(error);
    }
#endif

    if (ckt->CKTniState & NIUNINITIALIZED) {
        error = NIreinit(ckt);
        if (error) {
            return(error);
        }
    }

    for (;;) {
        ckt->CKTnoncon = 0;
#ifdef NEWPRED
        if (!(ckt->CKTmode & MODEINITPRED))
#endif
        {
            int firstGMRES = 0;
            if (arr->PrecNeedReset) {
                printf("prec needed reset\n");
                error = CKTloadPreconditioner(ckt, arr);
                startTime = SPfrontEnd->IFseconds();
                initPreconditoner(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTpivotRelTol ,arr);
                ckt->CKTstat->STATdecompTime += SPfrontEnd->IFseconds() - startTime;
                arr->origiters = 0;
                arr->totaliters = 0;
                arr->GMREStime = 0;
                firstGMRES = 1;
            } else {
                error = CKTload(ckt);
            }
            iterno++;
            if (error) {
                ckt->CKTstat->STATnumIter += iterno;
                FREE(OldCKTstate0);
                return (error);
            }

            if (!OldCKTstate0)
                OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1);
            memcpy(OldCKTstate0, ckt->CKTstate0,
                   (size_t) ckt->CKTnumStates * sizeof(double));

#ifdef orig
            startTime = SPfrontEnd->IFseconds();
                error = SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol,
                                 ckt->CKTdiagGmin);
                ckt->CKTstat->STATdecompTime +=
                    SPfrontEnd->IFseconds() - startTime;
            printf("decomp time = %g\n", SPfrontEnd->IFseconds() - startTime);
            startTime = SPfrontEnd->IFseconds();
            SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare);
            ckt->CKTstat->STATsolveTime +=
                SPfrontEnd->IFseconds() - startTime;
#endif
#ifndef orig
            printf("NI iterno = %d\n", iterno);
            startTime = SPfrontEnd->IFseconds();
            int iters = gmresSolvePreconditoned(arr, ckt->CKTmatrix, ckt->CKTdiagGmin, ckt->CKTrhs, ckt->CKTrhs);
            arr->totaliters += iters;
            if (firstGMRES) {
                arr->origiters = iters;
                arr->extraiters = 0;
                arr->totalrounds = 1;
                if (iters > 0.9 * GMRESmaxiter) {
                    arr->ratio += ratiodiff;
                    arr->PrecNeedReset = 1;
                }
            } else {
                ++arr->totalrounds;
                int estiters = iters * arr->totalrounds;
                int reducediters = estiters - arr->totaliters;
                double avgitertime = arr->GMREStime / (double)arr->totaliters;
                if (iters > 0.9 * GMRESmaxiter || 
                    arr->Prectime < reducediters * avgitertime) {
                    printf("preconditioner reset after %d iters and %d rounds\n", arr->totaliters, arr->totalrounds);
                    printf("total GMREStime = %g, Prectime = %g\n", arr->GMREStime, arr->Prectime);
                    printf("avgitertime = %g\n", avgitertime);
                    printf("reducediters = %d\n", reducediters);
                    printf("estiters = %d\n", estiters);
                    arr->PrecNeedReset = 1;
                }
            }
            ckt->CKTstat->STATsolveTime +=
                SPfrontEnd->IFseconds() - startTime;
#endif

            ckt->CKTrhs[0] = 0;
            ckt->CKTrhsSpare[0] = 0;
            ckt->CKTrhsOld[0] = 0;

            if (iterno % 5 == 0 && iterno > 4) {
                arr->ratio += ratiodiff;
                arr->PrecNeedReset = 1;
            }

            if (iterno > maxIter) {
                ckt->CKTstat->STATnumIter += iterno;
                /* we don't use this info during transient analysis */
                if (ckt->CKTcurrentAnalysis != DOING_TRAN) {
                    FREE(errMsg);
                    errMsg = copy("Too many iterations without convergence");
                }
                FREE(OldCKTstate0);
                return(E_ITERLIM);
            }
            if ((ckt->CKTnoncon == 0) && (iterno != 1))
                ckt->CKTnoncon = NIconvTest(ckt);
            else
                ckt->CKTnoncon = 1;
        }

        if ((ckt->CKTnodeDamping != 0) && (ckt->CKTnoncon != 0) &&
            ((ckt->CKTmode & MODETRANOP) || (ckt->CKTmode & MODEDCOP)) &&
            (iterno > 1))
        {
            CKTnode *node;
            double diff, maxdiff = 0;
            for (node = ckt->CKTnodes->next; node; node = node->next)
                if (node->type == SP_VOLTAGE) {
                    diff = fabs(ckt->CKTrhs[node->number] - ckt->CKTrhsOld[node->number]);
                    if (maxdiff < diff)
                        maxdiff = diff;
                }

            if (maxdiff > 10) {
                double damp_factor = 10 / maxdiff;
                if (damp_factor < 0.1)
                    damp_factor = 0.1;
                for (node = ckt->CKTnodes->next; node; node = node->next) {
                    diff = ckt->CKTrhs[node->number] - ckt->CKTrhsOld[node->number];
                    ckt->CKTrhs[node->number] =
                        ckt->CKTrhsOld[node->number] + (damp_factor * diff);
                }
                for (i = 0; i < ckt->CKTnumStates; i++) {
                    diff = ckt->CKTstate0[i] - OldCKTstate0[i];
                    ckt->CKTstate0[i] = OldCKTstate0[i] + (damp_factor * diff);
                }
            }
        }

        if (ckt->CKTmode & MODEINITFLOAT) {
            if ((ckt->CKTmode & MODEDC) && ckt->CKThadNodeset) {
                if (ipass)
                    ckt->CKTnoncon = ipass;
                ipass = 0;
            }
            if (ckt->CKTnoncon == 0) {
                ckt->CKTstat->STATnumIter += iterno;
                FREE(OldCKTstate0);
                return(OK);
            }
        } else if (ckt->CKTmode & MODEINITJCT) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFIX;
            ckt->CKTniState |= NISHOULDREORDER;
        } else if (ckt->CKTmode & MODEINITFIX) {
            if (ckt->CKTnoncon == 0)
                ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
            ipass = 1;
        } else if (ckt->CKTmode & MODEINITSMSIG) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else if (ckt->CKTmode & MODEINITTRAN) {
            if (iterno <= 1)
                ckt->CKTniState |= NISHOULDREORDER;
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else if (ckt->CKTmode & MODEINITPRED) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else {
            ckt->CKTstat->STATnumIter += iterno;
            FREE(OldCKTstate0);
            return(E_INTERN);
            /* impossible - no such INITF flag! */
        }

        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);
    }
    /*NOTREACHED*/
}