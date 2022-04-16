#ifdef __cplusplus
extern "C" {
#endif
#include "ngspice/gmres.h"
#include "ngspice/smpdefs.h"
#include "ngspice/ngspice.h"
#include "ngspice/devdefs.h"
#include "ngspice/sperror.h"
#include "spdefs.h"
#include "spconfig.h"
#ifdef XSPICE
#include "ngspice/enh.h"
/* gtri - add - wbk - 11/26/90 - add include for MIF global data */
#include "ngspice/mif.h"
/* gtri - end - wbk - 11/26/90 */
#endif
#ifdef __cplusplus
}
#endif

#include <math.h>
#include <assert.h>
#include <time.h>

//Everything is 1-indexed!
//calculate b = Ax
void Mult(MatrixPtr A, RealVector x, RealVector b) {
    spMultiply(A, b, x, NULL, NULL);
}

double Dot(RealVector x, RealVector y, int n) {
    double result = 0.0;
    for (int i = 1; i <= n; i++)
        result += x[i] * y[i];
    return result;
}

double Norm(RealVector x, int n) {
    double result = 0.0;
    for (int i = 1; i <= n; i++)
        result += x[i] * x[i];
    return sqrt(result);
}

//calculate y = cx
void VectorConstMult(RealVector x, double c, RealVector y, int n) {
    for (int i = 1; i <= n; i++)
        y[i] = c * x[i];
}

//calculate x = x + b*y
void VectorAdd(RealVector x, RealVector y, double b, int n) {
    for (int i = 1; i <= n; i++)
        x[i] += b * y[i];
}

void PrintVector(RealVector x, int n) {
    for (int i = 1; i <= n; i++)
        printf("%f ", x[i]);
    printf("\n");
    fflush(stdout);
}

void constructGMRES(GMRESarr *arr) {
    arr->n = 0;
    arr->hadPrec = 0;
    arr->PrecNeedReset = 1;
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

void initPreconditoner(MatrixPtr Matrix, GMRESarr *arr) {
    clock_t start = clock();

    int error = SMPpreOrder(arr->Prec);
    error = spOrderAndFactor(arr->Prec, NULL, Matrix->RelThreshold, Matrix->AbsThreshold, YES);

    if (!arr->hadPrec) {
        arr->Precarr = SP_MALLOC(double, arr->Prec->Elements);
        arr->idx = SP_MALLOC(int, arr->Prec->Elements);
    }
    continuify(arr);
    clock_t end = clock();
    arr->Prectime = (double) (end - start) / CLOCKS_PER_SEC;
    arr->hadPrec = 1;
    arr->PrecNeedReset = 0;
    printf("Preconditioner time: %f\n", arr->Prectime);
    fflush(stdout);
}

void initGMRES(GMRESarr *arr, int n) {
    printf("initGMRES\n");
    printf("n = %d\n", n);
    fflush(stdout);
    if (arr->n != 0)
        return;
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
    arr->n = 0;
}

void fastSolve(GMRESarr *arr, RealVector RHS, RealVector Solution) {
    RealVector Intermediate;
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

int gmresSolvePreconditoned(GMRESarr *arr, MatrixPtr Matrix, double *RHS, double *Solution)
{
    int n = Matrix->Size, iters = 0;
    double eps = 1e-6;
    double *x0 = arr->x0;

    clock_t start = clock();
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
            r0[i] = 0.0, w[i] = 0.0, q[i] = 0.0;
        for (int i = 1; i <= n; i++)
            v[1][i] = r0[i];
        
        Mult(Matrix, x0, r0); //r0 = Ax0
        for (int i = 1; i <= n; i++)
            r0[i] = RHS[i] - r0[i]; //r0 = b - Ax0
        //SMPsolve(Prec, r0, r0); //r0 = (Prec)^{-1} r0
        fastSolve(arr, r0, r0);
        double beta = Norm(r0, n);

        VectorConstMult(r0, 1.0 / beta, v[1], n); //v[1] = r0/beta
        q[1] = beta;

        double relres = 0;
        int m = 0;
        for (int j = 1; j <= maxiter; ++j) {
            v[j + 1] = SP_MALLOC(double, n + 1);
            Mult(Matrix, v[j], w); //w[j] = Av[j]
            //SMPsolve(Prec, w, w); //w[j] = (Prec)^{-1} w[j]
            fastSolve(arr, w, w);
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
            m = j;
        }
        iters = MAX(m, iters);
        double *y = arr->y;

        //solve h*y = q
        for (int i = m; i >= 1; --i) {
            y[i] = q[i];
            for (int j = i + 1; j <= m; ++j)
                y[i] -= h[i][j] * y[j];
            y[i] /= h[i][i];
        }

//        double *estb = SP_MALLOC(double, n + 1);
        for (int i = 1; i <= m; i++)
            VectorAdd(x0, v[i], y[i], n);
//        spSolveU(Prec, x0, x0);
//        Mult(Matrix, sol, estb);
//        double soldif = 0.0;
//        for (int i = 1; i <= n; i++)
//            soldif += (estb[i] - RHS[i]) * (estb[i] - RHS[i]);
//        soldif = sqrt(soldif);
//        printf("reboot = %d, soldif = %e\n", reboot, soldif);
//        printf("relres = %e\n", relres);
//        fflush(stdout);
        if (relres < eps)
            break;
    }
    for (int I = n; I > 0; I--)
        Solution[I] = x0[I];
    clock_t end = clock();
    arr->GMREStime += (double)(end - start) / CLOCKS_PER_SEC;
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
    return 1;
    if (strcmp(name, "Resistor") == 0 ||
        strcmp(name, "Vsource") == 0 ||
        strcmp(name, "Isource") == 0) {
        return 1;
    } else {
        return 0;
    }
}

//extract the linear components of the circuit while loading matrix
int CKTloadPreconditioner(CKTcircuit *ckt, GMRESarr *arr) {
    int i;
    int size;
    double startTime;
    CKTnode *node;
    int error;
#ifdef XSPICE
    /* gtri - begin - Put resistors to ground at all nodes */
    /*   SMPmatrix  *matrix; maschmann : deleted , because unused */

    double gshunt;
    int num_nodes;

    /* gtri - begin - Put resistors to ground at all nodes */
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
            printf("device type %s\n",
                       DEVices[i]->DEVpublic.name);
#ifdef STEPDEBUG
            if (noncon != ckt->CKTnoncon) {
                printf("device type %s nonconvergence\n",
                       DEVices[i]->DEVpublic.name);
                noncon = ckt->CKTnoncon;
            }
#endif /* STEPDEBUG */
            if (error) return(error);
        }
    }
    for (i = 0; i < DEVmaxnum; i++) {
        if (DEVices[i] && DEVices[i]->DEVload && ckt->CKThead[i] && !isLinear(DEVices[i]->DEVpublic.name)) {
            error = DEVices[i]->DEVload (ckt->CKThead[i], ckt);
            if (ckt->CKTnoncon)
                ckt->CKTtroubleNode = 0;
            printf("device type %s\n",
                       DEVices[i]->DEVpublic.name);
#ifdef STEPDEBUG
            if (noncon != ckt->CKTnoncon) {
                printf("device type %s nonconvergence\n",
                       DEVices[i]->DEVpublic.name);
                noncon = ckt->CKTnoncon;
            }
#endif /* STEPDEBUG */
            if (error) return(error);
        }
    }
    MatrixPtr Matrix = ckt->CKTmatrix;
    int n = Matrix->Size;
    if (!arr->hadPrec) {
        arr->Prec = spCreate(n, 0, &error);
    } else {
        SMPclear(arr->Prec);
    }
    for (int I = 1; I <= n; I++) {
        ElementPtr pElement = Matrix->FirstInCol[I];
        while (pElement != NULL)
        {
            int Row = Matrix->IntToExtRowMap[pElement->Row];
            int Col = Matrix->IntToExtColMap[I];
            SMPaddElt(arr->Prec, Row, Col, pElement->Real);
            pElement = pElement->NextInCol;
        }
    }
#ifdef XSPICE
    /* gtri - add - wbk - 11/26/90 - reset the MIF init flags */

    /* init is set by CKTinit and should be true only for first load call */
    g_mif_info.circuit.init = MIF_FALSE;

    /* anal_init is set by CKTdoJob and is true for first call */
    /* of a particular analysis type */
    g_mif_info.circuit.anal_init = MIF_FALSE;

    /* gtri - end - wbk - 11/26/90 */

    /* gtri - begin - Put resistors to ground at all nodes. */
    /* Value of resistor is set by new "rshunt" option.     */

    if (ckt->enh->rshunt_data.enabled) {
        gshunt = ckt->enh->rshunt_data.gshunt;
        num_nodes = ckt->enh->rshunt_data.num_nodes;
        for (i = 0; i < num_nodes; i++) {
            *(ckt->enh->rshunt_data.diag[i]) += gshunt;
        }
    }
    /* gtri - end - Put resistors to ground at all nodes */
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
                    /* DAG: Original CIDER fix. If above fix doesn't work,
                     * revert to this.
                     */
                    /*
                     *  ckt->CKTrhs[node->number] += 1.0e10 * node->nodeset;
                     *  *(node->ptr) += 1.0e10;
                     */
                }
            }
        }
        if ((ckt->CKTmode & MODETRANOP) && (!(ckt->CKTmode & MODEUIC))) {
            for (node = ckt->CKTnodes; node; node = node->next) {
                if (node->icGiven) {
                    printf("node %s\n", node->name);
                    if (ZeroNoncurRow(ckt->CKTmatrix, ckt->CKTnodes,
                                      node->number)) {
                        /* Original code:
                         ckt->CKTrhs[node->number] += 1.0e10 * node->ic;
                        */
                        ckt->CKTrhs[node->number] = 1.0e10 * node->ic *
                                                      ckt->CKTsrcFact;
                        *(node->ptr) += 1.0e10;
                    } else {
                        /* Original code:
                          ckt->CKTrhs[node->number] = node->ic;
                        */
                        ckt->CKTrhs[node->number] = node->ic*ckt->CKTsrcFact; /* AlansFixes */
                        *(node->ptr) = 1;
                    }
                    /* DAG: Original CIDER fix. If above fix doesn't work,
                     * revert to this.
                     */
                    /*
                     *  ckt->CKTrhs[node->number] += 1.0e10 * node->ic;
                     *  *(node->ptr) += 1.0e10;
                     */
                }
            }
        }
    }
    /* SMPprint(ckt->CKTmatrix, stdout); if you want to debug, this is a
    good place to start ... */
    ckt->CKTstat->STATloadTime += SPfrontEnd->IFseconds()-startTime;
    return(OK);
}