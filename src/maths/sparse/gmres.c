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
#include "trialmodel.hpp"
#include "gmresutils.h"
#include <math.h>
#include <assert.h>
#include <time.h>

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
            if (ABS(pElement->Real) != 0) {
                arr->idx[cnt] = pElement->Row;
                arr->Precarr[cnt++] = pElement->Real;
            }
            pElement = pElement->NextInCol;
        }
        arr->colind[I] = cnt;
    }

    /* Backward Substitution. Solves Ux = c.*/
    for (I = Size; I > 0; I--) {
        pElement = Matrix->Diag[I]->NextInRow;
        while (pElement != NULL)
        {
            if (ABS(pElement->Real) != 0) {
                arr->idx[cnt] = pElement->Col;
                arr->Precarr[cnt++] = pElement->Real;
            }
            pElement = pElement->NextInRow;
        }
        arr->rowind[I] = cnt;
    }
}

void initPreconditoner(MatrixPtr Matrix, double relthres, double absthres, GMRESarr *arr) {
    int error;
    clock_t start = clock();
    if (!arr->hadPrec) {
        error = SMPpreOrder(arr->Prec);
        clock_t endorder = clock();
        printf("Preconditioner order time: %f\n", (double)(endorder - start) / CLOCKS_PER_SEC);
        spSetReal(Matrix);
        error = spOrderAndFactor(arr->Prec, NULL, relthres, absthres, YES);
    } else {
        error = spFactor(arr->Prec);
    }
    clock_t endfactor = clock();
    printf("preconditioner LU factor time: %f\n", (double)(endfactor - start) / CLOCKS_PER_SEC);

    if (!arr->hadPrec && arr->LUsize == 0) {
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
    if (!arr->hadPrec) {
        arr->Prectime += (double) (end - start) / CLOCKS_PER_SEC;
        printf("Total Preconditioner construction time: %f\n", arr->Prectime);
    }
    arr->hadPrec = 1;
    arr->PrecNeedReset = 0;
}

//spSolve(arr->Prec, RHS, Solution, NULL, NULL);
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

int gmresSolvePreconditoned(GMRESarr *arr, CKTcircuit *ckt, MatrixPtr origMatrix, double Gmin, double *RHS, double *Solution) {
    MatrixPtr Matrix = origMatrix;
    int n = Matrix->Size, iters = 0;
    double eps = arr->eps;
    double *x0 = arr->x0;
    double *tmp = SP_MALLOC(double, n + 1);
    int maxiter = GMRESmaxiter;
    double **v = arr->v;
    double *h[GMRESmaxiter + 3];
    for (int i = 0; i < GMRESmaxiter + 3; i++)
        h[i] = arr->h[i];
    double *c = arr->c, *s = arr->s;
    double *r0 = arr->r0;
    double *w = arr->w;
    double *q = arr->q;
    for (int i = 1; i <= n; i++) {
        x0[i] = 0.0;
    }
    for (int reboot = 0; reboot < GMRESreboots; reboot++) {
        for (int i = 1; i <= n; i++) {
            r0[i] = 0.0, w[i] = 0.0, q[i] = 0.0;
        }
        
        Mult(Matrix, x0, r0); //r0 = Ax0
        for (int i = 1; i <= n; i++) {
            r0[i] = RHS[i] - r0[i]; //r0 = b - Ax0
        }
        double beta = Norm(r0, n);

        VectorConstMult(r0, 1.0 / beta, v[1], n); //v[1] = r0/beta
        q[1] = beta;

        double relres = 0;
        int m = 0;
        for (int j = 1; j <= maxiter; ++j) {
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
            m = j;
            if (relres < eps) {
                break;
            }
            if (j % 100 == 0) {
                printf("relres: %e iter: %d\n", relres, j);
            }
        }
        printf("relres after iters: %e\n", relres);
        iters += m;
        double *y = arr->y;

        //solve h*y = q
        for (int i = m; i >= 1; --i) {
            y[i] = q[i];
            for (int j = i + 1; j <= m; ++j) {
                y[i] -= h[i][j] * y[j];
            }
            y[i] /= h[i][i];
            if (ABS(h[i][i]) < 1e-16) {
                printf("h[%d][%d] = %e\n", i, i, h[i][i]);
            }
        }

        for (int i = 1; i <= n; ++i) {
            tmp[i] = 0;
        }
        //x_m = x_0 + Prec^{-1}V_my_m
        for (int i = 1; i <= m; i++) {
            VectorAdd(tmp, v[i], y[i], n);
        }
        fastSolve(arr, tmp, tmp);
        VectorAdd(x0, tmp, 1, n);

        if (relres < eps) {
            break;
        }
/*
        double absdiff = 0, reldiff = 0;
        double totabsdiff = 0, totreldiff = 0;
        for (int i = 1; i <= n; i++) {
            absdiff = MAX(absdiff, ABS(RHS[i] - tmp[i]));
            totabsdiff += ABS(RHS[i] - tmp[i]);
            if (MAX(ABS(RHS[i]), ABS(tmp[i]) > 1e-8)) {
                double rdiff = ABS(RHS[i] - tmp[i]) / MAX(ABS(RHS[i]), ABS(tmp[i]));
                totreldiff += rdiff;
                reldiff = MAX(reldiff, rdiff);
            }
        }
        printf("max absdiff: %e max reldiff: %e\n", absdiff, reldiff);
        printf("avg absdiff: %e avg reldiff: %e\n", totabsdiff / n, totreldiff / n);
*/
    }
    /*
    double diff = 0;
    Mult(Matrix, x0, tmp);
    for (int i = 1; i <= n; i++)
        diff += (RHS[i] - tmp[i]) * (RHS[i] - tmp[i]);
    diff = sqrt(diff);
    printf("GMRES end, real residual: %e\n", diff);
    diff = diff / Norm(RHS, n);
    printf("real relative residual: %e\n", diff);
*/

    for (int I = n; I > 0; I--)
        Solution[I] = x0[I];

    SP_FREE(tmp);
    printf("GMRES iters: %d\n", iters);
    return iters;
}

void initGMRES(GMRESarr *arr, int n) {
    if (arr->n != 0)
        return;
    printf("initGMRES\n");
    printf("n = %d\n", n);
    arr->n = n;
    arr->x0 = SP_MALLOC(double, n + 1);
    for (int i = 0; i <= n; i++)
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
    arr->eps = GMRESeps;
    arr->precChanged = 0;
    arr->trialno = 0;
    arr->NIitercnt = 0;
    arr->LUsize = 0;
    arr->stable = 0;
    arr->ratioset = 0;
    arr->precUpdate = 0;
    initGraph(&arr->G, n);
    initTrialModel(&arr->T, n);
    inittempMat(&arr->M, n);
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