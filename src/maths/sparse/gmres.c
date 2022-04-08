#include "gmres.h"
#include "spdefs.h"
#include "spconfig.h"
#include "ngspice/smpdefs.h"
#include <math.h>
#include <assert.h>

//Everything is 1-indexed!
//calculate b = Ax
void Mult(MatrixPtr A, RealVector x, RealVector b)
{
    spMultiply(A, b, x, NULL, NULL);
}

double Dot(RealVector x, RealVector y, int n)
{
    double result = 0.0;
    for (int i = 1; i <= n; i++)
        result += x[i] * y[i];
    return result;
}

double Norm(RealVector x, int n)
{
    double result = 0.0;
    for (int i = 1; i <= n; i++)
        result += x[i] * x[i];
    return sqrt(result);
}

//calculate y = cx
void VectorConstMult(RealVector x, double c, RealVector y, int n)
{
    for (int i = 1; i <= n; i++)
        y[i] = c * x[i];
}

//calculate x = x + b*y
void VectorAdd(RealVector x, RealVector y, double b, int n)
{
    for (int i = 1; i <= n; i++)
        x[i] += b * y[i];
}

void PrintVector(RealVector x, int n)
{
    for (int i = 1; i <= n; i++)
        printf("%f ", x[i]);
    printf("\n");
    fflush(stdout);
}

MatrixPtr getPreconditoner(MatrixPtr Matrix, double Gmin)
{
    int error = 0, n = Matrix->Size;
    MatrixPtr Prec = NULL;
    Prec = spCreate(n, 0, &error);

    //inefficient way of copying the matrix
    for (int I = 1; I <= n; I++)
    {
        int added = 0;
        ElementPtr pElement = Matrix->FirstInCol[I];
        while (pElement != NULL)
        {
            int Row = Matrix->IntToExtRowMap[pElement->Row];
            int Col = Matrix->IntToExtColMap[I];
            SMPaddElt(Prec, Row, Col, pElement->Real);
            added = 1;
            pElement = pElement->NextInCol;
        }
    }

    //SMPprint(Prec, "./Preconditoner.m");
    //SMPprint(Matrix, "./Matrix.m");

    error = SMPluFac(Prec, 0, 0);
    printf("luFac:%d\n", error);
    return Prec;
}

const int maxiter = 200;
int gmresSolvePreconditoned(MatrixPtr Matrix, MatrixPtr Prec, double *RHS, double *Solution, double Gmin)
{
    if (Matrix->Complex) {
        return 1;
    }

    int n = Matrix->Size;
    double eps = 1e-8;

    double *x0 = SP_MALLOC(double, n + 1);
    for (int i = 1; i <= n; i++)
        x0[i] = 0;

    for (int reboot = 0; reboot < 1; reboot++) {
        double h[maxiter + 3][maxiter + 3];
        
        double *r0 = SP_MALLOC(double, n + 1);
        for (int i = 1; i <= n; i++)
            r0[i] = 0.0;
        double *v[maxiter + 3];
        double c[maxiter + 3], s[maxiter + 3];
        double *w = SP_MALLOC(double, n + 1);
        for (int i = 1; i < maxiter + 3; i++)
            v[i] = SP_MALLOC(double, n + 1);
        double *q = SP_MALLOC(double, n + 1);
        for (int i = 1; i <= n; i++)
            v[1][i] = r0[i];
        for (int i = 1; i <= n; i++)
            q[i] = 0, w[i] = 0;
        
        Mult(Matrix, x0, r0); //r0 = Ax0
        for (int i = 1; i <= n; i++)
            r0[i] = RHS[i] - r0[i]; //r0 = b - Ax0
        SMPsolve(Prec, r0, r0); //r0 = (Prec)^{-1} r0
        double beta = Norm(r0, n);

        VectorConstMult(r0, 1.0 / beta, v[1], n); //v[1] = r0/beta
        q[1] = beta;

        double relres = 0;
        int m = 0;
        for (int i = 1; i <= m + 1; ++i) {
            for (int j = 1; j <= m + 1; ++j) {
                h[i][j] = 0;
            }
        }
        for (int j = 1; j <= maxiter; ++j) {
            Mult(Matrix, v[j], w); //w[j] = Av[j]
            SMPsolve(Prec, w, w); //w[j] = (Prec)^{-1} w[j]
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
        printf("finished GMRES iters = %d\n", m);
        fflush(stdout);
        double y[maxiter + 3];

        //solve h*y = q
        for (int i = m; i >= 1; --i) {
            y[i] = q[i];
            for (int j = i + 1; j <= m; ++j)
                y[i] -= h[i][j] * y[j];
            y[i] /= h[i][i];
        }

        double *sol = SP_MALLOC(double, n + 1);
        for (int i = 1; i <= n; i++)
            sol[i] = x0[i];
        for (int i = 1; i <= m; i++)
            VectorAdd(sol, v[i], y[i], n);
        double *estb = SP_MALLOC(double, n + 1);
        Mult(Matrix, sol, estb);
        double soldif = 0.0;
        for (int i = 1; i <= n; i++)
            soldif += (estb[i] - RHS[i]) * (estb[i] - RHS[i]);
        soldif = sqrt(soldif);
        printf("reboot = %d, soldif = %e\n", reboot, soldif);
        printf("relres = %e\n", relres);
        fflush(stdout);

        for (int i = 1; i <= n; i++)
            x0[i] = sol[i];
        if (relres < eps)
            break;
        SP_FREE(sol);
        SP_FREE(w);
        SP_FREE(q);
        for (int i = 1; i < maxiter + 3; ++i)
            SP_FREE(v[i]);
        SP_FREE(r0);
    }
    for (int I = n; I > 0; I--)
        Solution[I] = x0[I];
    SP_FREE(x0);
    return 0;
}