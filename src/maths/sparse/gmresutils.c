#include "gmresutils.h"
#include "trialmodel.hpp"

//Everything is 1-indexed
//calculate b = Ax
void Mult(MatrixPtr A, double *x, double *b) {
    spMultiply(A, b, x, NULL, NULL);
}

double Dot(double *x, double *y, int n) {
    double result = 0.0;
    for (int i = 1; i <= n; i++)
        result += x[i] * y[i];
    return result;
}

double Norm(double *x, int n) {
    double result = 0.0;
    for (int i = 1; i <= n; i++)
        result += x[i] * x[i];
    return sqrt(result);
}

//calculate y = cx
void VectorConstMult(double *x, double c, double *y, int n) {
    for (int i = 1; i <= n; i++)
        y[i] = c * x[i];
}

//calculate x = x + b*y
void VectorAdd(double *x, double *y, double b, int n) {
    for (int i = 1; i <= n; i++)
        x[i] += b * y[i];
}

void PrintVector(double *x, int n) {
    for (int i = 1; i <= n; i++)
        printf("%f ", x[i]);
    printf("\n");
    fflush(stdout);
}


void copyMatrix(MatrixPtr Matrix, MatrixPtr dest) {
    int error, n = Matrix->Size;
    for (int I = 1; I <= n; I++) {
        ElementPtr pElement = Matrix->FirstInCol[I];
        while (pElement != NULL)
        {
            int Row = Matrix->IntToExtRowMap[pElement->Row];
            int Col = Matrix->IntToExtColMap[I];
            if (ABS(pElement->Real) > 1e-16) {
                SMPaddElt(dest, Row, Col, pElement->Real);
            }
            pElement = pElement->NextInCol;
        }
    }
}

void loadMatrix(MatrixPtr Matrix, GMRESarr *arr) {
    clearMat(arr->M);
    int n = Matrix->Size;
    for (int I = 1; I <= n; I++) {
        ElementPtr pElement = Matrix->FirstInCol[I];
        while (pElement != NULL)
        {
            int Row = Matrix->IntToExtRowMap[pElement->Row];
            int Col = Matrix->IntToExtColMap[I];
            int isedge = 0;
            double nv = checkEdge(arr->G, Row, Col, pElement->Real, &isedge);
            if (!isedge || nv != 0 || Row == Col) {
                addElem(arr->M, Row, Col, nv);
            }
            pElement = pElement->NextInCol;
        }
    }

    MatrixPtr Prec = arr->Prec;
    for (int I = 1; I <= n; I++) {
        ElementPtr pElement = Prec->FirstInCol[I];
        while (pElement != NULL)
        {
            int Row = Prec->IntToExtRowMap[pElement->Row];
            int Col = Prec->IntToExtColMap[I];
            pElement->Real = checkElem(arr->M, Row, Col);
            pElement = pElement->NextInCol;
        }
    }
}

void LoadGmin(MatrixPtr Matrix, double Gmin)
{
    int I;
    ArrayOfElementPtrs Diag;
    ElementPtr diag;

    /* Begin `LoadGmin'. */

    if (Gmin != 0.0) {
	Diag = Matrix->Diag;
	for (I = Matrix->Size; I > 0; I--) {
	    if ((diag = Diag[I]) != NULL)
		diag->Real += Gmin;
	}
    }
    return;
}