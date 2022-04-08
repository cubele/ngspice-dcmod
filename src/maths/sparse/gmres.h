#include "ngspice/spmatrix.h"

int gmresSolvePreconditoned(MatrixPtr Matrix, MatrixPtr Prec, double *RHS, double *Solution, double Gmin);
MatrixPtr getPreconditoner(MatrixPtr Matrix, double Gmin);