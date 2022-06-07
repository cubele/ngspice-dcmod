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

int CKTloadUpdate(CKTcircuit *ckt, GMRESarr *arr) {
    MatrixPtr Matrix = ckt->CKTmatrix;
    int n = Matrix->Size;
    SMPclear(arr->Prec);
    loadMatrix(Matrix, arr);
    return(OK);
}


#define FREE(x) do { if(x) { txfree(x); (x) = NULL; } } while(0)
//extract the linear components of the circuit while loading matrix
//is called when the preconditioner structure has changed
int CKTloadPreconditioner(CKTcircuit *ckt, GMRESarr *arr) {
    int i;
    int size;
    double startTime;
    CKTnode *node;
    int error;
    MatrixPtr Matrix = ckt->CKTmatrix;
    int n = Matrix->Size;
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
    if (!arr->hadPrec) {
        clearGraph(arr->G);
        for (i = 0; i < DEVmaxnum; i++) {
            if (DEVices[i] && DEVices[i]->DEVload && ckt->CKThead[i] && isLinear(DEVices[i]->DEVpublic.name)) {
                error = DEVices[i]->DEVload (ckt->CKThead[i], ckt);
                if (ckt->CKTnoncon)
                    ckt->CKTtroubleNode = 0;
                if (error) return(error);
            }
        }

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
    } else {
        for (i = 0; i < DEVmaxnum; i++) {
            if (DEVices[i] && DEVices[i]->DEVload && ckt->CKThead[i]) {
                error = DEVices[i]->DEVload (ckt->CKThead[i], ckt);
                if (ckt->CKTnoncon)
                    ckt->CKTtroubleNode = 0;
                if (error) return(error);
            }
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
        SMPclear(arr->Prec);
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
            int isedge = 0;
            double nv = checkEdge(arr->G, Row, Col, pElement->Real, &isedge);
            if (!isedge || nv != 0 || Row == Col) {
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
    clock_t end = clock();
    arr->Prectime = (double)(end - start) / CLOCKS_PER_SEC;
    printf("preconditioner creation time: %f\n", arr->Prectime);
    fflush(stdout);

    ckt->CKTstat->STATloadTime += SPfrontEnd->IFseconds()-startTime;
    return(OK);
}

int NIiter_fast(CKTcircuit *ckt, GMRESarr *arr, int maxIter)
{
    arr->NIitercnt++;
    double startTime, *OldCKTstate0 = NULL;
    int error, i, j;

    int iterno = 0;
    int ipass = 0;
    int firsttran = 0;
    if (arr->NIitercnt == 2) {
        firsttran = 1;
    }

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
            iterno++;
            printf("NI iterno = %d\n", iterno);
            int firstGMRES = 0;
            int calcGMREStime = 0, updateGMREStime = 0;
            arr->keptPrec = 0;
            if (arr->PrecNeedReset) {
                printf("Preconditioner structure change\n");
                error = CKTloadPreconditioner(ckt, arr);
                initPreconditoner(ckt->CKTmatrix, ckt->CKTpivotRelTol, ckt->CKTpivotAbsTol ,arr);
                arr->precChanged = 1;
                arr->stable = 0;
                arr->LUtime = 0;
                arr->GMREStime = 0;
                arr->totalrounds = 0;
                firstGMRES = 1;
            } else {
                if (!arr->ratioset) {
                    error = CKTload(ckt);
                    error = CKTloadUpdate(ckt, arr);
                    clock_t start = clock();
                    initPreconditoner(ckt->CKTmatrix, ckt->CKTpivotRelTol, ckt->CKTpivotAbsTol, arr);
                    clock_t end = clock();
                    printf("preconditioner update time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
                    if (arr->stable > 1 && arr->NIitercnt > 1 && !firsttran) {
                        calcGMREStime = 1;
                        arr->LUtime += (double)(end - start) / CLOCKS_PER_SEC;
                    }
                    if (!arr->precChanged && arr->NIitercnt > 1 && !firsttran) {
                        ++arr->stable;
                    }
                    arr->precChanged = 0;
                } else {
                    error = CKTload(ckt);
                    //if (arr->precUpdate) {
                    error = CKTloadUpdate(ckt, arr);
                    clock_t start = clock();
                    initPreconditoner(ckt->CKTmatrix, ckt->CKTpivotRelTol, ckt->CKTpivotAbsTol, arr);
                    clock_t end = clock();
                    printf("preconditioner update time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
                    arr->precUpdate = 0;
                    arr->LUtime = (double)(end - start) / CLOCKS_PER_SEC;
                    /*
                    updateGMREStime = 1;
                    } else {
                        arr->keptPrec = 1;
                    }
                    */
                }
            }
            SMPpreOrder(ckt->CKTmatrix);//just in case
            if (error) {
                ckt->CKTstat->STATnumIter += iterno;
                FREE(OldCKTstate0);
                return (error);
            }

            if (!OldCKTstate0)
                OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1);
            memcpy(OldCKTstate0, ckt->CKTstate0,
                   (size_t) ckt->CKTnumStates * sizeof(double));
            
            arr->iterno = iterno;
            clock_t start = clock();
            int iters = gmresSolvePreconditoned(arr, ckt, ckt->CKTmatrix, ckt->CKTdiagGmin, ckt->CKTrhs, ckt->CKTrhs);
            clock_t end = clock();
            printf("gmresSolvetime: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
            if (!arr->keptPrec && iters > GMRESmaxiter * GMRESreboots * 0.95) {
                arr->ratio += ratiodiff / 2;
                arr->PrecNeedReset = 1;
                arr->totalrounds = 0;
            }
            if (calcGMREStime) {
                arr->totalrounds++;
                arr->GMREStime += (double)(end - start) / CLOCKS_PER_SEC;
            }
            if (arr->ratioset) {
                double GMREStime = (double)(end - start) / CLOCKS_PER_SEC;
                printf("---GMREStime: %f LUtime: %f totaltime: %f---\n", GMREStime, arr->LUtime, GMREStime + arr->LUtime);
                /*
                if (GMREStime > arr->finalest) {
                    arr->precUpdate = 1;
                    printf("preconditioner update\n");
                }
                if (updateGMREStime) {
                    arr->GMREStime = (double)(end - start) / CLOCKS_PER_SEC;
                    arr->finalest = arr->GMREStime + arr->LUtime;
                }
                */
            } else if (arr->totalrounds >= 4) {
                addTrial(arr->T, arr->ratio, arr->LUtime / arr->totalrounds, arr->GMREStime / arr->totalrounds);
                ++arr->trialno;
                printf("ratio = %f, LUtime = %f, GMREStime = %f\n", arr->ratio, arr->LUtime / arr->totalrounds, arr->GMREStime / arr->totalrounds);
                if (arr->trialno == 4) {
                    printf("------------starting ratio calculation------------\n");
                    arr->ratio = getRatio(arr->T);
                    SMPdestroy(arr->Prec);
                    arr->PrecNeedReset = 1;
                    arr->hadPrec = 0;
                } else {
                    if (arr->trialno < 4) {
                        //test run
                        arr->ratio += ratiodiff;
                        arr->PrecNeedReset = 1;
                    } else {
                        arr->ratioset = 1;
                        arr->finalest = arr->LUtime / arr->totalrounds + arr->GMREStime / arr->totalrounds;
                        printf("------------final ratio = %f------------\n", arr->ratio);
                        arr->precUpdate = 0;
                    }
                }
            }

            ckt->CKTrhs[0] = 0;
            ckt->CKTrhsSpare[0] = 0;
            ckt->CKTrhsOld[0] = 0;

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
            
            if (iterno > 95) {
                ckt->CKTnoncon = 0;
            }
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
        firsttran = 0;
    }
    /*NOTREACHED*/
}