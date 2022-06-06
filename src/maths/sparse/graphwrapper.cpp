#include "feGRASS.hpp"
#include "trialmodel.hpp"

double checkEdge(matGraph *g, int u, int v, double w, int *isedge) {
    return g->checkEdge(u, v, w, isedge);
}

void initGraph(matGraph **g, int n) {
    *g = new matGraph(n);
}

void addEdge(matGraph *g, int u, int v, double w) {
    g->addEdge(u, v, w);
}

int sparsify(matGraph *g, double p) {
    return g->sparsify(p);
}

double findRatio(matGraph *g) {
    for (double p = 0.1; p < 1; p += 0.01) {
        printf("ratio=%f ", p);
        printf("%d\n", sparsify(g, p));
    }
    return 0.1;
}

void clearGraph(matGraph *g) {
    g->clear();
}

void initTrialModel(trialModel **t, int n) {
    *t = new trialModel();
}

void addTrial(trialModel *T, double ratio, double lutime, double gmrestime) {
    T->addTrial(ratio, lutime, gmrestime);
}

double getRatio(trialModel *T) {
    return T->getRatio();
}

void inittempMat(tempMat **t, int n) {
    *t = new tempMat(n);
}
void addElem(tempMat *T, int i, int j, double val) {
    T->add(i, j, val);
}

void clearMat(tempMat *T) {
    T->clearElem();
}

double checkElem(tempMat *T, int u, int v) {
    return T->check(u, v);
}