#include "feGRASS.hpp"

double checkEdge(matGraph *g, int u, int v, double w) {
    return g->checkEdge(u, v, w);
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