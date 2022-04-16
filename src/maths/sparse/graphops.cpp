#include <algorithm>
#include <iostream>
#include "graphops.hpp"
//ngspice header files written in C
#ifdef __cplusplus
extern "C" {
#endif
    #include "ngspice/smpdefs.h"
#ifdef __cplusplus
}

void initGraph(graph **g, int n) {
    *g = new graph;
    (*g)->n = n;
    (*g)->m = 0;
    (*g)->adj.resize(n + 1);
    (*g)->w.resize(n + 1);
    (*g)->diag.resize(n + 1);
    (*g)->alle.reserve(n + 1);
    (*g)->nowe = 0;
}

void deleteGraph(graph *g) {
    delete g;
}

void addEdge(graph *g, int u, int v, int ou, int ov, double w) {
    g->alle.push_back(edge(u, v, ou, ov, w));
}

void sparsify(graph *g) {
    std::sort(g->alle.begin(), g->alle.end(), [](edge e1, edge e2) {
        return e1.w < e2.w;
    });
    for (auto e : g->alle) {
        g->insEdge(e);
    }
    g->m = g->alle.size();
    g->alle.clear();
}

void clearGraph(graph *g) {
    for (int i = 1; i <= g->n; ++i) {
        g->adj[i].clear();
        g->w[i].clear();
        g->diag[i] = 0;
    }
    g->alle.clear();
    g->m = 0;
    g->nowe = 0;
}

void setDiag(graph *g, int Row, double d) {
    g->diag[Row] = d;
}

void graphToMatrix(graph *g, MatrixPtr Prec) {
    int n = g->n;
    for (int i = 1; i <= n; i++) {
        double sum = 0;
        for (int j = 0; j < g->adj[i].size(); j++) {
            int v = g->adj[i][j];
            double w = g->w[i][j];
            sum += w;
            SMPaddElt(Prec, i, v, w);
        }
        SMPaddElt(Prec, i, i, -sum);
    }
}

double checkEdge(graph *g, int u, int v, double w) {
    edge now = g->alle[g->nowe];
    if (now.ou == u && now.ov == v) {
        g->nowe++;
        double diff = now.w - w;
        diff = diff < 0 ? -diff : diff;
        if (diff < 1e-16) {
            return 0;
        } else {
            return w;
        }
    } else {
        return w;
    }
}
#endif