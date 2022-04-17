#include <algorithm>
#include <iostream>
#include <cmath>
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
    (*g)->orige.reserve(n + 1);
    (*g)->nowe = 0;
}

void deleteGraph(graph *g) {
    delete g;
}

void addEdge(graph *g, int u, int v, double w) {
    g->alle.push_back(edge(u, v, w));
}

void addOrigEdge(graph *g, int u, int v, double w) {
    g->orige.push_back(edge(u, v, w));
}

void graph::insEdge(edge e) {
    adj[e.u].push_back(e.v);
    w[e.u].push_back(e.w);
}

int graph::sparsify() {
    std::sort(alle.begin(), alle.end(), [](edge a, edge b) {
        return a.w > b.w;
    });
    for (int i = 0; i < alle.size(); i++) {
        diag[alle[i].u] += alle[i].w;
    }
    unionset *us = new unionset(n);
    for (int i = 0; i < alle.size(); i++) {
        if (us->find(alle[i].u) != us->find(alle[i].v)) {
            insEdge(alle[i]);
            us->merge(alle[i].u, alle[i].v);
            ++m;
        } else {
            insEdge(alle[i]);
            ++m;
        }
    }
    printf("edges berore sparsify: %d\n", alle.size());
    printf("edges after sparsify: %d n: %d\n", m, n);
    fflush(stdout);
    delete us;
    return m;
}

int sparsify(graph *g) {
    return g->sparsify();
}

void clearGraph(graph *g) {
    for (int i = 1; i <= g->n; ++i) {
        g->adj[i].clear();
        g->w[i].clear();
        g->diag[i] = 0;
    }
    g->alle.clear();
    g->orige.clear();
    g->m = 0;
    g->nowe = 0;
}

void setDiag(graph *g, int Row, double d) {
    g->diag[Row] = d;
}

int graphToMatrix(graph *g, MatrixPtr Prec) {
    int n = g->n, nnz = 0;
    for (int i = 1; i <= n; i++) {
        g->diag[i] = 0;
    }
    for (int i = 1; i <= n; i++) {
        for (int j = 0; j < g->adj[i].size(); j++) {
            int v = g->adj[i][j];
            double w = g->w[i][j];
            g->diag[i] += w;
            g->diag[v] += w;
            nnz += 2;
            SMPaddElt(Prec, i, v, w);
            SMPaddElt(Prec, v, i, w);
        }
    }
    for (int i = 1; i <= n; ++i) {
        if (fabs(g->diag[i]) > 1e-16) {
            SMPaddElt(Prec, i, i, -g->diag[i]);
            ++nnz;
        }
    }
    return nnz;
}

double checkEdge(graph *g, int u, int v, double w, int *nnz) {
    edge now = g->orige[g->nowe];
    if (now.u == u && now.v == v) {
        g->nowe++;
        double diff = now.w - w;
        diff = diff < 0 ? -diff : diff;
        if (diff <= 1e-16) {
            return 0;
        } else {
            return w;
        }
    } else {
        ++*nnz;
        return w;
    }
}
#endif