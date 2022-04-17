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
    (*g)->cande.resize(n + 1);
    (*g)->w.resize(n + 1);
    (*g)->diag.resize(n + 1);
    (*g)->maxw.resize(n + 1);
    (*g)->wd.resize(n + 1);
    (*g)->deg.resize(n + 1);
    (*g)->list.resize(n + 1);
    (*g)->order.resize(n + 1);
    (*g)->swd.resize(n + 1);
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
    adj[e.v].push_back(e.u);
    w[e.v].push_back(e.w);
    maxw[e.u] = std::max(maxw[e.u], fabs(e.w));
    maxw[e.v] = std::max(maxw[e.v], fabs(e.w));
    diag[e.u] += fabs(e.w), diag[e.v] += fabs(e.w);
    deg[e.u]++, deg[e.v]++;
    swd[e.u] = diag[e.u] / maxw[e.u];
    swd[e.v] = diag[e.v] / maxw[e.v];
    ++m;
}

//edges are with negative weights
void graph::constructMST() {
    std::sort(alle.begin(), alle.end(), [](edge a, edge b) {
        return -a.w > -b.w;
    });
    unionset *us = new unionset(n);
    for (int i = 0; i < alle.size(); i++) {
        if (us->find(alle[i].u) != us->find(alle[i].v)) {
            insEdge(alle[i]);
            us->merge(alle[i].u, alle[i].v);
            alle[i].u = -1, alle[i].v = -1; //selected
        }
    }
    delete us;
}

int graph::sparsify(double p) {
    printf("edges berore sparsify: %d\n", alle.size());
    for (int i = 1; i <= n; ++i) {
        diag[i] = 0, maxw[i] = 0, wd[i] = 0, deg[i] = 0;
        list[i] = i;
    }
    for (int i = 0; i < alle.size(); i++) {
        diag[alle[i].u] += fabs(alle[i].w);
        diag[alle[i].v] += fabs(alle[i].w);
        deg[alle[i].u]++, deg[alle[i].v]++;
        maxw[alle[i].u] = std::max(maxw[alle[i].u], fabs(alle[i].w));
        maxw[alle[i].v] = std::max(maxw[alle[i].v], fabs(alle[i].w));
    }
    for (int i = 1; i <= n; ++i) {
        if (deg[i] > 0) {
            wd[i] = diag[i] / maxw[i];
        } else {
            wd[i] = -1;
        }
    }
    std:sort(list.begin() + 1, list.end(), [this](int a, int b) {
        return wd[a] > wd[b];
    });
    for (int i = 1; i <= n; ++i)  {
        order[list[i]] = i;
    }
    for (int i = 1; i <= n; ++i) {
        diag[i] = 0, maxw[i] = 0, deg[i] = 0, swd[i] = 0;
    }
    constructMST();
    //edges not in MST
    for (int i = 0; i < alle.size(); i++) {
        if (alle[i].u != -1) {
            if (order[alle[i].u] > order[alle[i].v]) {
                std::swap(alle[i].u, alle[i].v);
            }
            cande[alle[i].u].push_back(alle[i]);
        }
    }
    for (int i = 1; i <= n; ++i) {
        int u = list[i];
        std::sort(cande[u].begin(), cande[u].end(), [](edge a, edge b) {
            return -a.w > -b.w;
        });
        for (int j = 0; j < cande[u].size(); j++) {
            if (swd[u] < p * wd[u]) {
                insEdge(cande[u][j]);
            } else {
                //edges not selected may be selected in the next round
                std::swap(cande[u][j].u, cande[u][j].v);
                cande[cande[u][j].v].push_back(cande[u][j]);
            }
        }
    }
    printf("edges after sparsify: %d n: %d\n", m, n);
    return m;
}

int sparsify(graph *g, double p) {
    return g->sparsify(p);
}

void clearGraph(graph *g) {
    for (int i = 1; i <= g->n; ++i) {
        g->adj[i].clear();
        g->w[i].clear();
        g->cande[i].clear();
        g->diag[i] = 0;
    }
    g->alle.clear();
    g->orige.clear();
    g->nowe = 0;
    g->m = 0;
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
            nnz += 1;
            SMPaddElt(Prec, i, v, w);
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