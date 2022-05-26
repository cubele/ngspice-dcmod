#include <algorithm>
#include <iostream>
#include <cmath>
#include "graphops.hpp"
//ngspice header files written in C
#ifdef __cplusplus
extern "C" {
#endif
    #include "ngspice/smpdefs.h"
    #include "ngspice/ngspice.h"
#ifdef __cplusplus
}
#endif

using namespace old;

void initGraph(graph **g, int n) {
    *g = new graph(n);
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
    for (int i = 1; i <= n; ++i) {
        us->parent[i] = i;
        us->rank[i] = 0;
    }
    for (int i = 0; i < alle.size(); i++) {
        if (us->find(alle[i].u) != us->find(alle[i].v)) {
            insEdge(alle[i]);
            us->merge(alle[i].u, alle[i].v);
            alle[i].sel = true;
        }
    }
}

int graph::sparsify(double p) {
    printf("sparsifying graph\n");
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
        alle[i].sel = false;
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
        diag[i] = 0, maxw[i] = 0, deg[i] = 0, swd[i] = 0, del_diag[i] = 0;
    }
    constructMST();
    //edges not in MST
    for (int i = 0; i < alle.size(); i++) {
        if (!alle[i].sel) {
            if (order[alle[i].u] > order[alle[i].v]) {
                std::swap(alle[i].u, alle[i].v);
            }
            cande[alle[i].u].push_back(alle[i]);
        }
    }
    for (int i = 1; i <= n; ++i) {
        int u = list[i], cnt = 0;
        std::sort(cande[u].begin(), cande[u].end(), [](edge a, edge b) {
            return -a.w > -b.w;
        });
        for (int j = 0; j < cande[u].size(); j++) {
            if (swd[u] < p * wd[u] && cnt < 2) {
                insEdge(cande[u][j]);
                ++cnt;
            } else {
                //edges not selected may be selected in the next round
                if (order[cande[u][j].u] < order[cande[u][j].v]) {
                    std::swap(cande[u][j].u, cande[u][j].v);
                    cande[cande[u][j].u].push_back(cande[u][j]);
                } else {
                    del_adj[cande[u][j].u].push_back(cande[u][j].v);
                    del_adj[cande[u][j].v].push_back(cande[u][j].u);
                    del_w[cande[u][j].u].push_back(cande[u][j].w);
                    del_w[cande[u][j].v].push_back(cande[u][j].w);
                    del_diag[cande[u][j].u] += -cande[u][j].w;
                    del_diag[cande[u][j].v] += -cande[u][j].w;
                }
            }
        }
    }
    printf("edges after sparsify: %d n: %d\n", m, n);
    return m;
}

int sparsify(graph *g, double p) {
    return g->sparsify(p);
}

double graph::findBestRatio(int sampleNum) {
    for (int i = 1; i <= n; ++i) {
        diag[i] = 0, maxw[i] = 0, wd[i] = 0, deg[i] = 0;
    }
    for (int i = 0; i < alle.size(); i++) {
        diag[alle[i].u] += fabs(alle[i].w);
        diag[alle[i].v] += fabs(alle[i].w);
        deg[alle[i].u]++, deg[alle[i].v]++;
        maxw[alle[i].u] = std::max(maxw[alle[i].u], fabs(alle[i].w));
        maxw[alle[i].v] = std::max(maxw[alle[i].v], fabs(alle[i].w));
    }
    double awd = 0;
    for (int i = 1; i <= n; ++i) {
        if (deg[i] > 0) {
            wd[i] = diag[i] / maxw[i];
            awd += wd[i];
        } else {
            wd[i] = -1;
            awd += 1;
        }
    }
    awd /= n;
    double minr = 1 / awd, maxr = 1, ratio = minr;
    double maxdif = 0;
    std::vector<double> nnz(sampleNum + 1);
    for (int i = 1; i <= sampleNum; ++i) {
        double r = minr + (maxr - minr) * (i - 1) / (sampleNum - 1);
        nnz[i] = sparsify(r);
        nnz[i] = nnz[i] * nnz[i];
        for (int j = 1; j <= n; ++j) {
            adj[j].clear();
            w[j].clear();
            cande[j].clear();
            del_adj[j].clear();
            del_w[j].clear();
            m = 0;
        }
        if (i > 1 && nnz[i] - nnz[i - 1] > maxdif) {
            ratio = minr + (maxr - minr) * (i - 2) / (sampleNum - 1);
            maxdif = nnz[i] - nnz[i - 1];
        }
    }
    printf("\n");
    printf("%.6lf %.6lf %.6lf\n", minr, ratio, maxr);
    return ratio;
}

double findBestRatio(graph *g, int sampleNum) {
    return g->findBestRatio(sampleNum);
}

void clearGraph(graph *g) {
    for (int i = 1; i <= g->n; ++i) {
        g->adj[i].clear();
        g->w[i].clear();
        g->cande[i].clear();
        g->del_adj[i].clear();
        g->del_w[i].clear();
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
    printf("graphToMatrix\n");
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
            printf("%d %.3lf\n", i, -g->diag[i]);
            ++nnz;
        }
    }
    printf("nnz: %d\n", nnz);
    fflush(stdout);
    return nnz;
}

double checkEdge(graph *g, int u, int v, double w, int *nnz) {
    if (u == v) {
        return w - g->del_diag[u];
    } else {
        for (int i = 0; i < g->del_adj[u].size(); i++) {
            if (g->del_adj[u][i] == v) {
                w -= g->del_w[u][i];
            }
        }
        return fabs(w) < 1e-10 ? 0 : w;
    }
}