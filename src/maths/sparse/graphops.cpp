#include <algorithm>
#include <iostream>
#include <cmath>
#include "graphops.hpp"
using namespace old;

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
    std::vector<edge> tempe;
    for (auto e : alle) {
        tempe.push_back(e);
    }

    std::sort(tempe.begin(), tempe.end(), [](edge a, edge b) {
        return -a.w > -b.w;
    });
    for (int i = 1; i <= n; ++i) {
        us->parent[i] = i;
        us->rank[i] = 0;
    }
    for (int i = 0; i < tempe.size(); i++) {
        if (us->find(tempe[i].u) != us->find(tempe[i].v)) {
            insEdge(tempe[i]);
            us->merge(tempe[i].u, tempe[i].v);
            alle[tempe[i].id].sel = true;
        }
    }
}

int graph::sparsify(double p) {
    printf("sparsifying graph\n");
    int remm = alle.size();
    printf("%d edges\n", remm);
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
    double awd = 0;
    for (int i = 1; i <= n; ++i) {
        if (deg[i] > 0) {
            wd[i] = diag[i] / maxw[i];
            awd += wd[i];
        } else {
            wd[i] = -1;
        }
    }
    awd = awd / n;
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
    int tcnt = 0;
    std::vector<bool> nowsel;
    for (int i = 0; i < alle.size(); i++) {
        if (alle[i].sel) {
            nowsel.push_back(true);
            ++tcnt;
        } else {
            nowsel.push_back(false);
        }
    }
    
    printf("%d edges in MST\n", tcnt);
    int target = (int)(p * n);
    double y = 1 / awd;
    double ratioinc = (1 - y) / 20;
    int rcnt = 0;
    while (target > 0) {
        for (int i = 1; i <= n; ++i) {
            cande[i].clear();
        }
        for (int i = 0; i < alle.size(); i++) {
            if (!nowsel[i]) {
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
                if (swd[u] < y * wd[u] && cnt < 1) {
                    insEdge(cande[u][j]);
                    nowsel[cande[u][j].id] = true;
                    ++rcnt;
                    ++cnt, --target;
                } else {
                    //edges not selected may be selected in the next round
                    if (order[cande[u][j].u] < order[cande[u][j].v]) {
                        std::swap(cande[u][j].u, cande[u][j].v);
                        cande[cande[u][j].u].push_back(cande[u][j]);
                    }
                }
            }
            if (target < 0) {
                break;
            }
        }
        y += ratioinc;
    }

    for (int i = 0; i < alle.size(); ++i) {
        if (!nowsel[i]) {
            del_adj[alle[i].u].push_back(alle[i].v);
            del_adj[alle[i].v].push_back(alle[i].u);
            del_w[alle[i].u].push_back(alle[i].w);
            del_w[alle[i].v].push_back(alle[i].w);
            del_diag[alle[i].u] -= alle[i].w;
            del_diag[alle[i].v] -= alle[i].w;
        }
    }
    printf("%d edges recovered\n", rcnt);
    printf("edges after sparsify: %d n: %d\n", m, n);
    return m;
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