# ifndef _fegrass_h_
# define _fegrass_h_

#ifdef __cplusplus
# include <iostream>
# include <fstream>
# include <vector>
# include <algorithm>
# include <cmath>
# include <stdlib.h>
# include <cstring>
# include <queue>
# include <ctime>
# include <string>
# include <list>
# include <stack>
# include <unordered_map>

#ifdef __cplusplus
extern "C" {
#endif
    #include "ngspice/smpdefs.h"
    #include "ngspice/ngspice.h"
#ifdef __cplusplus
}
#endif

void feGRASS(int* ai_in,int* aj_in,double* av_in,int M_in,int N_in,int* insparsifier_in,double alpha_in);

struct edge {
    int u, v;
    double w;
    bool in;
    edge(int u, int v, double w) : u(u), v(v), w(w), in(false) {}
};
struct matGraph {
    int n, now;
    std::vector<std::vector<int>> del_adj;
    std::vector<std::vector<double>> del_w;
    std::vector<double> del_diag;
    std::vector<edge> edges;
    std::vector<std::unordered_map<int, double>> del_map;
    matGraph(int n): n(n) {
        del_adj.resize(n + 1);
        del_w.resize(n + 1);
        del_diag.resize(n + 1);
        del_map.resize(n + 1);
    }
    void init() {
        for (int i = 1; i <= n; ++i) {
            del_adj[i].clear();
            del_w[i].clear();
            del_diag[i] = 0;
            del_map[i].clear();
        }
        now = 0;
    }
    void clear() {
        edges.clear();
    }
    int sparsify(double p) {
        init();
        int m = edges.size();
        int *u = new int[m];
        int *v = new int[m];
        double *w = new double[m];
        int *in = new int[m];
        for (int i = 0; i < m; i++) {
            u[i] = edges[i].u - 1;
            v[i] = edges[i].v - 1;
            w[i] = fabs(edges[i].w);
        }
        feGRASS(u, v, w, m, n, in, p);
        int remm = 0;
        for (int i = 0; i < m; i++) {
            if (in[i] == 0) {
                del_map[edges[i].u][edges[i].v] = edges[i].w;
                del_map[edges[i].v][edges[i].u] = edges[i].w;
                /*
                del_adj[edges[i].u].push_back(edges[i].v);
                del_adj[edges[i].v].push_back(edges[i].u);
                del_w[edges[i].u].push_back(edges[i].w);
                del_w[edges[i].v].push_back(edges[i].w);
                */
                del_diag[edges[i].u] -= edges[i].w;
                del_diag[edges[i].v] -= edges[i].w;
            } else {
                remm++;
            }
        }
        delete[] u;
        delete[] v;
        delete[] w;
        delete[] in;
        printf("%d edges removed\n", m - remm);
        return remm;
    }
    void addEdge(int u, int v, double w) {
        if (u < v) {
            std::swap(u, v);
        }
        edges.push_back(edge(u, v, w));
    }
    double checkEdge(int u, int v, double w, int *isedge) {
        if (u == v) {
            return w - del_diag[u];
        } else {
            if (del_map[u].find(v) != del_map[u].end()) {
                *isedge = 1;
                w -= del_map[u][v];
            }
            return fabs(w) < 1e-20 ? 0.0 : w;
        }
    }
};
#else
typedef struct edge edge;
typedef struct matGraph matGraph;
#endif

#ifdef __cplusplus
extern "C" {
#endif
double checkEdge(matGraph *g, int u, int v, double w, int *isedge);
void initGraph(matGraph **g, int n);
void addEdge(matGraph *g, int u, int v, double w);
int sparsify(matGraph *g, double p);
double findRatio(matGraph *g);
void clearGraph(matGraph *g);
#ifdef __cplusplus
}
#endif

# endif