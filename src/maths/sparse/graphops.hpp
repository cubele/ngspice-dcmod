#ifndef GRAPHOPS_H
#define GRAPHOPS_H

#ifdef __cplusplus
extern "C" {
#endif
    #include "ngspice/spmatrix.h"
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#include <vector>
struct edge {
    int u, v;
    int ou, ov; //original outer row and col from matrix
    double w;
    edge(int u, int v, int ou, int ov, double w) : u(u), v(v), ou(ou), ov(ov), w(w) {}
};

struct graph {
    int n, m;
    std::vector<std::vector<int>> adj;
    std::vector<std::vector<double>> w;
    std::vector<edge> alle;
    int nowe;
    std::vector<double> diag;
    void insEdge(edge e) {
        adj[e.u].push_back(e.v);
        w[e.u].push_back(e.w);
        adj[e.v].push_back(e.u);
        w[e.v].push_back(e.w);
    }
};
#else
typedef struct edge edge;
typedef struct graph graph;
#endif

#ifdef __cplusplus
extern "C" {
#endif
void initGraph(graph **g, int n);
void addEdge(graph *g, int u, int v, int ou, int ov, double w);
void deleteGraph(graph *g);
void sparsify(graph *g);
void graphToMatrix(graph *g, MatrixPtr Prec);
void clearGraph(graph *g);
void setDiag(graph *g, int Row, double diag);
double checkEdge(graph *g, int u, int v, double w);
#ifdef __cplusplus
}
#endif

#endif