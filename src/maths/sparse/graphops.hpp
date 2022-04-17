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
    double w;
    edge(int u, int v, double w) : u(u), v(v), w(w) {}
};

struct unionset {
    int *parent;
    int *rank;
    int n;
    unionset(int n) : n(n) {
        parent = new int[n + 1];
        rank = new int[n + 1];
        for (int i = 1; i <= n; i++) {
            parent[i] = i;
            rank[i] = 0;
        }
    }
    ~unionset() {
        delete[] parent;
        delete[] rank;
    }
    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }
    void merge(int x, int y) {
        x = find(x);
        y = find(y);
        if (x == y) return;
        if (rank[x] < rank[y]) {
            parent[x] = y;
        } else {
            parent[y] = x;
            if (rank[x] == rank[y]) rank[x]++;
        }
    }
};

struct graph {
    int n, m;
    std::vector<std::vector<int>> adj;
    std::vector<std::vector<double>> w;
    std::vector<edge> alle, orige;
    int nowe;
    std::vector<double> diag;
    void insEdge(edge e);
    int sparsify();
};
#else
typedef struct edge edge;
typedef struct graph graph;
#endif

#ifdef __cplusplus
extern "C" {
#endif
void initGraph(graph **g, int n);
void addEdge(graph *g, int u, int v, double w);
void addOrigEdge(graph *g, int u, int v, double w);
void deleteGraph(graph *g);
int sparsify(graph *g);
int graphToMatrix(graph *g, MatrixPtr Prec);
void clearGraph(graph *g);
void setDiag(graph *g, int Row, double diag);
double checkEdge(graph *g, int u, int v, double w, int *nnz);
#ifdef __cplusplus
}
#endif

#endif