#ifndef GRAPHOPS_H
#define GRAPHOPS_H

#ifdef __cplusplus
extern "C" {
#endif
    #include "ngspice/spmatrix.h"
    #include "spdefs.h"
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#include <vector>
namespace old {
struct edge {
    int u, v;
    double w;
    bool sel;
    edge(int u, int v, double w) : u(u), v(v), w(w), sel(false) {}
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
    std::vector<std::vector<int>> adj, del_adj;
    std::vector<std::vector<double>> w, del_w;
    std::vector<edge> alle, orige;
    std::vector<std::vector<edge>> cande;
    int nowe;
    std::vector<double> diag, maxw, wd, swd, del_diag;
    std::vector<int> deg, list, order;
    unionset *us;
    graph(int n): n(n), m(0), nowe(0) {
        adj.resize(n + 1);
        w.resize(n + 1);
        del_adj.resize(n + 1);
        del_w.resize(n + 1);
        cande.resize(n + 1);
        deg.resize(n + 1);
        list.resize(n + 1);
        order.resize(n + 1);
        diag.resize(n + 1);
        maxw.resize(n + 1);
        wd.resize(n + 1);
        swd.resize(n + 1);
        alle.reserve(n + 1);
        orige.reserve(n + 1);
        del_diag.resize(n + 1);
        us = new unionset(n);
    }
    ~graph() {
        delete us;
    }
    void insEdge(edge e);
    int sparsify(double p);
    double findBestRatio(int sampleNum);
    void constructMST();
};

void addEdge(graph *g, int u, int v, double w);
}
#endif

#endif