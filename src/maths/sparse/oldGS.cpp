#include "feGRASS.hpp"
#include "graphops.hpp"

int matGraph::sparsify_old(double p) {
    old::graph g = old::graph(n);
    for (auto e : edges) {
        printf("%d %d %f\n", e.u, e.v, e.w);
        old::addEdge(&g, e.u, e.v, e.w);
    }
    int remm = g.sparsify(p);
    for (int i = 1; i <= n; i++) {
        for (int j = 0; j < del_adj[i].size(); ++j) {
            int v = del_adj[i][j];
            double w = del_w[i][j];
            del_map[i][v] = w;
            del_map[v][i] = w;
            del_diag[i] -= w;
            del_diag[v] -= w;
        }
    }
    return remm;
}