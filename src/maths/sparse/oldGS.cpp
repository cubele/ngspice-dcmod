#include "feGRASS.hpp"
#include "graphops.hpp"

int matGraph::sparsify_old(double p) {
    init();
    old::graph g = old::graph(n);
    for (auto e : edges) {
        g.addEdge(e.u, e.v, e.w);
    }
    int remm = g.sparsify(p);
    int cnt = 0;
    for (int i = 1; i <= n; i++) {
        for (int j = 0; j < g.del_adj[i].size(); ++j) {
            int v = g.del_adj[i][j];
            double w = g.del_w[i][j];
            del_map[i][v] = w;
            del_diag[i] -= w;
            ++cnt;
        }
    }
    printf("%d edges removed\n", cnt / 2);
    return remm;
}