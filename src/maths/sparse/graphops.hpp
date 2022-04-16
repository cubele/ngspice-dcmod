#ifndef GRAPHOPS_H
#define GRAPHOPS_H

#ifdef __cplusplus
#include <vector>
struct edge {
    int u, v;
    double w;
    edge(int u, int v, double w) : u(u), v(v), w(w) {}
};

struct graph {
    int n, m;
    std::vector<std::vector<edge>> adj;
    std::vector<edge> e;
};
#else
typedef struct edge edge;
typedef struct graph graph;
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif