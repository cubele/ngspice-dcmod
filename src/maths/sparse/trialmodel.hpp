# ifndef TRIALMODEL_HPP
# define TRIALMODEL_HPP

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

#ifdef __cplusplus
extern "C" {
#endif
    #include "ngspice/smpdefs.h"
    #include "ngspice/ngspice.h"
#ifdef __cplusplus
}
#endif

struct trial {
    double ratio, lutime, gmrestime;
    trial(double ratio, double lutime, double gmrestime) : ratio(ratio), lutime(lutime), gmrestime(gmrestime) {}
};

struct trialModel {
    std::vector<trial> trials;
    trialModel() {}
    void addTrial(double ratio, double lutime, double gmrestime) {
        trials.push_back(trial(ratio, lutime, gmrestime));
    }
    void calCurve(double *x, double *y, double *a) {
        for (int i = 0; i < 3; i++) {
            double d = 1;
            double t[3] = {1, 0, 0};
            for (int j = 0; j < 3; ++j) {
                if (i != j) {
                    d = d * (x[i] - x[j]);
                    for (int k = 2; k >= 0; --k) {
                        t[k] = t[k] * -x[j];
                        if (k != 0) {
                            t[k] += t[k - 1];
                        }
                    }
                }
            }
            for (int j = 0; j < 3; ++j) {
                t[j] = t[j] * y[i] / d;
                a[j] += t[j];
            }
        }
    }
    double getRatio() {
        int n = trials.size();
        if (n != 4) {
            return -1;
        }
        //calculate the quadratic function using the last three trials
        double x[3], y[3], a[3];
        double ratio, mineps = 1e9;
        for (int i = 0; i < n; ++i) {
            printf("ratio:%f lutime:%f GMREStime:%f totaltime:%f\n", trials[i].ratio, trials[i].lutime, trials[i].gmrestime, trials[i].lutime + trials[i].gmrestime);
            int now = 0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    x[now] = trials[j].ratio;
                    y[now] = trials[j].lutime + trials[j].gmrestime;
                    now++;
                }
            }
            a[0] = 0, a[1] = 0, a[2] = 0;
            calCurve(x, y, a);
            double x0 = trials[i].ratio, y0 = trials[i].lutime + trials[i].gmrestime;
            double eps = fabs(a[2] * x0 * x0 + a[1] * x0 + a[0] - y0);
            printf("eps:%f\n", eps);
            if (eps < mineps) {
                mineps = eps;
                ratio = -a[1] / (2 * a[2]);
            }
        }
        printf("final ratio: %f\n", ratio);
        ratio = std::max(ratio, x[0]);
        return ratio;
    }
};
#else
typedef struct trialModel trialModel;
#endif

#ifdef __cplusplus
extern "C" {
#endif
void addTrial(trialModel *T, double ratio, double lutime, double gmrestime);
void initTrialModel(trialModel **t, int n);
double getRatio(trialModel *T);
#ifdef __cplusplus
}
#endif

#endif