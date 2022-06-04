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
    double getRatio() {
        int n = trials.size();
        if (n < 3) {
            return -1;
        }
        //calculate the quadratic function using the last three trials
        double x[3] = {trials[n-1].ratio, trials[n-2].ratio, trials[n-3].ratio};
        double y[3] = {trials[n-1].lutime + trials[n - 1].gmrestime, trials[n-2].lutime + trials[n - 2].gmrestime, trials[n-3].lutime + trials[n - 3].gmrestime};
        double a[3] = {0, 0, 0};
        for (int i = 0; i < 3; i++) {
            printf("ratio: %f time: %f\n", x[i], y[i]);
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
        //calculate the ratio based on the smallest value of the quadratic function
        double ratio = -a[1] / (2*a[2]);
        printf("final ratio: %f\n", ratio);
        double minratio = std::min(x[0], std::min(x[1], x[2]));
        ratio = std::max(ratio, minratio);
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