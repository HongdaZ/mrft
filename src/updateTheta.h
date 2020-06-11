#ifndef UPDATETHETA_H
#define UPDATETHETA_H

#include <set>
#include <vector>

using std::set;
using std::vector;

vector<double> updateTheta( set<int> &region,
                            double mu,
                            double sigma2,
                            double lambda2,
                            const int *ptr_seg,
                            const int *ptr_nidx,
                            const double *ptr_intst,
                            const double *ptr_nintst );

#endif