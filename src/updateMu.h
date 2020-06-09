#ifndef UPDATEMU_H
#define UPDATEMU_H

#include <set>
#include <vector>

using std::set;
using std::vector;

double updateMu( set<int> &region,  
                 double sigma2,
                 double m,
                 double nu2,
                 vector<double> theta,
                 double *ptr_intst );

#endif