#ifndef UPDATEBETA_H
#define UPDATEBETA_H

#include <vector>
#include <map>
#include <set>

using std::vector;
using std::map;
using std::set;

// add new outlier
void updateBeta( double *ptr_beta, const double *ptr_alpha, 
                 map<int, vector<double>> &health_parm,
                 map<int, set<int>> &tumor_regions,
                 map<int, vector<double>> &tumor_parm );

#endif