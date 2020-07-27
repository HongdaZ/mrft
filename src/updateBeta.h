#ifndef UPDATEBETA_H
#define UPDATEBETA_H

#include <vector>
#include <map>
#include <list>

using std::vector;
using std::map;
using std::list;

// update beta for healthy and tumor regions
void updateBeta( double *ptr_beta, const double *ptr_alpha, 
                 map<int, vector<double>> &health_parm,
                 map<int, list<int>> &tumor_regions,
                 map<int, vector<double>> &tumor_parm );

#endif