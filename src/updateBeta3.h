#ifndef UPDATEBETA3_H
#define UPDATEBETA3_H

#include <vector>
#include <map>

using std::vector;
using std::map;

// update beta for t1ce or flair images
void updateBeta3( double *ptr_beta, const double *ptr_alpha, 
                 map<int, vector<double>> &health_parm );

#endif