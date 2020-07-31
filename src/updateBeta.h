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
                 const vector<double> &health_parm,  
                 const vector<int> &n_voxel, const int &n_tumor,
                 const vector<double> &tumor_parm );
// update beta for t1ce or flair images
void updateBeta3( double *ptr_beta, const double *ptr_alpha, 
                  const vector<double> &health_parm );

#endif