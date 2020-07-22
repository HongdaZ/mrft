#ifndef COPYPARM_H
#define COPYPARM_H

#include <map>
#include <vector>

using std::map;
using std::vector;

// copy estimated parameters to SEXP
void copyParm( map<int, vector<double>> &health_parm,
                  map<int, vector<double>> &tumor_parm,
                  map<int, vector<double>> &outl_parm,
                  double *ptr_res_parm, int nrow );

#endif