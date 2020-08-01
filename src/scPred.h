#ifndef SCPRED_H
#define SCPRED_H

#include <list>
#include <map>
#include <set>
#include <vector>

using std::list;
using std::map;
using std::set;
using std::vector;

// split or combine for prediction
int scPred( vector<int> &labels, vector<int> &regions,  vector<int> &front,
            vector<int> &region, const vector<int> &tumor_labels,
            const vector<int> &n_voxel, int *ptr_label, const int *ptr_nidx, 
            const int &len, int start );

#endif