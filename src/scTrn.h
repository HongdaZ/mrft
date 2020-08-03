#ifndef SCTRN_H
#define SCTRN_H

#include <vector>

using std::vector;

int scTrn( vector<int> &labels, vector<int> &regions,  vector<int> &front,
            vector<int> &region, const vector<int> &tumor_labels,
            const vector<int> &n_voxel, int *ptr_label, const int *ptr_nidx, 
            const int &len, int start );

#endif
