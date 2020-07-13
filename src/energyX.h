#include <R.h>
#include <Rinternals.h>

#include <vector>

using std::vector;



// region starts from 1
// calculate energyX for single voxel
double energyX( const int curr_label, const int curr_idx, bool region,
                const int *ptr_seg, const int *ptr_nidx,
                double delta, double gamma );
