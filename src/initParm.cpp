#include <R.h>
#include <Rinternals.h>

#include "initParm.h"
#include "updateParm.h"
#include "energyY.h"
#include "energyX.h"

// compare energy for training
void initParm( map<int, vector<double>> &health_parm,
               map<int, vector<double>> &tumor_parm,
               const int *ptr_seg, const double *ptr_m,
               const double *nu2, const double *ptr_intst, 
               const double *ptr_lambda2, const int *ptr_nidx,
               const double *ptr_nintst, const double *ptr_alpha,
               const double *ptr_beta,
               map<int, set<int>> tumor_regions, 
               double *ptr_a, double *ptr_b, int maxit ) {
  // update parameters for tumor regions
  for( map<int, set<int>>::iterator it = tumor_regions.begin();
       it != tumor_regions.end(); ++ it ) {
    int curr_label = it->first;
    set<int> set_region = it->second;
    map<int, int>  region;
    for( set<int>::iterator set_it = set_region.begin();
         set_it != set_region.end(); ++ set_it ) {
      region[ *set_it ] = curr_label;
    }
    double mu, sigma2;
    vector<double> theta[ 6 ];
    updateParm( mu, theta, sigma2, region, ptr_m[ 3 ], ptr_m[ 2 ], ptr_a[ 0 ],
                ptr_b[ 0 ], ptr_intst, curr_label, ptr_lambda2[ 3 ], ptr_seg,
                ptr_nidx, ptr_nintst, ptr_alpha[ 3 ], ptr_beta[ 3 ], maxit );
    
  }
}