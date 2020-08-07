#include <R.h>
#include <Rinternals.h>

#include <list>
#include <vector>

#include "scPred.h"
#include "findRegion.h"
#include "newTumorLabel.h"
#include "tumorNbr.h"
#include "nonTumor.h"
#include "clearVector.h"
#include "zeroVector.h"
#include "assignRegion.h"
#include "getRegion.h"

using std::list;
using std::vector;

// start = 1 to length array
// nidx also starts from 1

// split or combine for prediction
int scPred( int &n_region, list<list<int>> &tumor_regions, 
            vector<int> &front, vector<int> &region, 
            const vector<int> &tumor_labels,
            int *ptr_label, const int *ptr_nidx, 
            const int &len, const int &start, 
            vector<int> &regions_whole, vector<int> &regions_sub ) {
  n_region = 0;
  clearVector( regions_whole );
  clearVector( regions_sub );

  int current = ptr_label[ 2 * ( start - 1 ) ];
  list<int> tumor_nbr;
  list<int> tumor_label;
  // Add tumor neighbors to list
  tumorNbrLabel( tumor_label, tumor_nbr, start, ptr_label, ptr_nidx );
  if( tumor_nbr.size() < 2 ) {
    return 0;
    // possibly split
  } else if( current <= - 4 ) {
    regions_whole.push_back( current );
    regions_whole.push_back( start );
    // remember to reset seg[ 2, ] == 0
    ptr_label[ 2 * start - 1  ] = 1;
    // remember to reset seg[ 2, ] == 0
    for( list<int>::iterator it = tumor_nbr.begin();
         it != tumor_nbr.end(); ++ it ) {
      ptr_label[ 2 * *it - 1 ] = 2;
    }
    bool early_return = 0;
    // find sub-regions
    while( !tumor_nbr.empty() ) {
      int new_start = tumor_nbr.front();
      tumor_nbr.pop_front();
      ptr_label[ 2 * new_start - 1 ] = 0;
      findRegion( region, front, n_region, ptr_label, ptr_nidx, 
                  false, tumor_nbr, early_return, new_start );
      if( early_return ) {
        ptr_label[ 2 * start - 1  ] = 0;
        return 0;
      }
      int region_label;
      if( n_region == 0 ) {
        region_label = current;
      } else {
        region_label = newTumorLabel( n_region, tumor_labels );
        
      }
      assignRegion( regions_whole, regions_sub, region, region_label );
      ++ n_region;
    }
   
    regions_whole.push_back( NA_INTEGER );
    // restore the value of seg[ 2, ] to 0;
    ptr_label[ 2 * start - 1  ] = 0;
    if( n_region < 2 ) {
      return 0;
    } else {
      // count whole region
      ++ n_region;
      return 1;
    }
  } else {
    if( tumor_label.size() == 1 ) {
      return 0;
    }
    // possibly combine
    int combine_label = tumor_label.back();
    regions_whole.push_back( combine_label );
    regions_whole.push_back( start );
    // add whole region
    ++ n_region;
    // add sub-regions
    n_region += tumor_label.size();
    
    for( list<int>::iterator it = tumor_label.begin();
         it != tumor_label.end(); ++ it ) {
      assignRegion( regions_whole, regions_sub, tumor_regions,
                    *it );
    }
    regions_whole.push_back( NA_INTEGER );
    return 2;
  }
}