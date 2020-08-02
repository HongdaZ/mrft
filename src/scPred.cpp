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
int scPred( vector<int> &labels, vector<int> &regions,  vector<int> &front,
            vector<int> &region, const vector<int> &tumor_labels,
            const vector<int> &n_voxel, int *ptr_label, const int *ptr_nidx,
            const int &len, int start ) {
  clearVector( labels );
  zeroVector( regions );
  int current = ptr_label[ 2 * ( start - 1 ) ];
  if( current <= -4 ) {
    labels.push_back( current );
  }
  
  list<int> tumor_nbr;
  list<int> tumor_label;
  // Add tumor neighbors to list
  tumorNbrLabel( tumor_label, tumor_nbr, start, ptr_label, ptr_nidx );
  if( tumor_nbr.size() < 2 ) {
    return 0;
    // possibly split
  } else if( current <= - 4 ) {
    // remember to reset seg[ 2, ] == 0
    ptr_label[ 2 * start - 1  ] = 1;
    // remember to reset seg[ 2, ] == 0
    for( list<int>::iterator it = tumor_nbr.begin();
         it != tumor_nbr.end(); ++ it ) {
      ptr_label[ 2 * *it - 1 ] = 2;
    }
    bool early_return = 0;
    int n_region;
    while( !tumor_nbr.empty() ) {
      n_region = labels.size() - 1;
      int new_start = *tumor_nbr.begin();
      tumor_nbr.pop_front();
      ptr_label[ 2 * new_start - 1 ] = 0;
      findRegion( region, front, n_region, ptr_label, ptr_nidx, 
                  false, tumor_nbr, early_return, new_start );
      if( early_return ) {
        ptr_label[ 2 * start - 1  ] = 0;
        return 0;
      }
      int region_label;
      if( labels.size() == 1 ) {
        region_label = current;
      } else {
        region_label = newTumorLabel( regions.size(), tumor_labels );
        // region_label = ( *tumor_labels.begin() ) - regions.size();
        
      }
      labels.push_back( region_label );
      assignRegion( regions, region, 1, region_label );
    }
    // restore the value of seg[ 2, ] to 0;
    ptr_label[ 2 * start - 1  ] = 0;
    if( labels.size() < 3 ) {
      return 0;
    } else {
      getRegion( region, current, ptr_label, len );
      assignRegion( regions, region, 0, current );
      return 1;
    }
  } else {
    if( tumor_label.size() == 1 ) {
      return 0;
    }
    // possibly combine
    int combine_label = tumor_label.back();
    labels.push_back( combine_label );

    for( list<int>::iterator it = tumor_label.begin();
         it != tumor_label.end(); ++ it ) {
      getRegion( region, *it, ptr_label, len );
      assignRegion( regions, region, 1, *it );
      
      labels.push_back( *it );
      assignRegion( regions, region, 0, combine_label );
    }
    
    region.resize( 1 );
    region[ 0 ] = start;
    assignRegion( regions, region, 0, combine_label );
    
    return 2;
  }
}