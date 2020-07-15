#include <R.h>
#include <Rinternals.h>

#include <set>
#include <vector>

#include "scTrn.h"
#include "findRegion.h"
#include "newTumorLabel.h"

using std::set;
using std::vector;

// start = 1 to length array
// nidx also starts from 1
int scTrn( list< map<int, int>> &regions,  
           const set<int> &tumor_labels,
           map<int, set<int>> &tumor_regions,
           int *ptr_label, const int *ptr_nidx, int start ) {
  regions.clear();
  int current = ptr_label[ 2 * ( start - 1 ) ];
  if( current <= 0 && current >= -3 ) {
    return 0;
  } else {
    // possibly split
    set<int> tumor_nbr;
    set<int> nbr_label;
    // Add tumor neighbors to set
    for( int i = 0; i < 6; ++ i ) {
      int nbr_idx = ptr_nidx[ 6 * ( start - 1 ) + i ];
      if( nbr_idx != NA_INTEGER ) {
        int new_label = ptr_label[ 2 * ( nbr_idx - 1 ) ];
        if(  new_label <= - 4 ) {
          tumor_nbr.insert( nbr_idx );
          nbr_label.insert( new_label );
        }
      }
    }
    if( tumor_nbr.size() < 2 ) {
      return 0;
    } else if( current <= - 4 ) {
      // remember to reset seg[ 2, ] == 0
      ptr_label[ 2 * start - 1  ] = 1;
      // remember to reset seg[ 2, ] == 0
      for( set<int>::iterator it = tumor_nbr.begin();
           it != tumor_nbr.end(); ++ it ) {
        ptr_label[ 2 * *it - 1 ] = 2;
      }
      bool early_return = 0;

      while( !tumor_nbr.empty() ) {
        
        int new_start = *tumor_nbr.begin();
        tumor_nbr.erase( new_start );
        ptr_label[ 2 * new_start - 1 ] = 0;
        set<int> contaguous_region = findRegion( ptr_label, ptr_nidx, false,
                                                 tumor_nbr, early_return,
                                                 new_start );
        if( early_return ) {
          ptr_label[ 2 * start - 1  ] = 0;
          return 0;
        }
        map<int, int> sub_region;
        int region_label;
        if( regions.empty() ) {
          region_label = current;
        } else {
          region_label = newTumorLabel( regions.size(), tumor_labels );
          // region_label = ( *tumor_labels.begin() ) - regions.size();
          
        }
        for( set<int>::iterator it = contaguous_region.begin();
             it != contaguous_region.end(); ++ it ) {
          sub_region[ *it ] = region_label;
        }
        // sub_region[ start ] = NA_INTEGER;
        regions.push_back( sub_region );
      }
      // restore the value of seg[ 2, ] to 0;
      ptr_label[ 2 * start - 1  ] = 0;
      if( regions.size() < 2 ) {
        regions.clear();
        return 0;
      } else {
        map<int, int> whole;
        set<int> current_tumor = tumor_regions[ current ];
        for( set<int>::iterator it = current_tumor.begin();
             it != current_tumor.end(); ++ it ) {
          whole[ *it ] = current;
    
        }
        whole[ start ] = current;
        regions.push_front( whole );
        return 1;
      }
    } else {
      map<int, int> sub_region;
      map<int, int> whole;
      // possibly combine 
      set<int>::iterator it;
      set<int> contaguous_region;
      for( it = nbr_label.begin();
           it != nbr_label.end(); ++ it ) {
        sub_region.clear();
        // // debug
        // Rprintf( "neighbor label: %d \n", *it );
        // ///////////////////////////////////////
        contaguous_region = tumor_regions[ *it ];
        for( set<int>::iterator itr = contaguous_region.begin();
             itr != contaguous_region.end(); ++ itr ) {
          sub_region[ *itr ] = *it;
        }
        // sub_region[ start ] = NA_INTEGER;
        regions.push_back( sub_region );
        whole.insert( sub_region.begin(), sub_region.end() );
      }
      int combine_label = *( -- it );
      for( map<int, int>::iterator itr = whole.begin(); 
           itr!= whole.end(); ++ itr) {
        itr->second = combine_label;
      }
      whole[ start ] = combine_label;
      regions.push_front( whole );
      return 2;
    }
  }
}