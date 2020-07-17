#include <R.h>
#include <Rinternals.h>

#include <list>
#include <vector>
#include <algorithm>

#include "scTrn.h"
#include "findRegion.h"
#include "newTumorLabel.h"
#include "nbrLabel.h"

using std::list;
using std::vector;
using std::find;
// start = 1 to length array
// nidx also starts from 1
int scTrn( list<int> &labels, list<list<int>> &regions,  
           const list<int> &tumor_labels,
           map<int, list<int>> &tumor_regions,
           int *ptr_label, const int *ptr_nidx, int start ) {
  labels.clear();
  regions.clear();
  int current = ptr_label[ 2 * ( start - 1 ) ];
  if( current <= 0 && current >= -3 ) {
    return 0;
  } else {
    
    list<int> tumor_nbr;
    list<int> nbr_label;
    // Add tumor neighbors to set
    nbrLabel( nbr_label, tumor_nbr, start, ptr_label, ptr_nidx );
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
        n_region = regions.size();
        int new_start = *tumor_nbr.begin();
        tumor_nbr.pop_front();
        ptr_label[ 2 * new_start - 1 ] = 0;
        list<int> sub_region = findRegion( n_region, ptr_label, ptr_nidx, false,
                                                 tumor_nbr, early_return,
                                                 new_start );
        if( early_return ) {
          ptr_label[ 2 * start - 1  ] = 0;
          return 0;
        }
        int region_label;
        if( regions.empty() ) {
          region_label = current;
        } else {
          region_label = newTumorLabel( regions.size(), tumor_labels );
          // region_label = ( *tumor_labels.begin() ) - regions.size();
          
        }
        labels.push_back( region_label );
        regions.push_back( sub_region );
      }
      // restore the value of seg[ 2, ] to 0;
      ptr_label[ 2 * start - 1  ] = 0;
      if( regions.size() < 2 ) {
        regions.clear();
        return 0;
      } else {
        list<int> current_tumor = tumor_regions[ current ];
        labels.push_front( current );
        regions.push_front( current_tumor );
        return 1;
      }
    } else {
      if( nbr_label.size() == 1 ) {
        return 0;
      }
      list<int> whole;
      // possibly combine 
      list<int>::iterator it;
      list<int> sub_region;
      for( it = nbr_label.begin();
           it != nbr_label.end(); ++ it ) {
        sub_region = tumor_regions[ *it ];
        
        labels.push_back( *it );
        regions.push_back( sub_region );
        
        whole.insert( whole.end(), sub_region.begin(), sub_region.end() );
      }
      int combine_label = *( -- it );
      whole.push_front( start );
      labels.push_front( combine_label );
      regions.push_front( whole );
      return 2;
    }
  }
}