#include "addRegion.h"
// add new tumor regions
void addRegion( int *ptr_seg, 
                const list<vector<double>> &region_parm, 
                const list<list<int>> &regions, 
                vector<int> &tumor_labels,
                map<int, list<int>> &tumor_regions, 
                map<int, vector<double>> &tumor_parm,
                int &n_tumor ) {
  
  list<vector<double>>::const_iterator it = region_parm.begin();
  list<list<int>>::const_iterator it_region = regions.begin();
  for( ; it != region_parm.end(); ++ it, ++ it_region ) {
    vector<double>::const_iterator it_parm = it->begin();
    int new_label = *it_parm;
    vector<double> new_parm( ++ it_parm, it->end() );
    for( list<int>::const_iterator it_ = it_region->begin();
         it_ != it_region->end(); ++ it_ ) {
      // update ptr_seg
      ptr_seg[ 2 * ( *it_ - 1 ) ] = new_label;
    }
    tumor_labels[ - new_label - 4 ] = 1;
    tumor_regions[ new_label ] = *it_region;
    tumor_parm[ new_label ] = new_parm;
    ++ n_tumor;
  }
  return;
}