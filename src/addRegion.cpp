#include "addRegion.h"
// add new tumor regions
void addRegion( const list<vector<double>> &region_parm, 
                const list<map<int,int>> &regions, 
                set<int> &tumor_labels,
                map<int, set<int>> &tumor_regions, 
                map<int, vector<double>> &tumor_parm ) {
  
  list<vector<double>>::const_iterator it = region_parm.begin();
  list<map<int, int >>::const_iterator it_region = regions.begin();
  for( ; it != region_parm.end(); ++ it, ++ it_region ) {
    vector<double>::const_iterator it_parm = it->begin();
    int new_label = *it_parm;
    vector<double> new_parm( ++ it_parm, it->end() );
    set<int> new_region;
    for( map<int, int>::const_iterator it_map = it_region->begin();
         it_map != it_region->end(); ++ it_map ) {
      new_region.insert( it_map->first );
    }
    tumor_labels.insert( new_label );
    tumor_regions[ new_label ] = new_region;
    tumor_parm[ new_label ] = new_parm;
  }
  return;
}