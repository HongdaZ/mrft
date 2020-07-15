#include "search.h"
extern "C" {
  void search( const bool init, int *label, const int *nidx, 
               queue<int>& front, const int& l, set<int>& region,
               set<int> &tumor_nbr, bool &early_return ) {
    while( !front.empty() ) {
      int index = front.front();
      front.pop();
      int n_idx;
      for( int i = 0; i < 6; ++ i ) {
        n_idx = nidx[ 6 * ( index - 1 ) + i ];
        // Don't use R constants in recursive function, i.e. NA_INTEGER
        if( n_idx !=  NA_INTEGER ) { 
          
          if( label[ 2 * ( n_idx - 1 ) ] == l &&
              label[ 2 * n_idx - 1 ] != 1 ) {
            // if not initialize regions
            if( !init ) {
              // remove idx from tumor_nbr
              if( label[ 2 * n_idx - 1 ] == 2 ) {
                tumor_nbr.erase( n_idx );
                if( tumor_nbr.size() == 0 ) {
                  early_return = true;
                  label[ 2 * n_idx - 1 ] = 0;
                  return;
                }
              }
            }
            
            // Rprintf( "neighbor = %d \n", n_idx );
            front.push( n_idx );
            region.insert( n_idx );
            label[ 2 * n_idx - 1 ] = 1;
          }
        }
      }
    }
    return;
  }
} // extern "C"