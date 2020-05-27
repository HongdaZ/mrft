#include "search.h"
extern "C" {
  void search( int *label, int *nidx, 
               queue<int>& front, const int& l, stack<int>& region ) {
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
            // Rprintf( "neighbor = %d \n", n_idx );
            front.push( n_idx );
            region.push( n_idx );
            label[ 2 * n_idx - 1 ] = 1;
          }
        }
      }
    }
    return;
  }
} // extern "C"