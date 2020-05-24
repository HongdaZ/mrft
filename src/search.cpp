#include <R.h>
#include <Rinternals.h>
#include <stack>
#include <queue>
using std::queue;
using std::stack;
using std::size_t;
void search( int *label, int *nidx, 
             queue<int>& front, const int & l, stack<int>& region ) {
  if( !front.empty() ) {
    int index = front.front();
    front.pop();
    int n_idx;
    for( int i = 0; i < 6; ++ i ) {
      n_idx = nidx[ 6 * ( index - 1 ) + i ];
      // Rprintf( "neighbor = %d \n", n_idx );
      if( !ISNA( n_idx ) ) { // have to include an R internal function
        if( label[ 2 * ( n_idx - 1 ) ] == l &&
            label[ 2 * n_idx - 1 ] != 1 ) {
          front.push( n_idx );
          region.push( n_idx );
          label[ 2 * n_idx - 1 ] = 1;
        }
      }
    }
    return search( label, nidx, front, l, region );
  } else {
    return;
  }
  
}