#include <R.h>
#include <Rinternals.h>
#include <stack>
#include <queue>
using std::queue;
using std::stack;

stack<int>& search( int *label, int *nidx, 
                    queue<int>& front, const int & l, stack<int>& region ) {
  if( !front.empty() ) {
    int index = front.front();
    front.pop();
    if( label[ 2 * ( index - 1 ) ] == l & label[ 2 * index - 1 ] != 1 ) {
      label[ 2 * index - 1 ] = 1;
      region.push( index );
      for( int i = 0; i < 6; ++ i ) {
        if( !ISNA( nidx[ 6 * ( index - 1 ) + i ] ) ) {
          front.push( nidx[ 6 * ( index - 1 ) + i ] );
        }
      }
    }
    return search( label, nidx, front, l, region );
  } else {
    return region;
  }
  
}