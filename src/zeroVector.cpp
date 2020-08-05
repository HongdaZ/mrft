#include "zeroVector.h"

// assign zero to all elements of the vector
template < typename T >
void zeroVector( vector< T > &v ) {
  for( int i = 0; i < v.size(); ++ i ) {
    v[ i ] = 0;
  }
  return;
}