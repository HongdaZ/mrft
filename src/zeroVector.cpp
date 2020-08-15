#include "zeroVector.h"

// assign zero to all elements of the vector
void zeroVector( vector<double> &v ) {
  for( vector<double>::iterator it = v.begin();
       it != v.end(); ++ it ) {
    *it = 0;
  }
  return;
}