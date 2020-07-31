#include "clearVector.h"

// erase all elements from a vector
void clearVector( vector<int> &v ) {
  while( v.size() > 0 ) {
    v.pop_back();
  }
  return;
}