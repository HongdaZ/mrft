#include "clearVector.h"

// erase all elements from a container
template <typename T >
void clearVector( T &v ) {
  while( v.size() > 0 ) {
    v.pop_back();
  }
  return;
}