#include <cmath>

#include "radius.h"

double radius( const int &n ) {
  const double PI = 3.141592653589793;
  double radius = cbrt( ( double )n * 3 / 4 / PI  );
  return radius;
}