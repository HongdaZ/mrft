#ifndef __POINT__
#define __POINT__

#include <algorithm> // std::min std::max
#include <cassert> // assert macro
#include <cstddef> // std::size_t
#include <limits> // std::numeric_limits
#include <string> // std::string std::to_string

#if not defined(MULTI_PRECISION) and not defined(FLOATING_POINT_FILTER)
#define DOUBLE_PRECISION
#endif

#if defined(MULTI_PRECISION) or defined(FLOATING_POINT_FILTER)

#include "cphmpl/functions.h++" // cphmpl::width
#include "cphstl/integers.h++" // cphstl::bbbZ

#endif

template<typename T>
class point {
public:
  
  using self = point<T>;
  using coordinate = T;
  
  T x;
  T y;
  
  explicit point()
    : x(0), y(0) {
  }
  
  point(T x_coordinate, T y_coordinate)
    : x(x_coordinate), y(y_coordinate) {
  }
  
  std::string to_string() const {
    return "(" + std::to_string(x) + ", " + std::to_string(y) + ")";
  }
  
  bool operator==(point<T> const& other) const {
    return x == other.x and y == other.y;
  }
  
  bool operator!=(point<T> const& other) const {
    return not (*this == other);
  }
  
  bool operator<(point<T> const& other) const {
    return (x < other.x) or (x == other.x and y < other.y);
  }
  
  bool operator>(point<T> const& other) const {
    return other < *this;
  }
  
  friend std::ostream& operator<<(std::ostream& os, point<T> const& p) {
    return os << p.to_string();
  }
};

#if defined(MULTI_PRECISION) or defined(FLOATING_POINT_FILTER)

namespace multi_precision {

template<typename P>
bool left_turn(P const& p, P const& q, P const& r) {
  using N = std::size_t;
  using T = typename P::coordinate;
  constexpr N w = cphmpl::width<T>;
  using Z = cphstl::bbbZ<2 * w + 4>; 
  Z p_x = Z(p.x);
  Z p_y = Z(p.y);
  Z lhs = (Z(q.x) - p_x) * (Z(r.y) - p_y);
  Z rhs = (Z(r.x) - p_x) * (Z(q.y) - p_y);
  return lhs > rhs;
}

template<typename P>
bool no_turn(P const& p, P const& q, P const& r) {
  using N = std::size_t;
  using T = typename P::coordinate;
  constexpr N w = cphmpl::width<T>;
  using Z = cphstl::bbbZ<2 * w + 4>; 
  Z p_x = Z(p.x);
  Z p_y = Z(p.y);
  Z lhs = (Z(q.x) - p_x) * (Z(r.y) - p_y);
  Z rhs = (Z(r.x) - p_x) * (Z(q.y) - p_y);
  return lhs == rhs;
}

// compute the signed area of the parallelogram spanned by
// $\overrightarrow{p q}$ and $\overrightarrow{p r}$.

template<typename P>
auto signed_area(P const& p, P const& q, P const& r) {
  using N = std::size_t;
  using T = typename P::coordinate;
  constexpr N w = cphmpl::width<T>;
  using Z = cphstl::bbbZ<2 * w + 5>; 
  Z p_x = Z(p.x);
  Z p_y = Z(p.y);
  Z dx_1 = Z(q.x) - p_x;
  Z dx_2 = Z(r.x) - p_x;
  Z dy_1 = Z(q.y) - p_y;
  Z dy_2 = Z(r.y) - p_y;
  Z result = dx_1 * dy_2 - dx_2 * dy_1;
  return result;
}

// return the square of the distance between p and q

template<typename P>
auto distance_squared(P const& p, P const& q) {
  using N = std::size_t;
  using T = typename P::coordinate;
  constexpr N w = cphmpl::width<T>;
  using Z = cphstl::bbbZ<2 * w + 5>; 
  Z p_x = Z(p.x);
  Z p_y = Z(p.y);
  Z q_x = Z(q.x);
  Z q_y = Z(q.y);
  return (p_x - q_x) * (p_x - q_x) + (p_y - q_y) * (p_y - q_y);
}

template<typename P>
auto L_infty_distance(P const& p, P const& q) {
  using N = std::size_t;
  using T = typename P::coordinate;
  constexpr N w = cphmpl::width<T>;
  using Z = cphstl::bbbZ<w + 2>; 
  Z dx = cphstl::abs(Z(q.x) - Z(p.x));
  Z dy = cphstl::abs(Z(q.y) - Z(p.y));
  return std::max(dx, dy);
}

template<typename P>
auto L_1_distance(P const& p, P const& q) {
  using N = std::size_t;
  using T = typename P::coordinate;
  constexpr N w = cphmpl::width<T>;
  using Z = cphstl::bbbZ<w + 3>; 
  Z dx = cphstl::abs(Z(q.x) - Z(p.x));
  Z dy = cphstl::abs(Z(q.y) - Z(p.y));
  return dx + dy;
}

// return the orientation when moving from p to r via q
// +1: left turn
//  0: p, q, and r are collinear
// -1: right turn

template<typename P>
int orientation(P const& p, P const& q, P const& r) {
  using N = std::size_t;
  using T = typename P::coordinate;
  constexpr N w = cphmpl::width<T>;
  using Z = cphstl::bbbZ<2 * w + 4>; 
  Z p_x = Z(p.x);
  Z p_y = Z(p.y);
  Z lhs = (Z(q.x) - p_x) * (Z(r.y) - p_y);
  Z rhs = (Z(r.x) - p_x) * (Z(q.y) - p_y);
  if (lhs == rhs) {
    return 0;
  }
  return (lhs > rhs) ? +1 : -1;
}
}

#endif

#if defined(MULTI_PRECISION)

using namespace multi_precision;

#endif

#if defined(DOUBLE_PRECISION)

namespace double_precision {

template<typename P>
bool left_turn(P const& p, P const& q, P const& r) {
  
#ifdef MEASURE_TURNS
  
  ++turns;
  
#endif
  
  using T = typename P::coordinate;
  using Z = long long;
  Z mid_x = Z(q.x);
  Z mid_y = Z(q.y);
  Z dx_1 = mid_x - Z(p.x);
  Z dx_2 = Z(r.x) - mid_x;
  Z dy_1 = mid_y - Z(p.y);
  Z dy_2 = Z(r.y) - mid_y;
  
  /* We have 31 bits + 1 for the sign in int, and
   we may need (31 (+1 from minus)) * 2 (from mult) = 64 bits
   in total, but we only have 63 bits + 1 for the sign.
   If we overflow the same direction, the order is preserved. */
  bool check1 = Z(T(dx_1)) != dx_1 and Z(T(dy_2)) != dy_2;
  bool check2 = Z(T(dx_2)) != dx_2 and Z(T(dy_1)) != dy_1;
  if (check1 or check2) {
    // overflow happening
    if (check1) {
      if ((dx_1 > 0) == (dy_2 > 0)) {
        // first part overflows above
        if (check2 and ((dx_2 > 0) == (dy_1 > 0))) {
          // both overflow above
          return dx_1 * dy_2 > dx_2 * dy_1;
        }
        else {
          return true; // only first term overflows above
        }
      }
      else {
        // first term overflows below
        if (check2 and ((dx_2 > 0) != (dy_1 > 0))) {
          // second term also overflows below
          return dx_1 * dy_2 > dx_2 * dy_1;
        }
        else {
          // only first term overflows below
          return false;
        }
      }
    }
    else {
      // second term overflows; result is inverse of the direction
      return (dx_2 > 0) != (dy_1 > 0);
    }
    return dx_1 * dy_2 > dx_2 * dy_1;
  }
  else {
    // standard setting
    return dx_1 * dy_2 > dx_2 * dy_1;
  }
}

template<typename P>
bool no_turn(P const& p, P const& q, P const& r) {
  
#ifdef MEASURE_TURNS
  
  ++turns;
  
#endif
  
  using T = typename P::coordinate;
  using Z = long long;
  Z mid_x = Z(q.x);
  Z mid_y = Z(q.y);
  Z dx_1 = mid_x - Z(p.x);
  Z dx_2 = Z(r.x) - mid_x;
  Z dy_1 = mid_y - Z(p.y);
  Z dy_2 = Z(r.y) - mid_y;
  
  bool check1 = Z(T(dx_1)) != dx_1 and Z(T(dy_2)) != dy_2;
  bool check2 = Z(T(dx_2)) != dx_2 and Z(T(dy_1)) != dy_1;
  if (check1 or check2) {
    // if they do not both overflow, they must be different
    if (check1 and check2) {
      // if overflow direction not the same, they are different
      if (((dx_1 > 0) == (dy_2 > 0)) == ((dx_2 > 0) == (dy_1 > 0))) {
        // the non-overflow part must also be the same
        return dx_1 * dy_2 == dx_2 * dy_1;
      }
    }
    return false;
  }
  else {
    return dx_1 * dy_2 == dx_2 * dy_1;
  }
}

// may do controlled overflow on very large values

template<typename P>
auto signed_area(P const& p, P const& q, P const& r) {
  assert(not right_turn(p, q, r));
  
#ifdef MEASURE_TURNS
  
  ++turns;
  
#endif
  
  using Z = long long;
  Z mid_x = Z(q.x);
  Z mid_y = Z(q.y);
  Z dx_1 = mid_x - Z(p.x);
  Z dx_2 = Z(r.x) - mid_x;
  Z dy_1 = mid_y - Z(p.y);
  Z dy_2 = Z(r.y) - mid_y;
  return (unsigned long long) (dx_1 * dy_2 - dx_2 * dy_1);
}
}

using namespace double_precision;

#endif

#if defined(FLOATING_POINT_FILTER)

namespace floating_point_filter {

template<typename P>
bool left_turn(P const& p, P const& q, P const& r) {
  using R = double;
  constexpr R u = std::numeric_limits<R>::epsilon();
  R last_x = R(r.x);
  R last_y = R(r.y);
  R lhs = (R(p.x) - last_x) * (R(q.y) - last_y);
  R rhs = (R(p.y) - last_y) * (R(q.x) - last_x);
  R det = lhs - rhs;
  R detsum = 0.0;
  if (lhs > 0.0) {
    if (rhs <= 0.0) {
      return lhs > rhs;
    }
    else {
      detsum = lhs + rhs;
    }
  }
  else if (lhs < 0.0) {
    if (rhs >= 0.0) {
      return lhs > rhs;
    }
    else {
      detsum = -lhs - rhs;
    }
  }
  else {
    return lhs > rhs;
  }
  R errbound = (3 * u + 16 * u * u) * detsum;
  if (det >= errbound or -det >= errbound) {
    return lhs > rhs;
  }
  return multi_precision::left_turn(p, q, r);
}

template<typename P>
bool no_turn(P const& p, P const& q, P const& r) {
  using R = double;
  constexpr R u = std::numeric_limits<R>::epsilon();
  R last_x = R(r.x);
  R last_y = R(r.y);
  R lhs = (R(p.x) - last_x) * (R(q.y) - last_y);
  R rhs = (R(p.y) - last_y) * (R(q.x) - last_x);
  R det = lhs - rhs;
  R detsum = 0.0;
  if (lhs > 0.0) {
    if (rhs <= 0.0) {
      return false;
    }
    else {
      detsum = lhs + rhs;
    }
  }
  else if (lhs < 0.0) {
    if (rhs >= 0.0) {
      return false;
    }
    else {
      detsum = -lhs - rhs;
    }
  }
  else {
    return lhs == rhs;
  }
  R errbound = (3 * u + 16 * u * u) * detsum;
  if (det >= errbound or -det >= errbound) {
    return lhs == rhs;
  }
  return multi_precision::no_turn(p, q, r);
}
}

using namespace floating_point_filter;

#endif

template<typename P>
bool right_turn(P const& p, P const& q, P const& r) {
  return left_turn(r, q, p);
}

// is p on the line segment {q, r}?

template<typename P>
bool on_line_segment(P const& p, P const& q, P const& r) {
  P left = std::min(q, r);
  P right = std::max(q, r);
  if (p < left) {
    return false;
  }
  if (p > right) {
    return false;
  }
  return no_turn(p, q, r);
}

#endif