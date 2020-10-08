#include	<stdio.h>
#include	<math.h>
#define	ERROR	1
#define	X	0
#define	Y	1

#define DIM     2               /* Dimension of points */
typedef int     tPointi[DIM];   /* type integer point */
typedef double  tPointd[DIM];   /* type double point */
#define PMAX    1000            /* Max # of pts in polygon */

typedef tPointi tPolygoni[PMAX];/* type integer polygon */

void	PrintPoly( int n, tPolygoni P );
void	PrintPoint( tPointi p );
bool	InPoly( tPointi q, tPolygoni P, int n );
vector<int> inpoly( const vector<int> &p, const vector<int> &poly ) {
  int len = poly.size() / 2;
  int n_point = p.size() / 2;
  vector<Point> polygon( len );
  vector<int> inside( n_point );
  Point point;
  for( int i = 0; i < len; ++ i ) {
    polygon[ i ].x = poly[ 2 * i ];
    polygon[ i ].y = poly[ 2 * i + 1 ];
  }
  for( int i = 0; i < n_point; ++ i ) {
    point.x = p[ 2 * i ];
    point.y = p[ 2 * i + 1 ];
    inside[ i ] = isInside( polygon, len, point );
  }
  return inside;
}