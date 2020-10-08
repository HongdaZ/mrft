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

/*
 Returns true if q is inside polygon P.
 */
bool	InPoly( tPointi q, tPolygoni P, int n )
{
  int	i, i1;		/* point index; i1 = i-1 mod n */
int	d;		/* dimension index */
double	x;		/* x intersection of e with ray */
int	crossings = 0;	/* number of edge/ray crossings */

printf("\n==>In: u = "); PrintPoint(q); putchar('\n');

/* Shift so that q is the origin. */
for( i = 0; i < n; i++ ) {
  for( d = 0; d < DIM; d++ )
    P[i][d] = P[i][d] - q[d];
}

/* For each edge e=(i-1,i), see if crosses ray. */
for( i = 0; i < n; i++ ) {
  i1 = ( i + n - 1 ) % n;
  printf("e=(%d,%d)\t", i1, i);
  /* if e straddles the x-axis... */
  if( ( ( P[i] [Y] > 0 ) && ( P[i1][Y] <= 0 ) ) ||
  ( ( P[i1][Y] > 0 ) && ( P[i] [Y] <= 0 ) ) ) {
    /* e straddles ray, so compute intersection with ray. */
    x = (P[i][X] * P[i1][Y] - P[i1][X] * P[i][Y])
    / (double)(P[i1][Y] - P[i][Y]);
    printf("straddles: x = %g\t", x);
    /* crosses ray if strictly positive intersection. */
    if (x > 0) crossings++;
  }
  printf("crossings=%d\n", crossings);
}
/* q inside if an odd number of crossings. */
if( (crossings % 2) == 1 )
  return	true;
else	return	false;
}
void	PrintPoint( tPointi p )
{
  int	i;
  
  putchar('(');
  for ( i = 0; i < DIM; i++ ) {
    printf("%d", p[i]);
    if ( i != DIM-1 ) putchar(',');
  }
  putchar(')');
}
/*
 Reads in the coordinates of the vertices of a polygon from stdin,
 puts them into P, and returns n, the number of vertices.
 Formatting conventions: etc.
 */
int	ReadPoly( tPolygoni P )
{
  int	n = 0;
  
  printf("Polygon:\n");
  printf("  i   x   y\n");
  while ( (n < PMAX) && (scanf("%d %d",&P[n][0],&P[n][1]) != EOF) ) {
    printf("%3d%4d%4d\n", n, P[n][0], P[n][1]);
    ++n;
  }
  if (n < PMAX)
    printf("n = %3d vertices read\n",n);
  else	printf("Error in read_poly:  too many points; max is %d\n", PMAX);
  putchar('\n');
  
  return	n;
}

void	PrintPoly( int n, tPolygoni P )
{
  int	i;
  
  printf("Polygon:\n");
  printf("  i   x   y\n");
  for( i = 0; i < n; i++ )
    printf("%3d%4d%4d\n", i, P[i][0], P[i][1]);
}
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