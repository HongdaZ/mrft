## Find threshold for tumor
threshold <- function( res, delta ) {
  nc <- ncol( res )
  m <- res[ 1, nc ]
  sigma2 <- res[ 2, nc ]
  shift <- delta[ 3 ] - delta[ nc - 1 ]
  m + sqrt( sigma2 ) * shift
} 