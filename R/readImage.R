# read the zipped images
readImage <- function( patient ) {
  
  t1ce <- readNIfTI( patient[ 3 ], reorient = FALSE )@.Data
  flair <- readNIfTI( patient[ 1 ], reorient = FALSE )@.Data
  t2 <- readNIfTI( patient[ 4 ], reorient = FALSE )@.Data
  
  t1ce[ t1ce == 0 ] <- NaN
  flair[ flair == 0 ] <- NaN
  t2[ t2 == 0 ] <- NaN
  
  res <- list( t1ce = t1ce, flair = flair, t2 = t2 )
  return( res )
}