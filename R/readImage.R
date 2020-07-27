# read the zipped images
readImage <- function( patient ) {
  label <- readNIfTI( patient[ 2 ], reorient = FALSE )@.Data
  flair <- readNIfTI( patient[ 1 ], reorient = FALSE )@.Data
  t2 <- readNIfTI( patient[ 5 ], reorient = FALSE )@.Data
  t1ce <- readNIfTI( patient[ 4 ], reorient = FALSE )@.Data
  res <- list( label = label, t1ce = t1ce, flair = flair, t2 = t2 )
  return( res )
}