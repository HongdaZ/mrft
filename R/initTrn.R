## Initialize data for training
initTrn <- function( img, label_, modality ) {
  
  switch( modality,
          "flair" = {
            # label[ label == 1 | label == 4 ] <- NA_integer_
            label <- changeD( label_, 1L, 4L )
          },
          "t1ce" = {
            # label[ label == 1 | label == 2 ] <- NA_integer_
            label <- changeD( label_, 1L, 2L )
          },
          "t2" = {
            # label[ label == 4 ] <- NA_integer_
            label <- changeD( label_, 4L, NA_integer_ )
          } )
  info <- indexMat( img, label )
  seg <- label[ info$idx ]
  seg <- rbind( seg, rep( 0L, length( seg ) ) )

  list( info = info, seg = seg )
  
}