## Initialize data for training
initTrn <- function( img, label, modality ) {
  
  switch( modality,
          "flair" = {
            label[ label == 1 | label == 4 ] <- NA_integer_
          },
          "t1ce" = {
            label[ label == 1 | label == 2 ] <- NA_integer_
          },
          "t2" = {
            label[ label == 4 ] <- NA_integer_
          } )
  info <- indexMat( img, label )
  seg <- label[ info$idx ]
  seg <- rbind( seg, rep( 0L, length( seg ) ) )

  list( info = info, seg = seg )
  
  # regions <- NULL
  # # seg is numeric
  # n_region <- 0
  # for( i in 1 : dim( seg )[ 2 ] ) {
  #   if( seg[ 1, i ]  > 0 ) {
  #     region <- findRegion( seg, info$nidx, i )
  #     seg[ 1, region ] <- - 4 - n_region
  # 
  #     regions[[ n_region + 1 ]] <- list( label = - 4 - n_region,
  #                                         region = region )
  #     n_region <- n_region + 1
  #   }
  # }
  
}