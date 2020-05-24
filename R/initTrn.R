## Initialize data for training
initTrn <- function( patient, label, modality ) {
  
  switch( modality,
          "flair" = {
            file <- patient[ 1 ]
            label[ label == 1 | label == 4 ] <- NA
          },
          "t1ce" = {
            file <- patient[ 4 ]
            label[ label == 1 | label == 2 ] <- NA
          },
          "t2" = {
            file <- patient[ 5 ]
            label[ label == 4 ] <- NA
          } )
  img <- readIntst( file, label )
  seg <- label[ img$idx ]
  seg <- rbind( seg, rep( 0, length( seg ) ) )

  regions <- NULL
  # seg is numeric
  n_region <- 0
  for( i in 1 : dim( seg )[ 2 ] ) {
    if( seg[ 1, i ]  > 0 ) {
      region <- findRegion( seg, img$nidx, i )
      seg[ 1, region ] <- - 4 - n_region
      
      regions[[ n_region + 1 ]] <- list( label = - 4 - n_region,
                                          region = region )
      n_region <- n_region + 1
    }
  }
  

}