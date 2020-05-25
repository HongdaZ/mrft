## Initialize data for training
initTrn <- function( img, label, modality ) {
  
  switch( modality,
          "flair" = {
            label[ label == 1 | label == 4 ] <- NA
          },
          "t1ce" = {
            label[ label == 1 | label == 2 ] <- NA
          },
          "t2" = {
            label[ label == 4 ] <- NA
          } )
  info <- readIntst( img, label )
  seg <- label[ info$idx ]
  seg <- rbind( seg, rep( 0, length( seg ) ) )

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