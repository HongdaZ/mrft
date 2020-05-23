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
  label <- label[ img$idx ]
  label <- rbind( label, rep( 0, length( label ) ) )
  # Find initial regions of csf, wm and gm
}