## Initialize data for training
initTrn <- function( patient, modality ) {
  
  switch( modality,
          "flair" = {
            lbl <- splitCWG( patient )
            lbl[ lbl == 1 | lbl == 4 ] <- NA
          },
          "t1ce" = {
            lbl <- splitCWG( patient )
            lbl[ lbl == 1 | lbl == 2 ] <- NA
          },
          "t2" = {
            lbl <- splitCWG( patient )
            lbl[ lbl == 4 ] <- NA
          } )
  img <- readIntst( file, lbl )
  lbl <- lbl[ img$idx ]
  # Find initial regions of csf, wm and gm
}