# Split CSF & necrosis, white matter and grey matter in training data
# csf = 1, wm = 2, g = 3, TBD = 0, o.w. = NA 
splitCWG <- function( patient, modality ) {
  
  lbl <- readNifti( patient[ 2 ] )
  
  flair <- readNifti( patient[ 1 ] )
  t2 <- readNifti( patient[ 5 ] )
  t1ce <- readNifti( patient[ 4 ] )
  
  switch( modality,
          "flair" = {
            img <- flair
          },
          "t1ce" = {
            img <- t1ce
          },
          "t2" = {
            img <- t2
          } )
  
  flair[ lbl != 0 ] <- NA
  t2[ lbl != 0 ] <- NA
  t1ce[ lbl != 0 ] <- NA
  
  
  
  # Find csf & necrosis
  q_flair <- quantile( flair, probs = c( .20, .40, .50, .01 ) , na.rm = T )
  q_t2 <- quantile( t2, probs = c( .80, .01 ), na.rm = T )
  
  csf <- flair < q_flair[ 1 ] & t2 > q_t2[ 1 ]
  lbl[ csf ] <- -1
  
  # Find white matter
  q_t1ce <- quantile( t1ce, probs = c( .60, .50, .01 ), na.rm = T )
  
  wm <- flair < q_flair[ 2 ] & t1ce > q_t1ce[ 1 ]
  lbl[ wm ] <- -2
  
  # Find grey matter
  gm <- flair > q_flair[ 3 ] & t1ce < q_t1ce[ 2 ]
  lbl[ gm ] <- -3
  # Remove the darkest 1% of the voxels
  switch( modality,
          "flair" = {
            lbl[ flair < q_flair[ 4 ] ] <- NA
          },
          "t1ce" = {
            lbl[ t1ce < q_t1ce[ 3 ] ] <- NA
          },
          "t2" = {
            lbl[ t2 < q_t2[ 2 ] ] <- NA
          } )
  list( label = lbl, intensity = img )
}