## Initialize data for training
initTrn <- function( file ) {
  
  if( grepl( "flair.nii", file ) ) {
    
    label <- gsub( "flair", "seg", file )
    lbl <- readNifti( label )
    lbl[ lbl == 1 | lbl == 4 ] <- NA
    # No use of t1 images here
  } else if( grepl( "t1ce.nii", file ) ) {
    
    label <- gsub( "t1ce", "seg", file )
    lbl <- readNifti( label )
    lbl[ lbl == 1 | lbl == 2 ] <- NA
    
  } else if( grepl( "t2.nii", file ) ) {
    
    label <- gsub( "t2", "seg", file )
    lbl <- readNifti( label )
    lbl[ lbl == 4 ] <- NA
  }
  img <- readIntst( file, lbl )
  lbl <- lbl[ img$idx ]
}