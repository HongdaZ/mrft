## Initialize data for training
initTrn <- function( file ) {
  img <- readIntst( file )
  if( grepl( "flair.nii", file ) ) {
    label <- gsub( "flair", "seg",
        file )
  } else if( grepl( "t1.nii", file ) ) {
    label <- gsub( "t1", "seg",
                 file )
  } else if( grepl( "t1ce.nii", file ) ) {
    label <- gsub( "t1ce", "seg",
                 file )
  } else if( grepl( "t2.nii", file ) ) {
    label <- gsub( "t2", "seg",
                 file )
  }
  lbl <- readNifti( label )
}