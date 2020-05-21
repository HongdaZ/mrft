## Read intensity values
# Read a NIfTI file and return a intensity vector, 
# a neighboring voxel intensity matrix, a index vector and 
# a neighbor index matrix 

readIntst <- function( file ) {
  n_nbr <- 6
  img <- readNifti( file ) # colum major array
  indexMat( img )
} 