# Return voxel index vector, neighbor index matrix, intensity vector
# and neighbor intensity matrix
indexMat <- function(  img ) {
  .Call( "indexMat", img )
}