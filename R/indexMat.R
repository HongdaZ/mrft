# Return voxel index vector, neighbor index matrix, intensity vector
# and neighbor intensity matrix
indexMat <- function(  img, p = 12 ) {
  .Call( "indexMat", img, p )
}