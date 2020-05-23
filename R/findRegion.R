# Return voxel index vector, neighbor index matrix, intensity vector
# and neighbor intensity matrix
findRegion <- function(  label, nidx, start ) {
  .Call( "findRegion", label, nidx, start )
}