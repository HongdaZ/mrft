# m1, m2, m3, m4
priorMode <- function( model ) {
  intst <- model$info$intst
  seg <- model$seg
  .Call( "priorMode", intst, seg )
}