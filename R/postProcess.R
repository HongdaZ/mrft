## Postprocess the results
postProcess <- function( post_data ) {
  .Call( "postProcess", post_data )
}