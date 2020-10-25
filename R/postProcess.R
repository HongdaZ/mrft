## Postprocess the results
postProcess <- function( post_data, min_enh, min_tumor, min_prop_net ) {
  .Call( "postProcess", post_data, min_enh, min_tumor, min_prop_net )
}