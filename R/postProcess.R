## Postprocess the results
postProcess <- function( post_data, min_enh, 
                         max_prop_enh, min_tumor, spread,
                         min_prop_tumor_nbr ) {
  .Call( "postProcess", post_data, min_enh,
          max_prop_enh, min_tumor, spread,
          min_prop_tumor_nbr )
}