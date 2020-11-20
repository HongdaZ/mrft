## Postprocess the results
postProcess <- function( post_data, min_enh, min_enh_enc,
                         max_prop_enh, min_tumor, spread_add,
                         spread_rm ) {
  .Call( "postProcess", post_data, min_enh, min_enh_enc,
          max_prop_enh, min_tumor, spread_add, spread_rm )
}