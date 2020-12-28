## Postprocess the results
postProcess <- function( post_data, min_enh, min_enh_enc,
                         max_prop_enh_enc, max_prop_enh_slice,
                         min_tumor, spread_add,
                         spread_rm, spread_trim, round_trim ) {
  .Call( "postProcess", post_data, min_enh, min_enh_enc,
         max_prop_enh_enc, max_prop_enh_slice,
         min_tumor, spread_add, spread_rm, spread_trim, round_trim )
}