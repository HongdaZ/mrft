## Postprocess the results
postProcess <- function( post_data, min_enh, min_enh_enc,
                         max_prop_enh_enc, max_prop_enh_slice,
                         min_tumor, spread_add, spread_rm, 
                         trim1_spread, trim1_round, 
                         remove2d_spread, remove2d_round, 
                         spread_trim, round_trim,
                         on_flair_prop ) {
  .Call( "postProcess", post_data, min_enh, min_enh_enc,
         max_prop_enh_enc, max_prop_enh_slice,
         min_tumor, spread_add, spread_rm, 
         trim1_spread, trim1_round, 
         remove2d_spread, remove2d_round, 
         spread_trim, round_trim,
         on_flair_prop )
}