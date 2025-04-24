# connectUp ------------
#' Assign connected regions in a window
#'
#' This function assigns random connected regions on a square. Used within function `RandomRegionWindow`.
#' @param data nRegion is the No. of regions (e.g. nRegion=3)
#' @param method poly is a RasterLayer  (e.g. 20 by 20 square).
#' @return A list of the selected polygons for each region.

EstCCI = function(expr_data, spatial_data,
                  method=c("loc_loc", "expr_loc", "expr_expr")) {
  # load data

  # run a simple CCI analysis for co-localization
  if (method=="loc_loc") {


  }
  # run a simple CCI analysis for expression vs neighborhood cells

  # run a simple CCI analysis for L-R pairs in the neighborhood

  # save results in the  format that can be directly used in sCCIgen

}
load("Github/sCCIgen_data/InputData/expression_data/SeqFishPlusCortex_033023_expr.Rdata")
load("Github/sCCIgen_data/InputData/cell_feature_data/SeqFishPlusCortex_033023_cellfeature.Rdata")

expr_data=expr
spatial_data=cell_feature





