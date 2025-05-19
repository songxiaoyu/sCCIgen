# EstCCI ------------
#' Estimate three patterns for cell-cell interactions
#'
#' @param expr_data Input expression data.
#' @param spatial_data Input spatial data. 
#' @param paired If the expression and spatial data are paired. If not, analysis will not be run.
#' @param LR_database sCCIgen provides two ligand-recepter databases from CellTalkDB for mouse and human. Users
#' can select "external" and provide the path to external database at LR_database_path. Users can also select
#' "null", which will not run LR analysis.)
#' @param save_folder Provide a path to a folder that saves the CCI results for use in sCCIgen.
#' @return No return inside R. Results are saved in `save_folder`. 
#' File est_CCI_dist_dist.tsv will save results of CCI through cell location attraction and inhibition.
#' File est_CCI_expr_dist.tsv will save results of CCI through expression via cell location impact.
#' File est_CCI_expr_expr.tsv will save results of CCI through ligand and recepter expression.
#' @export

EstCCI = function(expr_data, spatial_data, paired=T, LR_database=c("mouse", "human", "external", "null"), 
                  LR_database_path=NULL, save_folder=getwd()) {

  if (paired==F) {stop("The expression and spatial data are not paired and analysis of CCI can not be done.")}
  
  
  # process data with Giotto
  anno=colnames(expr_data)
  colnames(expr_data) = paste0("cell_", 1:ncol(expr_data))
  expr_data=as.data.frame(expr_data)
  loc=spatial_data[,2:3]
  moreMeta=as.matrix(spatial_data[,-c(2,3)])
  colnames(moreMeta)[1]="anno"
  
  dat=preprocessGiotto(expr_data=expr_data, loc=loc, moreMeta=moreMeta)
  # run a simple CCI analysis for co-localization
  
  CCI1_Table=cellProximityTable(gobject=dat, abs_enrichm=0.3, p_adj = 0.05) 
  # run a simple CCI analysis for expression vs neighborhood 
 
  CCI2_Table=ExprDistanceTable(gobject=dat, abs_log2fc_ICG=0.25, p_adj = 0.05)
  # run a simple CCI analysis for L-R pairs in the neighborhood
  CCI3_Table=ExprExprTable(gobject=dat, LR_database=LR_database, LR_database_path=LR_database_path, 
                           abs_log2fc_LR=0.25, p_adj=0.05)
  # save results in the  format that can be directly used in sCCIgen
  readr::write_tsv(CCI1_Table, file=file.path(save_folder, "est_CCI_dist_dist.tsv"))
  
  readr::write_tsv(CCI2_Table, file=file.path(save_folder, "est_CCI_expr_dist.tsv"))
  
  if (is.null(CCI3_Table)==F){
    readr::write_tsv(CCI3_Table, file=file.path(save_folder, "est_CCI_expr_expr.tsv"))
  } 
  print(paste("Results are saved at", save_folder))

}

# library(dplyr)
# 
# setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_sCCIgen")
# load("Github/sCCIgen_data/InputData/expression_data/SeqFishPlusCortex_033023_expr.Rdata")
# load("Github/sCCIgen_data/InputData/cell_feature_data/SeqFishPlusCortex_033023_cellfeature.Rdata")
# 
# # This cleaning step will be added in data processing - colname is cell type; data type for expr is matrix. 
# expr_data=expr %>% as.matrix()
# spatial_data=cell_feature
# colnames(expr_data)=spatial_data$anno
# 
# 
# # Use Giotto - expr data frame; loc data; moreMeta data
# anno=colnames(expr_data)
# colnames(expr_data) = paste0("cell_", 1:ncol(expr_data))
# expr_data=as.data.frame(expr_data)
# loc=spatial_data[,2:3]
# moreMeta=as.matrix(spatial_data[,-c(2,3)])
# colnames(moreMeta)[1]="anno"


##%######################################################%##
#                                                          #
####            0) Preprocess in Giotto              ####
#                                                          #
##%######################################################%##
# preprocessGiotto ------------
#' Process SRT data using Giotto for CCI analysis. 
#'
#' This function takes input data and goes through data processing procedures for CCI analysis.

preprocessGiotto=function(expr_data, loc, moreMeta) {
  # create Giotto object
  dat <- createGiottoObject(expression = expr_data, spatial_locs = loc)
  
  # add additional annotation if wanted
  dat <- addCellMetadata(dat,new_metadata = moreMeta)
  
  # filter
  dat <- filterGiotto(gobject = dat,
                      expression_threshold = 1,
                      feat_det_in_min_cells = 10,
                      min_det_feats_per_cell = 10,
                      expression_values = "raw",
                      verbose = TRUE)
  
  # normalize
  dat <- normalizeGiotto(gobject = dat, 
                         scalefactor = 6000, 
                         verbose = TRUE)
  
  # add gene & cell statistics
  dat <- addStatistics(gobject = dat)
  
  # adjust expression matrix for technical or known variables
  dat <- adjustGiottoMatrix(gobject = dat, 
                            expression_values = "normalized",
                            covariate_columns = c("nr_feats", "total_expr"),
                            return_gobject = TRUE,
                            name = "custom")
  
  ## highly variable features (HVF)
  dat <- calculateHVF(gobject = dat)
  
  # networks
  dat <- createSpatialNetwork(gobject = dat,minimum_k = 2,maximum_distance_delaunay = 400)
  dat <- createSpatialNetwork(gobject = dat, method = "kNN", k = 5, name = "spatial_network")
  
  return(dat)
}



##%######################################################%##
#                                                          #
####            1) Cell-proximity analysis              ####
#                                                          #
##%######################################################%##

cellProximityTable = function (gobject, abs_enrichm=0.3, p_adj = 0.05) {
  
  # cellProximityEnrichment analysis
  
  cell_proximities <- cellProximityEnrichment(gobject = dat,
                                              cluster_column = "anno",
                                              spatial_network_name = "Delaunay_network",
                                              adjust_method = "fdr",
                                              number_of_simulations = 2000)
  # create table
  table_mean_results_dc <- CPscore$enrichm_res
  table_mean_results_dc_filter <- table_mean_results_dc[ abs(enrichm) >= abs_enrichm, ]
  table_mean_results_dc_filter <- table_mean_results_dc_filter[p.adj_higher <= p_adj | p.adj_lower <= p_adj, ]

  table_separated <- table_mean_results_dc_filter[, c(1,5)] %>%
    mutate(unified_int = as.character(unified_int)) %>%
    tidyr::separate(unified_int, into = c("CellTye1", "CellTye2"), sep = "--")
  
  return(table_separated)
}




##%######################################################%##
#                                                          #
#### 2) Cell Neighborhood: Interaction Changed Features ####
#                                                          #
##%######################################################%##
ExprDistanceTable = function (gobject, abs_log2fc_ICG=0.25, p_adj = 0.05) {
  plan(future::multisession)  
  # ## select genes based on highly variable features and gene statistics, both found in feature (gene) metadata
  gene_metadata <- fDataDT(dat)
  heg <- gene_metadata[hvf == "yes" & perc_cells > 4 & mean_expr_det > 0.5]$feat_ID
  
  ## identify genes that are associated with proximity to other cell types
  ICFsForesHighGenes <- findInteractionChangedFeats(
    gobject = dat,
    selected_feats = heg,
    spatial_network_name = "Delaunay_network",
    cluster_column = "anno",
    diff_test = "permutation",
    adjust_method = "fdr",
    nr_permutations = 2000,
    do_parallel = TRUE
  )
  
  
  df <- as_tibble(ICFsForesHighGenes$ICFscores)
  
  df_table <- df %>% filter(p.adj <= p_adj, abs(log2fc) >= abs_log2fc_ICG) 
  
  df_table2=df_table[,c("cell_type", "int_cell_type", "feats", "log2fc")]
  
  return(df_table2)

}



##%######################################################%##
#                                                          #
####            3) Ligand-Receptor analysis             ####
#                                                          #
##%######################################################%##
ExprExprTable = function(gobject, LR_database=c("mouse", "human", "external", "null"), LR_database_path=NULL,
                         p_adj=0.05, abs_log2fc_LR=0.25) {
  plan(future::multisession)  
  # load database
  
  if (LR_database=="null") {
    message("Ligand-Receptor analysis is not performed. Skipping to the next step.")
  } else{
    if (LR_database=="mouse") {
      LR_data <- readRDS(system.file("mouse_lr_pair.rds", package = "sCCIgen"))
    }
    if (LR_database=="human") {
      LR_data <- readRDS(system.file("human_lr_pair.rds", package = "sCCIgen"))
    }
    if (LR_database=="external") {
      LR_data <- data.table::fread(LR_database_path)
      message("Please ensure variables ligand_gene_symbol and receptor_gene_symbol exist in your external LR database.")
    }
    
    ## LR expression
    LR_data[, ligand_det := ifelse(LR_data$ligand_gene_symbol %in% dat@feat_ID$rna, TRUE, FALSE)]
    LR_data[, receptor_det := ifelse(LR_data$receptor_gene_symbol %in% dat@feat_ID$rna, TRUE, FALSE)]
    LR_data_det <- LR_data[ligand_det == TRUE & receptor_det == TRUE]
    
    select_ligands <- LR_data_det$ligand_gene_symbol
    select_receptors <- LR_data_det$receptor_gene_symbol
    
    ## get gene pair expression changes based on expression - only score
    expr_only_scores <- exprCellCellcom(gobject = dat,
                                        cluster_column = "anno",
                                        random_iter = 500,
                                        feat_set_1 = select_ligands,
                                        feat_set_2 = select_receptors,
                                        verbose = FALSE
    )
    
    ## get statistical significance of gene pair expression changes upon cell-cell interaction
    spatial_all_scores <- spatCellCellcom(dat, spatial_network_name = "spatial_network", 
                                          cluster_column = "anno",
                                          random_iter = 500, feat_set_1 = select_ligands, 
                                          feat_set_2 = select_receptors,
                                          adjust_method = "fdr", 
                                          do_parallel = TRUE, cores = 6, 
                                          verbose = "a little")
    
    ## select top LR table ##
    selected_spat <- spatial_all_scores[p.adj <= p_adj & abs(log2fc) > abs_log2fc_LR & lig_nr >= 5 & rec_nr >= 5]
    tb= selected_spat[, c("lig_cell_type", "rec_cell_type", "ligand", "receptor", "log2fc")]
    
    return(tb)
    
    
  }
  
 
} 