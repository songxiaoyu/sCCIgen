##%######################################################%##
#                                                          #
####            0) Preprocess in Giotto              ####
#                                                          #
##%######################################################%##

#' Process SRT data using Giotto for CCI analysis.
#'
#' This function takes input data and goes through data processing procedures for CCI analysis.
#'
#' #' @import Giotto
#' @export

preprocessGiotto=function(expr_data, spatial_data, run_hvg=T, run_kNN_network=T, run_Delaunay_network=T) {
  # harmonize input with Giotto requirement
  anno=colnames(expr_data)
  colnames(expr_data) = paste0("cell_", 1:ncol(expr_data))
  expr_data=as.data.frame(expr_data)
  loc=spatial_data[,2:3]
  moreMeta=as.matrix(spatial_data[,-c(2,3)])
  colnames(moreMeta)[1]="anno"

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
                      verbose = F)

  # normalize
  dat <- normalizeGiotto(gobject = dat,
                         scalefactor = 6000,
                         verbose = F)

  # add gene & cell statistics
  dat <- addStatistics(gobject = dat)

  # adjust expression matrix for technical or known variables
  dat <- adjustGiottoMatrix(gobject = dat,
                            expression_values = "normalized",
                            covariate_columns = c("nr_feats", "total_expr"),
                            return_gobject = T,
                            name = "custom")

  ## highly variable features (HVF)
  if(run_hvg==T) {dat <- calculateHVF(gobject = dat)}

  # networks
  if (run_Delaunay_network==T) {dat <- createSpatialNetwork(gobject = dat,minimum_k = 2,maximum_distance_delaunay = 400)}
  if (run_kNN_network==T) {dat <- createSpatialNetwork(gobject = dat, method = "kNN", k = 5, name = "spatial_network")}

  return(dat)
}



##%######################################################%##
#                                                          #
####            1) Cell-proximity analysis              ####
#                                                          #
##%######################################################%##
# EstCCI ------------
#' Estimate cell attraction and inhibition patterns and save the table for simulation
#'
#' @param expr_data Input expression data.
#' @param spatial_data Input spatial data.
#' @param output_file Provide a path and file name to save the CCI results for use in sCCIgen.
#' @param abs_enrichm Effect size threshold.
#' @param p_adj Adjusted p-value threshold (Default = 0.05).

#' File est_CCI_dist_dist.tsv will save results of CCI through cell location attraction and inhibition.
#' File est_CCI_expr_dist.tsv will save results of CCI through expression via cell location impact.
#' File est_CCI_expr_expr.tsv will save results of CCI through ligand and recepter expression.
#' @export

cellProximityTable = function (gobject, output_file=file.path(getwd(), "est_CCI_dist_dist.csv"), abs_enrichm=0.3, p_adj = 0.05) {

  # cellProximityEnrichment analysis
  cell_proximities <- cellProximityEnrichment(gobject = gobject,
                                              cluster_column = "anno",
                                              spatial_network_name = "Delaunay_network",
                                              adjust_method = "fdr",
                                              number_of_simulations = 2000)
  # create table
  dc <- cell_proximities$enrichm_res
  dc2 <- dc[abs(dc$enrichm) >= abs_enrichm, ]
  dc2 <- dc2[dc2$p.adj_higher <= p_adj | dc2$p.adj_lower <= p_adj, ]



  # format into sCCIgen input
  # <cell_type_A>,<cell_type_B>,<value>. Separate the entries by blank space (e.g. 'cell_type_A,cell_type_B,1.2 cell_type_B,cell_type_C,-0.8').
  # Alternatively, select a file separated by commas, where each line is an entry with the format described above.
  CCI1_Table <- dc2[, c(1,5)] %>%
    dplyr::mutate(unified_int = as.character(unified_int)) %>%
    dplyr::mutate(enrichm=round(enrichm, 2)) %>%
    tidyr::separate(unified_int, into = c("CellTye1", "CellTye2"), sep = "--")

  # save results in the  format that can be directly used in sCCIgen
  readr::write_csv(CCI1_Table, file=output_file, col_names = F)

  print(paste0("CCI co-localization analysis is saved at: ", output_file))

}




##%######################################################%##
#                                                          #
#### 2) Cell Neighborhood: Interaction Changed Features ####
#                                                          #
##%######################################################%##

#' Estimate cell expression - distance patterns and save the table for simulation
#'
#' @param expr_data Input expression data.
#' @import future
#' @export


ExprDistanceTable = function(gobject, in_hvg=F, output_file=file.path(getwd(), "est_CCI_dist_expr.csv"),
                              region_specific=F, abs_log2fc_ICG=0.25, p_adj = 0.05) {
  plan(future::multisession)

  # if in_hvg ==T, select genes based on highly variable features and gene statistics, both found in feature (gene) metadata
  if (in_hvg==T) {
    gene_meta <- fDataDT(gobject)
    heg <- gene_meta[which(gene_meta$hvf == "yes" & gene_meta$perc_cells > 4 & gene_meta$mean_expr_det > 0.5), ]$feat_ID
    gobject <- subsetGiotto(gobject = gobject, feat_ids = heg)
  }
  cell_meta =pDataDT(gobject)
  # Regions specific analysis
  if (region_specific==T) {

    R=unique(cell_meta[,3]) %>% as.matrix()

    res=NULL
    for (r in R) {
      print(r)
      cell_id_r=cell_meta[which(cell_meta[,3] == r), 1] %>% as.matrix()
      dat_r <- subsetGiotto(gobject = gobject, cell_ids = cell_id_r )
      res1=ExprDistanceTable_1region(gobject=dat_r, r=r, cell_meta=cell_meta,
                                     abs_log2fc_ICG=abs_log2fc_ICG, p_adj = p_adj)
      res=rbind(res, res1)
    }

  }
  if (region_specific==F) {

    res=ExprDistanceTable_1region(gobject=gobject, r="NULL", cell_meta=cell_meta,
                                   abs_log2fc_ICG=abs_log2fc_ICG, p_adj = p_adj)
  }
  # save results in the  format that can be directly used in sCCIgen
  readr::write_csv(res, file=output_file, col_names = F)

  print(paste0("CCI expression-distance association is saved at: ", output_file))

}


ExprDistanceTable_1region = function (gobject, r, cell_meta, abs_log2fc_ICG=0.25, p_adj = 0.05) {

  ## identify genes that are associated with proximity to other cell types
  ICFsForesHighGenes <- findInteractionChangedFeats(
    gobject = gobject,
    spatial_network_name = "Delaunay_network",
    cluster_column = "anno",
    diff_test = "permutation",
    adjust_method = "fdr",
    nr_permutations = 1000,
    do_parallel = TRUE
  )
  df <- dplyr::as_tibble(ICFsForesHighGenes$ICFscores)
  df_table <- df %>% dplyr::filter(p.adj <= p_adj, abs(log2fc) >= abs_log2fc_ICG)
  if (nrow(df_table) > 0) {
    df_table=df_table[,c("cell_type", "int_cell_type", "feats", "log2fc")]

    # add an approximate distance threshold by cell tyep pair
    delaunay_net <- getSpatialNetwork(gobject = gobject, name = 'Delaunay_network')
    dis_table=delaunay_net@networkDT
    dis_table2 <- merge(dis_table, cell_meta, by.x = "from", by.y = "cell_ID", all.x = TRUE)
    dis_table2 <- merge(dis_table2, cell_meta, by.x = "to", by.y = "cell_ID", all.x = TRUE)
    tb=tapply(dis_table2$distance, list(dis_table2$anno.x, dis_table2$anno.y), median)
    dist_value=reshape2::melt(tb) %>%
      dplyr::filter(is.na(value)==F)%>%
      dplyr::rename(cell_type = Var1, int_cell_type = Var2, threshold = value)
    df_merged <- df_table %>% dplyr::inner_join(dist_value, by = c("cell_type", "int_cell_type"))
    # <Region>,<Perturbed cell type>,<Adjacent cell type>,<Interaction distance threshold (default 0.1)>,
    # <Gene ID (optional)>,<Gene proportion (optional)>,<Mean effect at log(count) scale (default = 0.5)>,
    # <SD of effect at log(count) scale (default = 0)>
    print(dim(df_merged))
    res=data.frame(r, df_merged[,c("cell_type", "int_cell_type", "threshold","feats")], "NULL", df_merged$log2fc, 0)

    return(res)
  } else {return(NULL)}

}



##%######################################################%##
#                                                          #
####      3) Expression-Expression (e.g. LR) analysis   ####
#                                                          #
##%######################################################%##
#' Estimate cell expression -expression patterns and save the table for simulation
#'
#' @param expr_data Input expression data.
#' @import future
#' @export
#'
ExprExprTable = function(gobject, database=c("mouse", "human", "external"), external_database_path=NULL, region_specific=F,
                         p_adj=0.05, abs_log2fc_LR=0.25, output_file=file.path(getwd(), "est_CCI_expr_expr.csv")) {
  plan(future::multisession)
  # load database

  if (database=="external" & is.null(external_database_path)) {
    message("Ligand-Receptor analysis is not performed due to a lack of LR database. ")
  } else{
    if (database=="mouse") {
      LR_data <- readRDS(system.file("mouse_lr_pair.rds", package = "sCCIgen"))
    }
    if (database=="human") {
      LR_data <- readRDS(system.file("human_lr_pair.rds", package = "sCCIgen"))
    }
    if (database=="external") {
      LR_data <- data.table::fread(external_database_path)
      message("Please ensure variables ligand_gene_symbol and receptor_gene_symbol exist in your external LR database.")
    }

    ## LR expression
    LR_data$ligand_det <- LR_data$ligand_gene_symbol %in% gobject@feat_ID$rna
    LR_data$receptor_det <- LR_data$receptor_gene_symbol %in% gobject@feat_ID$rna
    LR_data_det <- LR_data[LR_data$ligand_det & LR_data$receptor_det, ]

    select_ligands <- LR_data_det$ligand_gene_symbol
    select_receptors <- LR_data_det$receptor_gene_symbol

    cell_meta =pDataDT(gobject)
    # Regions specific analysis
    if (region_specific==T) {
      R=unique(cell_meta[,3]) %>% as.matrix()
      res=NULL
      for (r in R) {
        print(r)
        cell_id_r=cell_meta[which(cell_meta[,3] == r), 1] %>% as.matrix()
        dat_r <- subsetGiotto(gobject = gobject, cell_ids = cell_id_r )
        res1=ExprExprTable_1region(gobject=dat_r, r=r, cell_meta=cell_meta, select_ligands=select_ligands,
                                   select_receptors=select_receptors,
                             abs_log2fc_LR=abs_log2fc_LR, p_adj = p_adj)
        res=rbind(res, res1)
      }
    }
    if (region_specific==F) {
      res=ExprExprTable_1region(gobject=gobject, r="NULL", cell_meta=cell_meta, select_ligands=select_ligands,
                                select_receptors=select_receptors,
                          abs_log2fc_LR=abs_log2fc_LR, p_adj = p_adj)
    }
    # save results in the  format that can be directly used in sCCIgen
    readr::write_csv(res, file=output_file, col_names = F)

    print(paste0("CCI expression-expression association is saved at: ", output_file))
  }

}

# ExprExprTable_1region ------------
ExprExprTable_1region = function(gobject, r, cell_meta, select_ligands, select_receptors,
                   p_adj=0.05, abs_log2fc_LR=0.25) {

  plan(future::multisession)


  ## get statistical significance of gene pair expression changes upon cell-cell interaction
  sc <- spatCellCellcom(gobject = gobject, spatial_network_name = "spatial_network",
                                          cluster_column = "anno",
                                          random_iter = 500, feat_set_1 = select_ligands,
                                          feat_set_2 = select_receptors,
                                          adjust_method = "fdr",
                                          do_parallel = TRUE, cores = 6,
                                          verbose = "a little")

  ## select top LR table ##
  selected_spat <- sc[sc$p.adj <= p_adj & abs(sc$log2fc) > abs_log2fc_LR & sc$lig_nr >= 5 & sc$rec_nr >= 5,]
  if (nrow(selected_spat) > 0) {
    selected_spat2= selected_spat[, c("lig_cell_type", "rec_cell_type",  "ligand", "receptor", "log2fc")]



    # add an approximate distance threshold by cell tyep pair
    spatial_net <- getSpatialNetwork(gobject = gobject, name = 'spatial_network')
    dis_table=spatial_net@networkDT
    dis_table2 <- merge(dis_table, cell_meta, by.x = "from", by.y = "cell_ID", all.x = TRUE)
    dis_table2 <- merge(dis_table2, cell_meta, by.x = "to", by.y = "cell_ID", all.x = TRUE)
    tb=tapply(dis_table2$distance, list(dis_table2$anno.x, dis_table2$anno.y), median)
    dist_value=reshape2::melt(tb) %>%
      dplyr::filter(is.na(value)==F)%>%
      dplyr::rename(lig_cell_type = Var1, rec_cell_type = Var2, threshold = value)
    df_merged <- selected_spat2 %>% dplyr::inner_join(dist_value, by = c("lig_cell_type", "rec_cell_type"))



    # <Region><Perturbed cell type>,<Adjacent cell type>, <Interaction distance threshold (default 0.1)>,
    # <Gene ID 1 (optional)>, <Gene ID 2 (optional)>,<Gene proportion (optional)>, <Bi-directional association
    # (TRUE or FALSE, default = TRUE)>,<Mean effect at log(count) scale (default = 0.5)>,<SD of effect at log(count)
    # scale (default = 0)>.
    tb=data.frame(r=r,  df_merged[,c("lig_cell_type",  "rec_cell_type", "threshold", "ligand", "receptor")], "NULL",
                  bi_direct=F, df_merged$log2fc, 0 )
    return(tb)

  } else {return(NULL)}

}

#' SpatialTable
#'
#' @param gobject Giotto object.
#' @param top_num Keep top K number of genes for each cell type
#' @param fdr_cut Keep genes whose FDR is less than this cutoff.
#' @param output_file output_file
#' @export
#'
SpatialTable=function(gobject, top_num=2, fdr_cut=0.05,
                      output_file= file.path(getwd(), "est_region_specific_genes.csv")){

  cell_meta=pDataDT(gobject)
  colnames(cell_meta)[3]

  # Separate analysis for each cell type

  topgenes=NULL
  for (c in unique(cell_meta$anno)) {

    #
    cell_id_c=cell_meta[which(cell_meta$anno == c), 1] %>% as.matrix()
    dat_c <- subsetGiotto(gobject = gobject, cell_ids = cell_id_c)
    Rname=colnames(cell_meta)[3]


    markers = findMarkers_one_vs_all(gobject = dat_c,
      method = "scran",
      expression_values = "custom",
      cluster_column =  Rname      # column name in cell metadata
    )

    # Filter by FDR first
    markers_filt <- markers[!is.na(markers$FDR) & markers$FDR < fdr_cut, ]

    # Split by cluster and get top N from each
    split_markers <- split(markers_filt, markers_filt$cluster)
    topgenes1 <- do.call(rbind, lapply(split_markers, function(df) head(df, top_num)))

    # plotMetaDataHeatmap(dat_c,
    #                     selected_feats = unique(topgenes),
    #                     metadata_cols = colnames(cell_meta)[3],
    #                     x_text_size = 10, y_text_size = 10)

    topgenes=rbind(topgenes, data.frame(anno=c, topgenes1))
  }

  tb=data.frame(topgenes[, c("cluster", "anno", "feats")], "NULL", topgenes$logFC, 0)


  # Specify the spatial patterns in the format <Region>,<Cell type>,<Gene ID (optional)>,<Gene proportion (optional)>,
  # <Mean effect at log(count + 1) scale>,<SD of effect at log(count + 1) scale>. NOTE: if you don't provide a Gene ID,
  # the gene proportion must be provided. If you're providing a Gene ID, the proportion will be automatically calculated
  # and you can set the gene proportion = NULL. Separate the entries by blank space
  # (e.g. '1,cell_type_A,NULL,0.1,0.5,0 1,cell_type_A,gene_A,NULL,0.5,0') Alternatively, click the button to select a file
  # separated by commma where each line contains an entry with the format mentioned above.
  # save results in the  format that can be directly used in sCCIgen
  readr::write_csv(tb, file=output_file, col_names = F)

  print(paste0("The estimation of region specific genes are saved at: ", output_file))
}


