
# generate_marginal_params ------
#' Generate model parameters from data.
#'
#' This function generates model parameters from input data for all cell types in regions.
#' @param expr Gene expression level (count).
#' @param feature Cell features (e.g. cell type, spatial coordinates, regions) of reference data.
#' @param region_specific_model Whether estimation model differ in different regions.
#' @param ncores No of cores
#' @return Provide model parameters including marginal distributions and copula
#' (if not NULL) for all cell types in all regions.
#' @export

generate_marginal_params=function(expr,
                                  feature,
                                  region=NULL,
                                  ncores=1) {
  # clean data format
  expr=as.matrix(expr)
  feature=as.matrix(feature)
  if (identical(colnames(expr), feature[,1])==F) {colnames(expr)=feature[,1]}
  Genes=rownames(expr)

  if (region_specific_model!="TRUE") { # not region specific
    cell_type_sel=names(table(colnames(expr)))
    model_params1=  fit_model_scDesign2(data_mat=expr,
                                        cell_type_sel=cell_type_sel,
                                        sim_method = 'ind',
                                        marginal='zinb',
                                        ncores = ncores)


    model_params=list(model_params1)
  }

  #  region specific model
  if (region_specific_model=="TRUE") { #  region specific

    Region=feature[,4]
    Runiq=unique(Region)
    R=length(Runiq)

    model_params= foreach (r = 1:R) %dopar% {
      idx=which(Region==Runiq[r])
      cell_type_sel=names(table(colnames(expr[,idx])))
      model_params1=fit_model_scDesign2(data_mat=expr[,idx],
                                        cell_type_sel=cell_type_sel,
                                        sim_method = 'ind',
                                        marginal='zinb',
                                        ncores = ncores)
      model_params1
    }
  }
  return(model_params)
}

