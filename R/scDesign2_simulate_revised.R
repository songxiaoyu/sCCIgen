
# simulate_count_copula ---------------
#' Simulate a count matrix for a single cell type based on a copula model
#'
#' @param copula_result A list that contains the parameters of a copula model.
#' @param n             An integer value that indicates the number of cells to generate.
#' @param marginal      A character string that indicates whether the generated values should
#'                      stay as discrete or switch to continuous. Default value is 'nb', which
#'                      should be used for generating a count marix. The alternative 'Gamma' is
#'                      only needed when this function is being called by other functions that
#'                      generate data with a user-specified sequencing depth. Normally, users
#'                      do not need to change this value.
#' @return A matrix of shape p by n that contains the simulated count values. p is derived from
#' \code{copula_result}
#'

simulate_count_copula <- function(copula_result, n = 100,
                                  marginal = c('nb', 'Gamma')){
  marginal <- match.arg(marginal)

  p1 <- length(copula_result$gene_sel1)
  if(p1 > 0){
    result1 <- MASS::mvrnorm(n = n, mu = rep(0.0, p1), Sigma = copula_result$cov_mat)
    result1 <- matrix(result1, nrow = n)
    result2 <- apply(result1, 2, pnorm)
    result2 <- matrix(result2, nrow = n)
  }


  p2 <- length(copula_result$gene_sel2)
  if(marginal == 'nb'){
    if(p1 > 0){
      result31 <- t(sapply(1:p1, function(iter){
        param <- copula_result$marginal_param1[iter, ]
        qnbinom(pmax(0.0, result2[, iter] - param[1]) / (1-param[1]),
                size = param[2], mu = param[3]) # simulate negative binomial dist
      }))
    }
    if(p2 > 0){
      result32 <- t(sapply(1:p2, function(iter){
        param <- copula_result$marginal_param2[iter, ]
        rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
      }))
    }
  }else if(marginal == 'Gamma'){
    if(p1 > 0){
      result31 <- t(sapply(1:p1, function(iter){
        param <- copula_result$marginal_param1[iter, ]
        qgamma(max(0.0, result2[, iter] - param[1]), shape = param[2], scale = param[3] / param[2])
      }))
    }
    if(p2 > 0){
      result32 <- t(sapply(1:p2, function(iter){
        param <- copula_result$marginal_param2[iter, ]
        rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
      }))
    }
  }

  result <- matrix(0, nrow = p1 + p2 + length(copula_result$gene_sel3), ncol = n)
  if(p1 > 0){
    result[copula_result$gene_sel1, ] <- result31
  }
  if(p2 > 0){
    result[copula_result$gene_sel2, ] <- result32
  }
  result
}


# simulate_count_ind ---------------
#' Simulate a count matrix for a single cell type based on a (w/o copula model)
#'
#' @param model_params A list that contains the model parameters (can be either the copula
#'                     model or the (w/o copula) model).
#' @inheritParams simulate_count_copula
#' @return A matrix of shape p by n that contains the simulated count values. p is derived from
#' \code{model_params}.
#'

simulate_count_ind <- function(model_params, n = 100,
                               marginal = c('nb', 'Gamma')){
  marginal <- match.arg(marginal)

  if(model_params$sim_method == 'copula' || 'gene_sel3' %in% names(model_params)){

    p1 <- length(model_params$gene_sel1)
    p2 <- length(model_params$gene_sel2)
    if(marginal == 'nb'){
      if(p1 > 0){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }
      if(p2 > 0){
        result32 <- t(sapply(1:p2, function(iter){
          param <- model_params$marginal_param2[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }
    }else if(marginal == 'Gamma'){

      if(p1 > 0){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
      if(p2 > 0){
        result32 <- t(sapply(1:p2, function(iter){
          param <- model_params$marginal_param2[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
    }
    result <- matrix(0, nrow = p1 + p2 + length(model_params$gene_sel3), ncol = n)
    if(p1 > 0){
      result[model_params$gene_sel1, ] <- result31
    }
    if(p2 > 0){
      result[model_params$gene_sel2, ] <- result32
    }
  }else{

    p1 <- length(model_params$gene_sel1)
    p2 <- length(model_params$gene_sel2)
    result <- matrix(0, nrow = p1 + p2, ncol = n)
    if(p1 > 0){
      if(marginal == 'nb'){

        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }else if(marginal == 'Gamma'){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
      result[model_params$gene_sel1, ] <- result31
    }
  }
  result
}



# scDesign2.revised ---------------
#' Simulate a count matrix for a single cell type
#'
#' Revise the initial scDesign2 functions to allow different sequencing depth
#'
#' @param model_params    A list with the same length as \code{cell_type_prop} that contains
#'                        the fitted model as each of its element (can be either the copula
#'                        model or the (w/o copula) model).
#' @param depth_simu_ref_ratio The (expected) sequencing depth ratio between simulated and refernece data.
#' @param n_cell_new      The total number of cells in the simulated count matrix.
#' @param cell_type_prop  The cell type proportion in the simulated count matrix.
#' @param sim_method Simulation method. Simulate genes independently 'ind' or considering
#'                    their correlations ('copula').
#' @return A matrix of shape p by n that contains the simulated count values. p is derived from
#' \code{model_params}.
#'

scDesign2.revised <- function(model_params,
                              n_cell_new,
                              cell_type_prop,
                              depth_simu_ref_ratio = NULL,
                              sim_method = c('copula', 'ind')){


  n_cell_vec <- sapply(model_params, function(x) x$n_cell)
  n_read_vec <- sapply(model_params, function(x) x$n_read)
  total_count_old <- sum(n_read_vec)
  n_cell_old      <- sum(n_cell_vec)

  if(length(model_params)!=length(cell_type_prop)){
    stop('Cell type proportion should have the same length as the number of models.')
  }

  n_cell_type <- length(cell_type_prop)
  cell_type_prop <- cell_type_prop / sum(cell_type_prop)
  n_cell_each <- round(cell_type_prop * n_cell_new)
  if(sum(n_cell_each) != n_cell_new){
    idx <- sample(n_cell_type, size = 1)
    n_cell_each[idx] <- n_cell_each[idx] + n_cell_new - sum(n_cell_each)
  }


  p <- length(model_params[[1]]$gene_sel1) + length(model_params[[1]]$gene_sel2) +
    length(model_params[[1]]$gene_sel3)
  new_count <- matrix(0, nrow = p, ncol = n_cell_new)


  r <- rep(depth_simu_ref_ratio, n_cell_type)

  for(iter in 1:n_cell_type)
    if(n_cell_each[iter] > 0){
      ulim <- sum(n_cell_each[1:iter])
      llim <- ulim - n_cell_each[iter] + 1
      params_new <- model_params[[iter]]
      params_new$marginal_param1[, 3] <- params_new$marginal_param1[, 3] * r[iter]
      params_new$marginal_param2[, 3] <- params_new$marginal_param2[, 3] * r[iter]
      if(sim_method == 'copula'){
        new_count[, llim:ulim] <- simulate_count_copula(copula_result=params_new,
                                                        n = n_cell_each[iter],
                                                        marginal = 'nb')
      }else if(sim_method == 'ind'){
        new_count[, llim:ulim] <- simulate_count_ind(model_params=params_new,
                                                     n = n_cell_each[iter],
                                                     marginal = 'nb')
      }

    }
  if(is.null(names(model_params))){
    colnames(new_count) <- unlist(lapply(1:n_cell_type, function(x){rep(x, n_cell_each[x])}))
  }else{
    colnames(new_count) <- unlist(lapply(1:n_cell_type, function(x){
      rep(names(model_params)[x], n_cell_each[x])}))
  }
  return(new_count)

}



