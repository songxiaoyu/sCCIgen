################### List of Functions #########################
# Est_ModelPara
# generate_marginal_params
# Use_scDesign2_1region
# Use_scDesign2
# Find.Neighbor.Pairs
# Add.Spatial.Expr.Pattern
# Add.Distance.Asso.Pattern
# Add.Expr.Asso.Pattern
# ExprPattern
# Pattern.adj.1region
# Pattern.Adj
# MergeRegion
##########################################




# Est_ModelPara ---------------
#' Estimate Marginal distribution and Gaussian Copula for Gene Expression Matrix
#' @param expr Expression levels of input data
#' @param anno Cell type annotation of input data
#' @param ncores No of cores for parallel computing.
#' @param region description
#' @param sim_method c('ind', 'copula')
#' @return Estimated Gaussian Copula
#' @export

Est_ModelPara <- function(expr, anno, sim_method=c('ind', 'copula'), region=NULL, ncores = 1) {
  expr=as.matrix(expr)

  if (is.null(region)==T) {
    ct=names(table(colnames(expr)))
    copula1=fit_model_scDesign2(data_mat=expr,
                               cell_type_sel=ct, sim_method =sim_method,
                               marginal = 'zinb',
                               ncores = ncores)
    copula=list(copula1)
  } else {
    R=unique(region)%>% as.character()
    copula=vector(mode = "list", length = length(R))
    names(copula)=R
    for (r in R) {
      expr2=expr[, region==r]
      ct=names(table(colnames(expr2)))
      copula1=fit_model_scDesign2(data_mat=expr2,
                                 cell_type_sel=ct, sim_method = sim_method,
                                 marginal = 'zinb',
                                 ncores = ncores)
      copula[[r]]=copula1
    }

  }

  return(copula)

}



# Use_scDesign2_1region ------
Use_scDesign2_1region=function(ppp.obj1, Genes, model_params,
                       depth_simu_ref_ratio=1, cell_type_sel, seed,
                       sim_method = c('copula', 'ind')) {
  # cell types in simulated and reference data
  n.ordered=table(ppp.obj1$marks)
  exist.cell.type=names(n.ordered)
  cell_type_prop=n.ordered/ppp.obj1$n
  model_params_exist=model_params[exist.cell.type]
  set.seed(seed)
  sim_count <- scDesign2.revised(model_params=model_params_exist,
                                 n_cell_new=ppp.obj1$n,
                                 cell_type_prop = cell_type_prop,
                                 depth_simu_ref_ratio=depth_simu_ref_ratio,
                                 sim_method =sim_method)

  # Update the order of sim_count to match the cell type of ppp.obj1
  sim_count2=matrix(NA, ncol=ncol(sim_count), nrow=nrow(sim_count))
  colnames(sim_count2)=ppp.obj1$marks
  rownames(sim_count2) = Genes
  for (f in levels(ppp.obj1$marks)) {
    sim_count2[, which(colnames(sim_count2)==f)]=
      sim_count[,which(colnames(sim_count)==f)]
  }

  return(sim.count=sim_count2)
}



# Use_scDesign2 ------
#' Generate expression profile under no spatial patterns.
#'
#' This function uses input parameters to simulate expression for all regions.
#' @param ppp.obj Cells as points on a spatial map for all regions.
#' @param model_params Provide model parameters including marginal distributions and copula (if not NULL).
#' @param expr Gene expression (count) in reference data.
#' @param region Cell regions in reference data.
#' @param depth_simu_ref_ratio Relative sequencing depth in comparison of reference data.
#' @param sim_method Simulate independent genes using'ind' or correlated genes using 'copula'.
#' @param region_specific_model Whether estimation model differ in different regions.
#' @param seed Seed
#' @return Simulated expression count data for all cells in all regions.
#' @export

Use_scDesign2=function(ppp.obj,
                       model_params,
                       expr,
                       region=NULLL,
                       depth_simu_ref_ratio=1,
                       sim_method = c('copula', 'ind'),
                       region_specific_model,
                       seed) {

  expr=as.matrix(expr)
  R=length(ppp.obj) # simulated regions
  cell_type_sel=names(table(colnames(expr)))
  Genes=rownames(expr)


  if (region_specific_model!="TRUE") { # not region specific


    if (R==1) {

      sim.count =list(Use_scDesign2_1region(ppp.obj1=ppp.obj[[1]],
                              Genes=Genes,
                              model_params=model_params[[1]],
                              depth_simu_ref_ratio=depth_simu_ref_ratio,
                              cell_type_sel=cell_type_sel,
                              seed=seed*31+1*931,
                              sim_method = sim_method))
      } else { # end if R==1

        sim.count= foreach (r = 1:R) %dopar%{
          Use_scDesign2_1region(ppp.obj1=ppp.obj[[r]],
                                Genes=Genes,
                                model_params=model_params[[1]],
                                depth_simu_ref_ratio=depth_simu_ref_ratio,
                                cell_type_sel=cell_type_sel,
                                seed=seed*31+r*931,
                                sim_method = sim_method)
          }
        } # end else of R==1
    } # end not region specific

  #  region specific model
  if (region_specific_model=="TRUE") { #  region specific

    Runiq=unique(region)

    sim.count= foreach (r = 1:length(Runiq)) %dopar% {

      Use_scDesign2_1region(ppp.obj1=ppp.obj[[r]],
                                           Genes=Genes,
                                           model_params=model_params[[r]],
                                           depth_simu_ref_ratio=depth_simu_ref_ratio,
                                           cell_type_sel=cell_type_sel,
                                           seed=seed*31+r*931,
                                           sim_method = sim_method)
    }

  }

  return(sim.count)
}

# Find.Neighbor.Pairs ----------

Find.Neighbor.Pairs=function(ppp.obj,
                             interacting.cell.type.pair,
                             int.dist.threshold) {
  cell.loc=cbind(ppp.obj$x, ppp.obj$y)
  cell1.idx=which(ppp.obj$marks==interacting.cell.type.pair[1])
  cell2.idx=which(ppp.obj$marks==interacting.cell.type.pair[2])
  m=spatstat.geom::crossdist(cell.loc[cell1.idx,1], cell.loc[cell1.idx,2],
                             cell.loc[cell2.idx,1],cell.loc[cell2.idx,2])
  # in neighbor or not?
  dmax=max( max(ppp.obj$x)-min(ppp.obj$x), max(ppp.obj$y)-min(ppp.obj$y))

  neighbo.loc.idx=which(m< (int.dist.threshold*dmax), arr.ind = TRUE)

  # index in original data
  nbr.idx=cbind(cell1.idx[neighbo.loc.idx[,1]],   cell2.idx[neighbo.loc.idx[,2]])

  colnames(nbr.idx)=interacting.cell.type.pair
  return(neighbor.idx=nbr.idx)
}

# Add.Spatial.Expr.Pattern -----------------------------------------------
#' Adds spatial differential expressed pattern (region-specific effects) to
#' a cell type.
#'
#' This function adds one type of spatial differential expressed patterns. This
#' function can be repeated used to add region-specific effects  in different
#' regions for different cell types.
#' @param sim.count Cells as points on a spatial map.
#' @param r Region index.
#' @param CellType Cell type index.
#' @param GeneID Gene(s) index. Default = NULL, and then a random subset
#' of genes will be perturbed based on the defined spatial patterns.
#' @param PropOfGenes Proportion of genes with this pattern if GeneID is not
#' provided.
#' @param delta.mean Mean effect (on the log count scale).
#' @param delta.sd SD of effect.
#' @param seed Seed
#' @return
#' \item{SignalSummary:}{Summary of this spatial pattern, including the type
#' of spatial patterns, impacted cell types, perturbed genes, and effect sizes.}
#' \item{beta.matrix:}{Effect size on each gene in each cell. }

#' @export
#'
Add.Spatial.Expr.Pattern= function(sim.count,
                                   r,
                                   CellType,
                                   GeneID=NULL,
                                   PropOfGenes=0.1,
                                   delta.mean=1,
                                   delta.sd=0.01, seed) {
  set.seed(938*seed-142)
  R=length(sim.count)
  idx=which(names(sim.count)==r)
  sim.count1=sim.count[[idx]]
  G=nrow(sim.count1)
  N=ncol(sim.count1)
  GeneAll=rownames(sim.count1)

  # key matrix
  beta.matrix=vector("list", length=R)
  for (i in 1:R) {beta.matrix[[i]]=matrix(0, nrow=G, ncol=ncol(sim.count[[i]]))}
  colnames(beta.matrix[[idx]])=colnames(sim.count1)
  rownames(beta.matrix[[idx]])=GeneAll

  # GeneID
  if (is.null(GeneID)) {GeneID=sample(GeneAll, round(PropOfGenes * G))}

  CellID=which(colnames(beta.matrix[[idx]])== CellType )
  beta=stats::rnorm(length(GeneID), delta.mean, delta.sd)
  SignalSummary=data.frame(Type="SpatialChange", Region=r, CellType, GeneID,
                           AdjCellType="NA",
                           AdjGene="NA", beta)
  beta.matrix[[idx]][GeneID,CellID] = beta +  beta.matrix[[idx]][GeneID,CellID]

  return(list(SignalSummary=SignalSummary, beta.matrix=beta.matrix))
}



MergePPP=function(points.list) {
  K=length(points.list)

  # points
  x.combine=unlist(lapply(1:K, function(f) points.list[[f]]$x))
  y.combine=unlist(lapply(1:K, function(f) points.list[[f]]$y))
  annotation=unlist(lapply(1:K, function(f) points.list[[f]]$marks))

  win1=points.list[[1]]$window
  if (K>1) {
    for (k in 2:K) {
      win1=spatstat.geom::union.owin(win1, points.list[[k]]$window)
    }
  }


  points1=spatstat.geom::ppp(x.combine, y.combine, window=win1)
  spatstat.geom::marks(points1)=annotation

  return(points1)
}


# Add.Distance.Asso.Pattern -----------------------------------------------
#' Add cell-cell expr-distance interaction to a pair of cell types
#'
#' This function add a type of cell-cell interactions to a pair of cell types:
#' the expression in a cell type associated with the proximity of
#' the other cell type. One can repeat this function for multiple times to
#' add cell-cell interactions for many cell type pairs and regions.
#' @param ppp.obj An object of class "ppp" representing simulated cell locations
#' @param sim.count Simulated expression counts from single-cell expression
#' data, before adding in additional spatial patterns.
#' @param r Which region to add in the spatial pattern.
#' @param perturbed.cell.type Which cell type is perturbed from this cell-cell
#' interaction (e.g. microglia).
#' @param adjacent.cell.type Which cell type in the neighbor perturbs from
#' the cell-cell interaction (e.g. neuron).
#' @param int.dist.threshold The minimal cell-cell distance for the interaction.
#' @param delta.mean Expected effect.size (at the log scale of the counts).
#' @param delta.sd Standard deviation of the effect.size
#' @param GeneID Affected genes.
#' @param PropOfGenes Proportion of genes impacted by the cell-cell interaction.
#' It is used if GenePairIDMatrix is NULL, and a random subset of genes with
#' specified proportion will be perturbed.
#' @param seed Seed
#' @return
#' \item{SignalSummary:}{Summary of this spatial pattern, including the type
#' of spatial patterns, impacted cell types, perturbed genes, and effect sizes.}
#' \item{beta.matrix:}{Effect size on each gene in each cell. }
#' @export

Add.Distance.Asso.Pattern = function(ppp.obj,
                                   sim.count, r,
                                   perturbed.cell.type,
                                   adjacent.cell.type,
                                   int.dist.threshold=0.1,
                                   delta.mean=1,
                                   delta.sd=0.001,
                                   GeneID=NULL, # Cell A Gene 1--> Cell B
                                   PropOfGenes=NULL,
                                   seed=NULL) {


  set.seed(seed*478-50194)

  R=length(sim.count)
  GeneAll=rownames(sim.count[[1]])
  G=length(GeneAll)
  N=sapply(1:R, function(f) ncol(sim.count[[f]]))

  # key matrix
  beta.matrix=vector("list", length=R)
  for (i in 1:R) {
    beta.matrix[[i]]=matrix(0, nrow=G,  ncol=ncol(sim.count[[i]]))
    colnames(beta.matrix[[i]])=colnames(sim.count[[i]])
    rownames(beta.matrix[[i]])=GeneAll
  }


  # GeneID
  if (is.null(GeneID)) {GeneID=sample(GeneAll, round(PropOfGenes * G))}
  beta=stats::rnorm(length(GeneID), delta.mean, delta.sd)

  if (r == "NULL") {  # add interaction in all regions
    N0=c(0, N)
    for (i in 1:R) {
      # spatial info
      ppp.obj.combined=MergePPP(ppp.obj)
      nbr.idx=Find.Neighbor.Pairs(ppp.obj=ppp.obj.combined,
                                  interacting.cell.type.pair=c(perturbed.cell.type,
                                                               adjacent.cell.type),
                                  int.dist.threshold=int.dist.threshold)
      N_b=sum(N0[1:i])
      N_a=sum(N0[1:(i+1)])
      idx=seq(N_b+1, N_a)
      idx1=nbr.idx[,1][which(nbr.idx[,1] %in% idx)]-N_b
      temp= beta +  beta.matrix[[i]][GeneID, idx1]
      colnames(temp)=idx1
      temp2=sapply(1:nrow(temp), function(f) tapply(temp[f,], idx1, sum))
      beta.matrix[[i]][GeneID, as.numeric(rownames(temp2))] =t(temp2)
    }
  } else {  # add interaction in one region
    idx_r=which(names(sim.count)==r)
    nbr.idx=Find.Neighbor.Pairs(ppp.obj=ppp.obj[[idx_r]],
                                interacting.cell.type.pair=c(perturbed.cell.type,
                                                             adjacent.cell.type),
                                int.dist.threshold=int.dist.threshold)
    idx1=nbr.idx[,1]
    temp= beta +  beta.matrix[[idx_r]][GeneID,idx1]
    colnames(temp)=idx1
    temp2=sapply(1:nrow(temp), function(f) tapply(temp[f,], idx1, sum))
    beta.matrix[[i]][GeneID, as.numeric(rownames(temp2))] =t(temp2)

  }

  SignalSummary=data.frame(Type="DistanceAssoGenes",
                           Region=r,
                           CellType=perturbed.cell.type,
                           GeneID,
                           AdjCellType=adjacent.cell.type,
                           AdjGene="NA",
                           beta)


  return(list(SignalSummary=SignalSummary, beta.matrix=beta.matrix))
}

# Add.Expr.Asso.Pattern --------
#' Add cell-cell expr-expr interaction pattern to a pair of cell types
#'
#' This function add cell-cell interactions to a pair of cell types (e.g.
#' neuron-microglia) for expression in a cell type associated with expression of
#' the neighboring other cell type. One can repeat this function for
#' multiple times to add cell-cell interactions for many cell types.
#' @param ppp.obj An object of class "ppp" representing simulated cell locations
#' @param sim.count Simulated expression counts from single-cell expression
#' data, before adding in additional spatial patterns.
#' @param r Which region to add in the spatial pattern. If simulated data do
#' not have multiple regions, r=1.
#' @param perturbed.cell.type Which cell type is perturbed from this cell-cell
#' interaction (e.g. microglia).
#' @param adjacent.cell.type Which cell type in the neighbor perturbs from
#' the cell-cell interaction (e.g. neuron).
#' @param Bidirectional Whether the perturbation is both directional.
#' @param int.dist.threshold The minimal cell-cell distance for the interaction.
#' @param delta.mean Expected effect.size (at the log scale of the counts).
#' @param delta.sd Standard deviation of the effect.size
#' @param GenePairIDMatrix Affected gene pairs.
#' @param PropOfGenes Proportion of genes impacted by the cell-cell interaction.
#' It is used if GenePairIDMatrix is NULL, and a random subset of genes with
#' specified proportion will be perturbed.
#' @param seed Seed
#' @return
#' \item{SignalSummary:}{Summary of this spatial pattern, including the type
#' of spatial patterns, impacted cell types, perturbed genes, and effect sizes.}
#' \item{beta.matrix:}{Effect size on each gene in each cell. }

Add.Expr.Asso.Pattern = function(ppp.obj, sim.count, r,
                           perturbed.cell.type,
                           adjacent.cell.type,
                           Bidirectional=T,
                           int.dist.threshold=0.1,
                           delta.mean=1,
                           delta.sd=0.001,
                           GenePairIDMatrix=NULL,
                           PropOfGenes=NULL,
                           seed=NULL) {


  set.seed(seed*3+194)

  R=length(sim.count)
  GeneAll=rownames(sim.count[[1]])
  G=length(GeneAll)
  N=sapply(1:R, function(f) ncol(sim.count[[f]]))

  # key matrix
  beta.matrix=vector("list", length=R)
  for (i in 1:R) {
    beta.matrix[[i]]=matrix(0, nrow=G,  ncol=ncol(sim.count[[i]]))
    colnames(beta.matrix[[i]])=colnames(sim.count[[i]])
    rownames(beta.matrix[[i]])=GeneAll
  }

  # GeneID
  if (is.null(GenePairIDMatrix)) {
    GenePairIDMatrix=matrix(sample(GeneAll, 2*round(PropOfGenes * G)),ncol=2)
    }
  # effect size
  beta=stats::rnorm(nrow(GenePairIDMatrix), delta.mean, delta.sd)

  if (r == "NULL") { # add interaction to all regions
    # spatial info
    ppp.obj.combined=MergePPP(ppp.obj)
    sim.count.combined=rlist::list.cbind(sim.count)
    beta.matrix.combined=rlist::list.cbind(beta.matrix)
    nbr.idx=Find.Neighbor.Pairs(ppp.obj=ppp.obj.combined,
                                interacting.cell.type.pair=c(perturbed.cell.type,
                                                             adjacent.cell.type),
                                int.dist.threshold=int.dist.threshold)

    # 1 --> 2
    count2=sim.count.combined[GenePairIDMatrix[,2], nbr.idx[,2]]

    beta.matrix.combined[GenePairIDMatrix[,1], nbr.idx[,1]]=
      beta.matrix.combined[GenePairIDMatrix[,1], nbr.idx[,1]]+
      beta*log(count2+1)

    SignalSummary=data.frame(Type="ExprAssoGenes", Region=r, CellType=perturbed.cell.type,
                             GeneID=GenePairIDMatrix[,1], AdjCellType=adjacent.cell.type,
                             AdjGene=GenePairIDMatrix[,2], beta)
    # 2 --> 1
    if (Bidirectional==T) {
      count1=sim.count.combined[GenePairIDMatrix[,1], nbr.idx[,1]]

      beta.matrix.combined[GenePairIDMatrix[,2], nbr.idx[,2]]=
        beta.matrix.combined[GenePairIDMatrix[,2], nbr.idx[,2]]+
        beta*log(count1+1)

      SignalSummary=rbind(SignalSummary,
                          data.frame(Type="ExprAssoGenes", Region=r, CellType=adjacent.cell.type ,
                                     GeneID=GenePairIDMatrix[,2], AdjCellType=perturbed.cell.type,
                                     AdjGene=GenePairIDMatrix[,1], beta))
    }
    N0=c(0, N)
    for (i in 1:R) {
      N_b=sum(N0[1:i])+1
      N_a=sum(N0[1:(i+1)])
      beta.matrix[[i]]=beta.matrix.combined[, (N_b:N_a)]
    }


  } else { # add interaction to one region
    idx_r=which(names(sim.count)==r)
    # spatial info
    nbr.idx=Find.Neighbor.Pairs(ppp.obj=ppp.obj[[idx_r]],
                                interacting.cell.type.pair=c(perturbed.cell.type, adjacent.cell.type),
                                int.dist.threshold=int.dist.threshold)


    idx1=nbr.idx[,1]
    idx2=nbr.idx[,2]
    # 1 --> 2
    count2=sim.count[[idx_r]][GenePairIDMatrix[,2], idx2]
    temp=beta.matrix[[idx_r]][GenePairIDMatrix[,1], idx1]+beta*log(count2+1)
    colnames(temp)=idx1
    temp2=sapply(1:nrow(temp), function(f) tapply(temp[f,], idx1, sum))
    beta.matrix[[idx_r]][GeneID, as.numeric(rownames(temp2))] =t(temp2)

    SignalSummary=data.frame(Type="ExprAssoGenes", Region=r, CellType=perturbed.cell.type,
                             GeneID=GenePairIDMatrix[,1], AdjCellType=adjacent.cell.type,
                             AdjGene=GenePairIDMatrix[,2], beta)


    # 2 --> 1
    if (Bidirectional==T) {
      count1=sim.count[[idx_r]][GenePairIDMatrix[,1], idx1]
      temp=beta.matrix[[idx_r]][GenePairIDMatrix[,1], idx2]+beta*log(count1+1)
      colnames(temp)=idx2
      temp2=sapply(1:nrow(temp), function(f) tapply(temp[f,], idx2, sum))
      beta.matrix[[idx_r]][GeneID, as.numeric(rownames(temp2))] =t(temp2)

      SignalSummary=rbind(SignalSummary,
                          data.frame(Type="ExprAssoGenes", Region=r, CellType=adjacent.cell.type ,
                                     GeneID=GenePairIDMatrix[,2], AdjCellType=perturbed.cell.type,
                                     AdjGene=GenePairIDMatrix[,1], beta))
    }
    }

  return(list(SignalSummary=SignalSummary, beta.matrix=beta.matrix))
}


# ExprPattern --------
ExprPattern=function(pattern.list.i){
  L=length(pattern.list.i)
  res=NULL
  for (l in 1:L) {
    res=rbind(res, pattern.list.i[[l]]$SignalSummary)
  }
  return(res)
}
# Pattern.adj.1region --------

Pattern.adj.1region= function(sim.count1, combined.beta.matrix,
                    bond.extreme=T, keep.total.count=F,
                    integer=T) {
  # w vs w/o combined.beta.matrix
  if (is.null(combined.beta.matrix)) {sim.count1.update=sim.count1} else{
    sim.count1.update= exp(log(sim.count1+1) +combined.beta.matrix)-1
  }

  # bond extreme values
  if (bond.extreme==T) {

    cut=quantile(sim.count1.update, probs=0.975, na.rm=T)*10
    sim.count1.update[which(sim.count1.update>cut)]=cut

    ten_pct=apply(sim.count1.update, 2, mean, na.rm=T)*100
    for (i in 1:ncol(sim.count1.update)) {
      sim.count1.update[,i][which(sim.count1.update[,i]>ten_pct[i])]=ten_pct[i]
    }
  }

  # scale by total count
  if (keep.total.count==T) {
    ratio=sum(sim.count1, na.rm=T)/sum(sim.count1.update, na.rm=T)
    sim.count1.update=sim.count1.update*ratio
  }

  # integer
  if (integer==T) {
    sim.count1.update=round(sim.count1.update)
  }

  return(sim.count1.update)
}


# Pattern.Adj --------
#' Adjust the count data for all cells in all regions based on the
#' input spatial patterns
#'
#' Adjust the count data for all regions based on the input spatial patterns
#' @param sim.count Spatial info for cell type 1 (e.g. neuron)
#' @param pattern.list A list of spatial patterns, which can be generated
#' from `Add.Spatial.Expr.Pattern`, `Add.Distance.Asso.Pattern`, or `Add.Expr.Asso.Pattern`.
#' @param bond.extreme Whether to bond extreme high values generated from
#' large effect sizes of spatial patterns (default = TRUE). If TRUE, no
#' gene can have more than 10% of the total counts across all genes.
#' @param keep.total.count If additional spatial patterns are added, whether
#' to rescale expression levels of all genes to keep the sequencing depth
#' (default = FALSE).
#' @param integer Whether to keep counts as integer (default=TRUE).
#' @return Updated simulated gene expression counts by taking into consideration
#' of spatial patterns.
#' @export

Pattern.Adj= function(sim.count, pattern.list=NULL,
                            bond.extreme=T, keep.total.count=F,
                            integer=T) {
  R=length(sim.count)

  sim.count.update=vector("list", length=R)
  for (i in 1:R) {
     if (is.null(pattern.list)) {
       combined.beta.matrix=NULL
     } else {
       K=length(pattern.list)
       beta.matrix.list=lapply(1:K, function(k)
         pattern.list[[k]]$beta.matrix[[i]])
       combined.beta.matrix=Reduce("+", beta.matrix.list)
     }

     sim.count.update[[i]]=Pattern.adj.1region(sim.count1=sim.count[[i]],
                                                combined.beta.matrix=combined.beta.matrix,
                                                bond.extreme=bond.extreme, keep.total.count=keep.total.count,
                                                integer=integer)

  }

  return(sim.count.update=sim.count.update)
}



# MergeRegion --------
#' Merge spatial and expression data from multiple regions
#'
#' Merge spatial and expression data from multiple regions
#' @param points.list points.list is a list of points from multiple regions
#' @param expr.list expr.list is a list of expressions from multiple regions
#' @return
#' \item{meta:}{meta}
#' \item{count:}{count}
#' @import spatstat
#' @export
#'
#'


MergeRegion=function(points.list, expr.list) {
  K=length(points.list)
  # points
  x.combine=unlist(lapply(1:K, function(f) points.list[[f]]$x))
  x.combine2=round(x.combine, digits=4)
  y.combine=unlist(lapply(1:K, function(f) points.list[[f]]$y))
  y.combine2=round(y.combine, digits=4)
  annotation=unlist(lapply(1:K, function(f) points.list[[f]]$marks))
  Cell=paste0("Cell", seq(1, length(x.combine)))

  n=sapply(1:K, function(f) points.list[[f]]$n)
  if (K>1) {
    region=rep(names(points.list), times=n)
    meta=data.frame(Cell=Cell, annotation=annotation,  x.loc=x.combine2,
                    y.loc=y.combine2,
                    region=region)
  } else {
    meta=data.frame(Cell=Cell, annotation=annotation,
                    x.loc=x.combine2, y.loc=y.combine2)
  }

  # expr
  expr.combine=Reduce(cbind, expr.list)
  colnames(expr.combine)=Cell

  return(list(meta=meta,count=as.data.frame(expr.combine)))
}


