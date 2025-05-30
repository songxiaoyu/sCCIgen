
# simu.window ---------------
#' Simulate the window of spatial data
#'
#' This function esimates window using spatial location data of existing cells.
#' @param PointLoc The location of input cells on x, y axis.
#' @param method Options include method=c("network", "convex", "convex2", "convex3",
#' "convex5", "rectangle"). Network is preferred unless the input data is very large.
#' @return The spatial window of input cells.
#' @export
# window=simu.window(PointLoc=NULL)
# window=simu.window(PointLoc=PointLoc, method="convex5")
simu.window=function(PointLoc,
                     method=c("network", "convex", "convex2", "convex3",
                              "convex5", "rectangle")) {

  x0=PointLoc[,1]
  y0=PointLoc[,2]
  n=nrow(PointLoc)

  if (is.null(PointLoc)==F &(method=="convex" | method=="rectangle")) {
    res=spatstat.geom::ripras(x0, y0, shape=method)
  }


  if (is.null(PointLoc)==F & method %in% c("convex2", "convex3", "convex5") ) {
    K=as.numeric(substr(method, 7,7))
    xmin=min(x0);xmax=max(x0);ymin=min(y0);ymax=max(y0)
    dx=xmax-xmin
    dy=ymax-ymin
    dmax=max(dx, dy)
    res=list()
    for (i in 1: K) {
      idx= which(x0>=(xmin+ dx/K*(i-1)) & x0< (xmin + dx/K*i*1.05))
      for (j in 1:K) {
        idy= which(y0>=(ymin+ dy/K*(j-1)) & y0< (ymin + dy/K*j*1.05))
        id=intersect(idx, idy)
        res1=spatstat.geom::ripras(x0[id], y0[id], shape="convex")
        res=suppressWarnings(spatstat.geom::union.owin(res, res1))
      }
    }
  }

  if (is.null(PointLoc)==F & method=="network") {

    # generate  perturbation

    dx=max(x0)-min(x0)
    dy=max(y0)-min(y0)
    dmax=max(dx, dy)
    x=x0+0.005*dmax*cos(seq(0, 2*pi, length=n*4))
    y=y0+0.005*dmax*sin(seq(0, 2*pi, length=n*4))
    PointLoc_noise=cbind(x,y)
    delaunay_triangle = geometry::delaunayn(PointLoc_noise)
    max_dist=sapply(1:nrow(delaunay_triangle), function(f)
      max(stats::dist(PointLoc_noise[delaunay_triangle[f,],])))
    #trimesh(delaunay_triangle,cbind(x,y))

    ## get rid of the edges that connect distant points
    max_dist_standize=max_dist/max(max_dist)
    log_dist=log(max_dist_standize+1)
    idx=which(log_dist> mean(log_dist) + 3*sd(log_dist))
    delaunay_triangle2=delaunay_triangle[-idx,]
    #trimesh(delaunay_triangle2,PointLoc_noise)

    ## get the outer window
    delaunay_edges <- rbind(delaunay_triangle2[ ,c(1,2)],
                            delaunay_triangle2[ ,c(1,3)],
                            delaunay_triangle2[ ,c(2,3)])
    delaunay_edges_ordered=cbind(apply(delaunay_edges, 1, min),
                                 apply(delaunay_edges, 1, max))
    dup1= duplicated(delaunay_edges_ordered)
    dup2=duplicated(delaunay_edges_ordered, fromLast = TRUE)
    outer_edges=delaunay_edges[which(dup1==F & dup2 ==F),]

    points_ordered=outer_edges[1,]
    edges_remained=outer_edges[-1,]
    n=nrow(edges_remained)
    for (i in 1:(n-1)) {
      #print(i)
      idx=which(apply(edges_remained, 1, function(r)
        any(r==points_ordered[i+1])))
      points_ordered=c(points_ordered,
                       setdiff(edges_remained[idx,], points_ordered[i+1]))
      edges_remained=edges_remained[-idx,]
      #print(points_ordered)
    }

    a <- try(spatstat.geom::owin(poly=list(x=PointLoc_noise[points_ordered,1],
                            y=PointLoc_noise[points_ordered,2])), silent =T)
    if (inherits(a, "try-error")) {
      points_ordered2=rev(points_ordered)
      res=spatstat.geom::owin(poly=list(x=PointLoc_noise[points_ordered2,1],
                         y=PointLoc_noise[points_ordered2,2]))
    } else {res=a}

  }
  return(res)
}



# cell.loc.1region.model.fc ---------------
# Simulate ST data for one region location based on parametric model
cell.loc.1region.model.fc=function(n,
                           PointLoc,
                           PointAnno,
                           window_method,
                           seed=NULL) {

  #
  if(is.null(seed)==F) {set.seed(seed)}
  cell_win=simu.window(PointLoc=PointLoc, method=window_method)
  p=spatstat.geom::as.ppp(PointLoc, W=cell_win)
  spatstat.geom::marks(p)=as.factor(PointAnno)
  # if too many cells
  if (p$n>10000) {
    idx=rbinom(p$n, 1, prob=10000/p$n)
    p2=subset(p, idx==1)
    p=p2
  }

  x=p$x
  y=p$y
  fit=spatstat.model::ppm(p, ~marks+ marks:polynom(x,y,2),
                          spatstat.model::Poisson())
  nsim=ceiling(n/p$n)
  if (nsim>1) {
    a=spatstat.random::rmh(model=fit, nsim=nsim)
    b=spatstat.geom::superimpose(a)
    spatstat.geom::marks(b)=spatstat.geom::marks(b)[,1]
  } else{b=spatstat.random::rmh(model=fit, nsim=nsim)}

  # get rid of cells on the same location
  while(b$n-n>10) {
    delete.n=(b$n-n)
    dis=spatstat.geom::pairdist(b)
    dis[lower.tri(dis, diag=T)]=NA
    r= quantile(dis, probs=delete.n*2/b$n/(b$n-1), na.rm=T)
    dis2=dis< r
    same.loc.idx=which(dis2 == T, arr.ind = TRUE)
    b=b[setdiff(1:b$n, same.loc.idx[,1]), ]
  }

  return(b)
}




# cell.loc.model.fc ---------------
#' Generate cell location data by modeling the spatial information of existing data.
#' @param n No. of cells
#' @param PointLoc The location of input cells on x, y axis.
#' @param PointAnno The cell type annotation of input cells.
#' @param PointRegion The spatial regions of input cells.
#' @param window_method Method for estimating window of cells.
#' @param seed Random seed.
#' @import spatstat
#' @export
cell.loc.model.fc=function(n,
                           PointLoc,
                           PointAnno,
                           PointRegion,
                           window_method,
                           seed=NULL) {
  Rcat=unique(PointRegion)
  bb=vector("list", length(Rcat))
  names(bb)=Rcat
  for ( i in 1:length(Rcat)) {
    idx= which(PointRegion %in% Rcat[i])
    n.sim.region=round(n*length(idx)/length(PointAnno))
    bb[[i]]=cell.loc.1region.model.fc(n=n.sim.region,
                              PointLoc=PointLoc[idx,],
                              PointAnno=PointAnno[idx],
                             window_method=window_method,
                             seed=seed)
  }

  return(bb)
}



# cell.loc.1region.existing.fc ---------------
# simulate ST data location for one region using existing cell location
cell.loc.1region.existing.fc=function(PointLoc,
                           PointAnno,
                           window_method="rectangle") {


  cell_win=simu.window(PointLoc=PointLoc, method=window_method)
  p=spatstat.geom::as.ppp(PointLoc, W=cell_win)
  spatstat.geom::marks(p)=as.factor(PointAnno)

  return(p)
}


# cell.loc.existing.fc ---------------

#' Generate cell location data by using existing SRT data directly.
#' @param PointLoc The location of input cells on x, y axis.
#' @param PointAnno The cell type annotation of input cells.
#' @param PointRegion The spatial regions of input cells.
#' @param window_method Method for estimating window of cells.
#' @import spatstat
#' @export
cell.loc.existing.fc=function(PointLoc,
                              PointAnno,
                              PointRegion,
                              window_method="rectangle") {
  Rcat=unique(PointRegion)
  pp=vector("list", length=length(Rcat))
  for (i in 1:length(Rcat)) {
    idx= which(PointRegion %in% Rcat[i])
    pp[[i]]=cell.loc.1region.existing.fc(PointLoc=PointLoc[idx,],
                         PointAnno=PointAnno[idx],
                         window_method=window_method)
  }
  names(pp)=Rcat

  return(pp)
}
