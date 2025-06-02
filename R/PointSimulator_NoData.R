# connectUp ------------
#' Assign connected regions in a window
#'
#' This function assigns random connected regions on a square. Used within function `RandomRegionWindow`.
#' @param nRegion nRegion is the No. of regions (e.g. nRegion=3)
#' @param r poly is a RasterLayer  (e.g. 20 by 20 square).
#' @param seed Random seed.
#' @return A list of the selected polygons for each region.


connectUp <- function(r, nRegion, seed=NULL){
  if (is.null(seed)==F) {set.seed(seed)}

  n=length(r)
  nb = lapply(1:n, function(f) raster::adjacent(r, cells=f)[,2]) # find neighbors
  g = igraph::graph.adjlist(nb) # change to network
  nb.length=sapply(nb, length)
  selected = as.list(sample(which(nb.length<4),nRegion))
  selected.sum=unlist(selected)
  n0=length(selected.sum)

  while(n0 < n){
    for (p in 1:nRegion) {
      nbrs = unlist(igraph::ego(g, 1, selected[[p]], mindist=1)) # find neighbors
      selected.sum=unlist(selected)
      newnbrs = nbrs[!nbrs %in% selected.sum]
      if (length(newnbrs)!=0) {
        tb=table(newnbrs)
        cell.id=as.numeric(names(tb))
        prob=as.numeric(tb)/nb.length[cell.id]
        newnbr.freq=round((10^(4*prob-1)))
        newbr.seq=rep(cell.id, newnbr.freq)
        if (length(newbr.seq)>1) {
          selected[[p]] = c(selected[[p]], sample(newbr.seq, 1))
        } else { selected[[p]] = c(selected[[p]], newbr.seq)}
        n0=sum(sapply(selected, length))
      }
    }
  }
  length(unlist(selected))
  length(unique(unlist(selected)))
  return(selected)
}

# RandomRegionWindow ------------
#' Generate random connected region in a window
#'
#' This function generates random regions on a unit square.
#' @param nRegion nRegion is the No. of regions (e.g. nRegion=3)
#' @param nGrid nGrid is the No. of spots on x and y.
#' @param seed Random seed.
#' @return
#' \item{window:}{The window of each region.}
#' \item{area:}{Area proportion of the region in the window.}
#' @export
#'
RandomRegionWindow <- function(nRegion=3, nGrid=20, seed=NULL){
  if (nRegion==1) {
    win=vector(mode = "list", length = 1); win[[1]] = spatstat.geom::unit.square()
  } else {
    # generate polygon
    r <- raster::raster(ncols=nGrid, nrows=nGrid,xmn=0, xmx=1, ymn=0, ymx=1)

    # connect
    connected = connectUp(r=r, nRegion=nRegion, seed=seed)
    poly <- raster::rasterToPolygons(r)

    region=lapply(1:nRegion, function(f) terra::aggregate(poly[connected[[f]],]))
    coords=lapply(region, function(f) data.frame(f@polygons[[1]]@Polygons[[1]]@coords))
    win=lapply(1:length(coords), function(f)
      spatstat.geom::owin(poly=list(x=rev(coords[[f]])[,1],
                               y=rev(coords[[f]])[,2])))
    }

  # plot(win[[1]], xlim=c(0, 1), ylim=c(0,1), col=1)
  # plot(win[[2]], add=T, col=2)
  # plot(win[[3]], add=T, col=3)
  area=sapply(win, spatstat.geom::area.owin)
  return(list(window=win, area=area))
}

# get.n.vec.raw ------------
# Generate cell location in one region
# This function generates cell pools allowing selection due
# to cell overlaps, inhibitions and attractions.

get.n.vec.raw=function(n, cell.prop,
                       cell.inh.attr.input=NULL,
                       same.dis.cutoff=0.05) {
  n.vec.target=n.vec.use=round(n*cell.prop)
  # if exists cell-cell inhibition & attraction, inflate n.vec.use
  if (is.null(cell.inh.attr.input)==F) {
    for (i in 1:nrow(cell.inh.attr.input)) {
      if (cell.inh.attr.input[i,1] ==cell.inh.attr.input[i,2]) {
        n.vec.use[,cell.inh.attr.input[i,1]]=n.vec.use[,cell.inh.attr.input[i,1]]*
          (1+abs(cell.inh.attr.input[i,3]))
      }

      if (cell.inh.attr.input[i,1] !=cell.inh.attr.input[i,2]) {
        temp=ifelse(cell.inh.attr.input[i,3]>0, 1, 3) # more inflation for attraction
        n.vec.use[,cell.inh.attr.input[i,1]]=
          n.vec.use[,cell.inh.attr.input[i,1]]* (1+abs(cell.inh.attr.input[i,3]) * temp)
        n.vec.use[,cell.inh.attr.input[i,2]]=
          n.vec.use[,cell.inh.attr.input[i,2]]*(1+abs(cell.inh.attr.input[i,3]) * temp)
      }
    }
  }

  # continue inflating cell number to allow the deletion of same location
  n.vec.use.adj=round(n.vec.use^(1+same.dis.cutoff))
  return(list(n.vec.target=n.vec.target, n.vec.raw=n.vec.use.adj))
}

# cell.loc.1region.fc ------------
# Generate cell location in one region

cell.loc.1region.initiate <- function(n.inflation, window1,
                                     same.dis.cutoff = 0.05, seed) {

  set.seed(seed)
  # inflate cell numbers for deleting
  n.vec.raw=n.inflation$n.vec.raw;
  n.vec.target=n.inflation$n.vec.target
  n.vec.target.nonzero=n.vec.target
  K=length(n.vec.target)
  KP=colnames(n.vec.target)

  # random generation
  gen1=spatstat.random::rmpoint(n.vec.raw, win=window1, types=KP)

  # calculate cell-cell distance for deletion of same loc
  dis=spatstat.geom::pairdist(gen1)
  dis[lower.tri(dis, diag=T)]=1
  ratio= sqrt(spatstat.geom::area.owin(window1)/sum(n.vec.target))
  dis2=dis<same.dis.cutoff * ratio
  same.loc.idx=which(dis2 == T, arr.ind = TRUE)
  same.loc.delete.no=nrow(same.loc.idx)
  # get rid of the cells that are too close to each other (on same loc)
  mean.delete.prop=(n.vec.raw-n.vec.target)/n.vec.raw
  cell.mark=gen1$marks
  delete.prop1=mean.delete.prop[cell.mark[same.loc.idx[,1]]]
  delete.prop2=mean.delete.prop[cell.mark[same.loc.idx[,2]]]
  delete.fraction=delete.prop1/(delete.prop1+delete.prop2)
  delete.cell1=rbinom(same.loc.delete.no, 1, delete.fraction)
  delete.idx=c(same.loc.idx[delete.cell1==1,1],
               same.loc.idx[delete.cell1==0,2])

  gen2=gen1[setdiff(1:gen1$n, delete.idx), ]

  return(gen2)
}

# three factors that affect mu (logit prob changes from population average):
# (1) inhibitory/attraction cells of the same cell type
# (2) inhibitory/attraction cells of different cell type
# (3) even distribution
cell.loc.1region.refine <- function(pt.initial, n.inflation,
                                   cell.inh.attr.input1,
                                   even.distribution.coef = 0.1,
                                   grid.size.small = 19, grid.size.large = 45, seed) {
  set.seed(seed)

  KP=names(table(pt.initial$marks) )
  n2.vec=summary(pt.initial)$marks$frequency
  delete.idx=rep(0, pt.initial$n) # fill in this index to tell if each cell should be deleted or not
  mean.delete.prop=(n2.vec-n.vec.target)/n2.vec # average deletion prop.
  mean.delete.logit=mean.delete.prop/(1-mean.delete.prop) # logit
  # two resolutions
  r1 <- raster::raster(ncols=grid.size.small, nrows=grid.size.small,
               xmn=0, xmx=1, ymn=0, ymx=1)
  r2 <- raster::raster(ncols=grid.size.large, nrows=grid.size.large,
               xmn=0, xmx=1, ymn=0, ymx=1)
  # in total - for even distribution
  pt.den=raster::rasterize(cbind(pt.initial$x, pt.initial$y), r2, fun=function(x,...)length(x))
  value=raster::values(pt.den)
  value[is.na(value)]=0
  pixel.density=value/mean(value)
  pixel.idx=raster::cellFromXY(pt.den, cbind(pt.initial$x, pt.initial$y))
  den.total=pixel.density[pixel.idx]
  mu1=den.total*even.distribution.coef


  # by cell  - for inhibitory/attraction cells
  if (is.null(cell.inh.attr.input1)) {
    cell.inh.attr.input1=cbind(Cell1=1, Cell2=1, Strength= 0)
    }
  # local density by cell type 1
  pt.cell.den1=lapply(KP, function(k)
    raster::rasterize(cbind(pt.initial[which(pt.initial$marks==k),]$x,
                            pt.initial[which(pt.initial$marks==k),]$y),
              r1, fun=function(x,...) length(x)))

  value.cell.r1=lapply(1:K, function(f) raster::values(pt.cell.den1[[f]]))
  for (k in 1:K) {value.cell.r1[[k]][is.na(value.cell.r1[[k]])]=0}
  pixel.density.cell1=lapply(1:K, function(f)
    value.cell.r1[[f]]/mean(value.cell.r1[[f]]))

  # local density by cell type 2
  pt.cell.den2=lapply(KP, function(k) raster::rasterize(cbind(pt.initial[which(pt.initial$marks==k),]$x,
                                                              pt.initial[which(pt.initial$marks==k),]$y), r2, fun=function(x,...)length(x)))
  value.cell.r2=lapply(1:K, function(f) raster::values(pt.cell.den2[[f]]))

  for (k in 1:K) {value.cell.r2[[k]][is.na(value.cell.r2[[k]])]=0}
  pixel.density.cell2=lapply(1:K, function(f) value.cell.r2[[f]]/mean(value.cell.r2[[f]]))
  names(pt.cell.den1)=names(pixel.density.cell1)=
    names(pt.cell.den2)=names(pixel.density.cell2)=KP

  # by cell type
  pt.initial.cell=split(pt.initial)
  for (k in which(mean.delete.prop>0)) { # these are the cell types that need deletion; do 1 by 1

    pt.initial.cell.loc=cbind(pt.initial.cell[[k]]$x, pt.initial.cell[[k]]$y)
    n0.k=nrow(pt.initial.cell.loc)

    # find cell types that impact the deletion probability of cell type k
    row.idx= apply(cell.inh.attr.input1, 1, function(f) KP[k] %in%f[1:2] ) # any
    # same cell type
    row.same.idx= apply(cell.inh.attr.input1, 1, function(f) KP[k] %in%f[1] & KP[k] %in%f[2])
    # different cell types
    row.diff.idx= row.idx ==T &  row.same.idx ==F
    # calculate mu2 the contribution of same cell inhibition/attraction
    if (sum(row.same.idx)==0) {mu2=matrix(rep(0, n0.k))} else {
       pixel.idx.cell1=raster::cellFromXY(pt.cell.den1[[k]], pt.initial.cell.loc)
       den.same1=pixel.density.cell1[[k]][pixel.idx.cell1]

       pixel.idx.cell2=raster::cellFromXY(pt.cell.den2[[k]], pt.initial.cell.loc)
       den.same2=pixel.density.cell2[[k]][pixel.idx.cell2]

       den.same=(den.same1+den.same2)/2
       mu2=den.same* cell.inh.attr.input1[row.same.idx, 3]
     }

     # calculate mu3 the contribution of different cell inhibition/attraction
     if (sum(row.diff.idx)==0) {mu3=matrix(rep(0, n0.k))} else {
       diff.cell.density=NULL
       pair.idx=unlist(apply(as.matrix(cell.inh.attr.input1[row.diff.idx, 1:2]),1,
                             function(f) setdiff(f, KP[k])))
       for (pk in pair.idx) {
          pixel.idx.cell1=raster::cellFromXY(pt.cell.den1[[pk]], pt.initial.cell.loc)
          den.diff1=pixel.density.cell1[[pk]][pixel.idx.cell1]

          pixel.idx.cell2=raster::cellFromXY(pt.cell.den2[[pk]], pt.initial.cell.loc)
          den.diff2=pixel.density.cell2[[pk]][pixel.idx.cell2]

          diff.cell.density=cbind(diff.cell.density,(den.diff1+den.diff2)/2)
       }
       mu3=diff.cell.density%*% matrix(cell.inh.attr.input1[row.diff.idx, 3])
     }

    mu=  mu1[which(pt.initial$marks==KP[k])]-mu2-mu3

    est.mean.delete.prop=1
    c=10
    while (est.mean.delete.prop>mean.delete.prop[k]) {
      c=c-0.1
      mu.logit=mean.delete.logit[[k]] +scale(mu, center=T, scale=F) +c
      delete.prop=exp(mu.logit)/(1+exp(mu.logit))
      est.mean.delete.prop=mean(delete.prop)
    }

    # plot(den.same, delete.prop)
    # plot(diff.cell.density, delete.prop)

    delete.idx.k=  rbinom(n0.k, 1, delete.prop)

    delete.idx[pt.initial$marks==KP[k]]=delete.idx.k
  }

  gen3=pt.initial[setdiff(1:pt.initial$n, which(delete.idx==1)), ]

  return(final.ppp=gen3)
}



# cell.loc.fc ------------
#' Generate cell location data from parameters.
#'
#' This function generates cell locations for all regions.
#' @param N No. of cells to generate.
#' @param win Spatial window within which cells will be generated.
#' @param cell.prop Proportion of cell in each cell type to generated for each region.
#' @param cell.inh.attr.input (Default=NULL). A matrix providing cell-cell location
#' attraction and inhibition parameters. Example input is like: cell.inh.attr.input=
#' cbind(Cell1=c("A", "B", "C"), Cell2=c("B","D","E), Strength= c(2, 0, -2)).
#' @param same.dis.cutoff (Default = 0) Cells with distance less than this cutoff
#' will be considered as overlapping cells non-realistic in real ST data,
#'  and only one will be kept.
#' @param even.distribution.coef (Default =0). The higher the value, the more we
#' require evenly distributed locations of cells, as opposed to randomly generated
#' locations (Poisson process).
#' @param grid.size.small (Default =19). If some levels of even distributions of cells are imposed,
#' this parameter is needed to cut the simulation window into grids of different
#' sizes to smooth cell densities.
#' @param grid.size.large (Default =45). If some levels of even distributions of cells are imposed,
#' this parameter is needed to cut the simulation window into grids of different
#' sizes to smooth cell densities.
#' @param seed Random seed
#' @return Cell location
#' @export


cell.loc.fc=function(N, win, cell.prop, cell.inh.attr.input=NULL,
                     same.dis.cutoff =0,
                     even.distribution.coef=0,
                     grid.size.small=19, grid.size.large=45,
                     seed) {

  R=length(cell.prop)
  cell.loc=vector("list", R);

  for(r in 1:R) {
    # inflate cell numbers for deleting
    n.inflation=get.n.vec.raw(n=round(win$area[r]*N),
                              cell.prop=cell.prop[[r]],
                              cell.inh.attr.input=cell.inh.attr.input[[r]],
                              same.dis.cutoff =same.dis.cutoff)
    # Simulate initial cell locations for an inflated number of cells
    pt1=cell.loc.1region.initiate(n.inflation=n.inflation,
                              window1=win$window[[r]],
                              same.dis.cutoff = same.dis.cutoff, seed=seed*11+r*141)

   # Delete the inflation from cell-cell interactions
    cell.loc[[r]]=cell.loc.1region.refine(pt.initial=pt1, n.inflation=n.inflation,
                                        cell.inh.attr.input1=cell.inh.attr.input[[r]],
                                        even.distribution.coef=even.distribution.coef,
                                        grid.size.small=grid.size.small, grid.size.large=grid.size.large,
                                        seed=seed*11+r*141)
  }
  names(cell.loc)=paste0("Region", 1:R)

  return(cell.loc)
}
