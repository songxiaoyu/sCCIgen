#' Assigns cells into spatial spots for multi-cell resolution data
#'
#' This function assigns cells into spatial spots. Each spot may contain zero, one, or
#' multiple cells. Spots with zero cells won't be in the output.
#' @param expr expression profile of the cells
#' @param spatial spatial features like their spatial X and Y coordinates.
#' @param NoSpot Number of targeted spots.
#' @param cl No. of cores.
#' @return
#' \item{count:}{Expression profile of the spots.}
#' \item{spot_feature:}{Spot feature such as spot center coordinates and cell allocations
#' within the spots .}
#' @import dplyr
#' @export

multicell=function(expr, spatial, NoSpot=500, cl=1) {
  cell_loc=spatial[,c("x.loc", "y.loc")]
  xrange=range(cell_loc[,1])
  yrange=range(cell_loc[,2])
  m=round(sqrt(NoSpot))
  mm=m^2

  r <- raster::raster(ncols=m, nrows=m,xmn=xrange[1],
                      xmx=xrange[2], ymn=yrange[1], ymx=yrange[2])
  spot.idx=raster::cellFromXY(r, cbind(cell_loc[,1], cell_loc[,2]))

  # Spot expression
  if (cl==1) {
    expr2=matrix(0, ncol=mm, nrow=nrow(expr))
    for (i in 1: mm) {
      expr2[,i]=apply(expr, 1, function(f) sum(as.numeric(f)[spot.idx==i], na.rm=T))
    }
  } else {
    expr2=foreach (i = 1: mm, .combine="cbind") %dopar% {
      apply(expr, 1, function(f) sum(as.numeric(f)[spot.idx==i], na.rm=T))
    }
  }


  rownames(expr2)=rownames(expr)
  colnames(expr2)=paste0("Spot", 1:mm)

  # Spot Coordinates
  spot_col=raster::colFromCell(r, 1:mm)
  spot_row=raster::rowFromCell(r, 1:mm)
  spot_coordinates=raster::xyFromCell(r, 1:mm) %>% round(., digits = 3)
  spot_loc=data.frame(Spot=colnames(expr2), col=spot_col,
                      row=spot_row, spot_coordinates)


  # Spot Region
  if (is.null(spatial$region)==F) {
    dat1=data.frame(spot.idx, region=spatial$region) %>%
      group_by(spot.idx, region) %>%
      mutate(count=1) %>%
      summarise(abundance = sum(count))
    # add zero
    dat2=data.frame(spot.idx=setdiff(1:mm, spot.idx),
                    region=dat1$region[1],
                    abundance=0)
    dat=rbind(dat1, dat2) %>%
      tidyr::pivot_wider(names_from = region, values_from = abundance,
                         values_fill = 0) %>%
      ungroup() %>%
      dplyr::select(-spot.idx) %>%
      rowwise() %>%
      mutate(region = names(.)[which.max(c_across(everything()))])

    spot_loc=data.frame(spot_loc, region=dat$region)
  }

  # Spot's Cell Type Count
  if (is.null(spatial$annotation)==F) {

    dat=data.frame(spot.idx, spatial) %>%
      group_by(spot.idx, annotation) %>%
      mutate(count=1) %>%
      summarise(abundance = sum(count)) %>%
      tidyr::pivot_wider(names_from = annotation,
                         values_from = abundance,
                         values_fill = 0)

    spot_loc=spot_loc%>% rownames_to_column("spot.idx")%>%
      left_join(dat %>% mutate(spot.idx=as.character(spot.idx)))
  }

  return(list(count=expr2, spot_feature=spot_loc))
}
