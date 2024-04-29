#' Assigns cells into spatial spots for multi-cell resolution data
#'
#' This function assigns cells into spatial spots. Each spot may contain zero, one, or
#' multiple cells. Spots with zero cells won't be in the output.
#' @param expr expression profile of the cells
#' @param cell_feature cell features like their spatial coordinates.
#' @param NoSpot Number of targeted spots.
#' @return
#' \item{count:}{Expression profile of the spots.}
#' \item{spot_feature:}{Spot feature such as spot center coordinates and cell allocations
#' within the spots .}
#' @import dplyr
#' @export

multicell=function(expr, cell_feature, NoSpot=500) {
  cell_loc=cell_feature[,c("x.loc", "y.loc")]
  xrange=range(cell_loc[,1])
  yrange=range(cell_loc[,2])
  m=round(sqrt(NoSpot))

  r <- raster::raster(ncols=m, nrows=m,xmn=xrange[1],
                      xmx=xrange[2], ymn=yrange[1], ymx=yrange[2])
  spot.idx=raster::cellFromXY(r, cbind(cell_loc[,1], cell_loc[,2]))
  #
  # deal with spots with no cell
  expr2=sapply(1:nrow(expr), function(f1) sapply(
    1: (m^2), function(f2) sum(as.numeric(expr[f1, spot.idx==f2]), na.rm=T)))
  expr3=t(expr2)
  rownames(expr3)=rownames(expr)
  colnames(expr3)=paste0("Spot", 1:ncol(expr3))

  # Spot Coordinates
  spot_col=raster::colFromCell(r, 1:nrow(expr2))
  spot_row=raster::rowFromCell(r, 1:nrow(expr2))
  spot_coordinates=raster::xyFromCell(r, 1:nrow(expr2))
  spot_loc=data.frame(Spot=colnames(expr3), col=spot_col,
                      row=spot_row, spot_coordinates)


  # Spot Region
  if (is.null(cell_feature$region)==F) {
    dat1=data.frame(spot.idx, region=cell_feature$region) %>%
      group_by(spot.idx, region) %>%
      mutate(count=1) %>%
      summarise(abundance = sum(count))
    # add zero
    dat2=data.frame(spot.idx=setdiff(1: (m^2), spot.idx),
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
  if (is.null(cell_feature$annotation)==F) {
    dat1=data.frame(spot.idx, cell_feature) %>%
      dplyr::select(spot.idx, annotation) %>%
      group_by(spot.idx, annotation) %>%
      mutate(count=1) %>%
      summarise(abundance = sum(count))
    # add zero
    dat2=data.frame(spot.idx=setdiff(1: (m^2), spot.idx),
                    annotation=dat1$annotation[1],
                    abundance=0)
    dat=rbind(dat1, dat2) %>%
      tidyr::pivot_wider(names_from = annotation, values_from = abundance,
                         values_fill = 0)
    spot_loc=data.frame(spot_loc, dat[,-1])
  }

  return(list(count=expr3, spot_feature=spot_loc))
}
