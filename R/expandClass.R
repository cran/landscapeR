#' Expand an existing class of patches.
#
#' @inheritParams makePatch
#' @inheritParams makeClass
#' @param class The raster value of class (or patch) to expand.
#' @param size integer. Size of expansion, as number of raster cells.
#' @param bgr integer. The background available where expansion is allowed (i.e. shrinking classes).
#' @return A SpatRaster object. If \code{rast=FALSE} returns a list of vectors, each containing the \code{context} cells assigned to each patch.
#' @examples
#' library(terra)
#'
#' m = matrix(0, 33, 33)
#' r = rast(m)
#' ext(r) = c(0, 10, 0, 10)
#' r = makeClass(r, 5, 10)
#' plot(r)
#'
#' rr = expandClass(r, 1, 100)
#' plot(rr)
#'
#' ## This function can be used to mimic shapes, by providing a skeleton:
#' m[,17] = 1
#' r = rast(m)
#' ext(r) = c(0, 10, 0, 10)
#' plot(r)
#'
#' rr = expandClass(r, 1, 100)
#' plot(rr)
#' @export
expandClass <- function(context, class, size, bgr=0, pts = NULL) {
  if(length(class) > 1){ stop('A single value only is admitted for "class" argument.') }
  if(class %in% bgr){ warning('Value to attribute to patches same as background cells value (arg. "class" equals "bgr").') }
  if(any(is.na(size) | size <=0)){ stop('Invalid "size" argument provided.') }
  bd <- terra::boundaries(context, inner=FALSE, classes=TRUE, directions=8)
  #----- LEM: the as.numeric is required for type matching in Rcpp -----------#
  bd <- terra::as.matrix(bd, wide=T)
  bd <- t(matrix(as.numeric(bd), ncol=ncol(bd), nrow=nrow(bd)))
  if(!is.matrix(context)) {
    #----- LEM: the as.numeric is required for type matching in Rcpp ---------#
    mtx <- terra::as.matrix(context, wide=T)
    mtx <- t(matrix(as.numeric(mtx), ncol=ncol(mtx), nrow=nrow(mtx)))
  } else {
    mtx <- context
  }
  if(is.null(pts)){
    edg <- which(bd==1 & mtx==class)
  } else {
    edg <- pts
  }

  if(length(bgr) > 1){
    p_bgr <- which(is.element(mtx, bgr[-1]))
    vals <- mtx[p_bgr]
    bgr <- bgr[1]
    bgrCells <- c(which(mtx == bgr), p_bgr)
    .assignValues(bgr, p_bgr, mtx) # mtx[p_bgr] <- bgr
  } else {
    bgrCells <- which(mtx == bgr)
  }
  if(length(bgrCells) == 0){stop('No cells available, landscape full')}
  if(size > (length(bgrCells))){stop('Expansion size bigger than available landscape')}
  if(length(edg) == 1){
    pts <- edg
  } else {
    pts <- sample(edg, 1)
  }
  cg <- 1
  while(cg < size){
    ad <- .contigCells(pts, bgr, mtx)
    if(length(ad) == 0) {
      edg <- edg[edg != pts]
      if(length(edg) <= 1) {
        if(cg == 1){
          warning('Expanding classes do not touch shrinking classes. Input raster returned.')
          break
        } else {
          warning('Maximum patch size reached: ', cg)
          break
        }
      }
      pts <- sample(edg, 1)
      next
    }
    .assignValues(class, ad, mtx) # mtx[ad] <- class
    cg <- cg + length(ad)
    edg <- c(edg[edg != pts], ad)
    if(length(edg) == 1){
      pts <- edg
    } else {
      pts <- sample(edg, 1)
    }
  }
  if(exists('p_bgr')){
#    id <- mtx[p_bgr] != class
    id <- .subset(mtx, p_bgr) != class
#    mtx[p_bgr[id]] <- vals[id]
    mtx[.subset(p_bgr, id)] <- .subset(vals, id)
  }
  terra::values(context) <- t(mtx)
  return(context)
}
