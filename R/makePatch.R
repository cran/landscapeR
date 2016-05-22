#' Create a single patch
#'
#' @param context Raster object or matrix, an empty landscape raster or a mask indicating where the patch cannot be generated (see bgr below).
#' @param size integer. Size of the patch to be generated, as number of raster cells.
#' @param spt integer. The seed point location around which the patch is generated (a random point is given by default)
#' @param bgr integer. Value of background cells, where a patch can be generated (default is zero). Cells/classes which cannot be changed must have a different value.
#' @param edge logical. Should the vector of edge cells of the patch be returned?
#' @param rast logical. If TRUE returns a Raster object, otherwise a vector of cell numbers where the patch occurs
#' @param val integer. The value to be assigned to patch cells, when \code{rast=TRUE}
#' @return A vector of matrix cell numbers, or a RasterLayer object if \code{rast=TRUE}.
#' @details The patch is created starting from the seed point and iteratively sampling randomly neighbouring cells at the edge of the patch.
#' There is a tolerance of +/- 3 cells from the patch size declared in \code{size} argument.
#' @examples
#' library(raster)
#' m = matrix(0, 33, 33)
#' r = raster(m, xmn=0, xmx=10, ymn=0, ymx=10)
#' patchSize = 500
#' rr = makePatch(r, patchSize, rast=TRUE)
#' plot(rr)
#'
#' ## Create a patch with value 3, starting from the centre cell
#' m = matrix(0, 33, 33)
#' r = raster(m, xmn=0, xmx=10, ymn=0, ymx=10)
#' newVal = 3
#' centre = 545
#' cells = makePatch(r, patchSize, centre)
#' r[cells] = newVal
#' plot(r)
#'
#' ## Now create a new patch with value 10 and size 100 inside the existing patch
#' rr = makePatch(r, 100, bgr=newVal, rast=TRUE, val=10)
#' plot(rr)
#' @export
makePatch <- function(context, size, spt = sample.int(length(context), 1), bgr=0, edge=FALSE, rast=FALSE, val=1) {
  if(!is.matrix(context)) {
    m <- t(raster::as.matrix(context))
  } else {
    m <- context
  }
  bgr <- which(.is.elem(m, bgr))
  if(length(bgr) == 0){
    stop('No cells available, landscape full')
  }
  if(size > (length(bgr))){
    warning('Patch size bigger than available landscape')
  }
  if(!spt %in% bgr){
    warning('Seed point within mask. Started from another sampled outside')
    spt <- ifelse(length(bgr) > 1, sample(bgr, 1), bgr)
  }
  cg <- edg <- spt
  dim1 <- dim(m)[1]
  dim2 <- dim(m)[2]
  while(length(cg) < size){
    if(spt %% dim1 == 0){
      rr <- dim1
      cc <- spt/dim1
    } else {
      cc <- trunc(spt/dim1) + 1
      rr <- spt - (cc-1) * dim1
    }
    ad <- c(rr-1, rr+1, rr, rr, cc, cc, cc-1, cc+1)
    ad[ad <= 0 | c(ad[1:4] > dim1, ad[5:8] > dim2)] <- NA
    ad <- ad[1:4] + (ad[5:8]-1)*dim1
    ad <- ad[!is.na(ad)]
    ad <- ad[which(!.is.elem(ad, c(cg, edg)))]
    ## The following stands for {ad <- bgr[which(bgr %in% ad)]}
    d <- fastmatch::fmatch(ad, bgr, nomatch = 0)
    ad <- bgr[d]
    if(length(ad) == 0) {
      edg <- edg[edg!=spt]
      if(length(edg) <= 1) {
        warning('Maximum patch size reached: ', length(cg))
        break
      }
      spt <- sample(edg, 1)
      next
    }
    cg <- c(cg, ad)
    edg <- c(edg[which(edg != spt)], ad)
    spt <- ifelse(length(edg) == 1, edg, sample(edg, 1) )
  }
  if(rast == TRUE) {
    context[cg]<-val
    return(context)
  } else if (edge == TRUE) {
    return(list(cg, edg))
  } else {
    return(cg)
  }
}

##
.is.elem <- function(el, set) {
  fastmatch::fmatch(el, set, 0L) > 0L
}
