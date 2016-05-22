#' Expand an existing class of patches.
#
#' @inheritParams makePatch
#' @inheritParams makeClass
#' @param class (or patch),
#' @param size integer. Size of expansion, as number of raster cells.
#' @param bgr integer. The background available where expansion is allowed (i.e. shrinking classes).
#' @return A RasterLayer object. If \code{rast=FALSE} returns a list of vectors, each containing the \code{context} cells assigned to each patch.
#' @examples
#' require(raster)
#'
#' m = matrix(0, 33, 33)
#' r = raster(m, xmn=0, xmx=10, ymn=0, ymx=10)
#' r = makeClass(r, 5, 10)
#' plot(r)
#'
#' rr = expandClass(r, 1, 100)
#' plot(rr)
#'
#' ## This function can be used to mimic shapes, by providing a skeleton:
#' m[,17] = 1
#' r = raster(m, xmn=0, xmx=10, ymn=0, ymx=10)
#' plot(r)
#'
#' rr = expandClass(r, 1, 100)
#' plot(rr)
#' @export
expandClass <- function(context, class, size, bgr=0, pts = NULL) {
  bd <- raster::boundaries(context, type='outer', classes=TRUE, directions=8)
  bd <- t(raster::as.matrix(bd))
  if(!is.matrix(context)) {m <- t(raster::as.matrix(context))} else {m <- context}
  if(is.null(pts)){edg <- which(bd==1 & m==class)} else {edg <- pts}
  bgr <- which(.is.elem(m, bgr))
  if(length(bgr) == 0){stop('No cells available, landscape full')}
  if(size > (length(bgr))){stop('Expansion size bigger than available landscape')}
  pts <- ifelse(length(edg) == 1, edg, sample(edg, 1) )
  cg <- pts
  dim1 <- dim(m)[1]
  dim2 <- dim(m)[2]
  while(length(cg) < size){
    if(pts %% dim1 == 0){
      rr <- dim1
      cc <- pts/dim1
    } else {
      cc <- trunc(pts/dim1) + 1
      rr <- pts - (cc-1) * dim1
    }
    ad <- c(rr-1, rr+1, rr, rr, cc, cc, cc-1, cc+1)
    ad[ad <= 0 | c(ad[1:4] > dim1, ad[5:8] > dim2)] <- NA
    ad <- ad[1:4] + (ad[5:8]-1)*dim1
    ad <- ad[!is.na(ad)]
    ad <- ad[which(!.is.elem(ad, c(cg, edg)))]
    ## The following stands for {ad <- bgr[which(bgr %in% ad)]}
    ##d <- vapply(ad, function(x) .fastSearch(as.numeric(bgr), x, "=="), 1); ad <- bgr[d]
    d <- fastmatch::fmatch(ad, bgr, nomatch = 0)
    ad <- bgr[d]
    if(length(ad) == 0) {
      edg <- edg[edg!=pts]
      if(length(edg) <= 1) {
        if(length(cg) == 1){
          warning('Expanding classes do not touch shrinking classes. Input raster returned')
          break
        } else {
          warning('Maximum patch size reached: ', length(cg))
          break
        }
      }
      pts <- sample(edg, 1)
      next
    }
    cg <- c(cg, ad)
    edg <- c(edg[which(edg != pts)], ad)
    pts <- ifelse(length(edg) == 1, edg, sample(edg, 1) )
  }
  context[cg] <- class
  return(context)
}
