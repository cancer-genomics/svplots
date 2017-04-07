#' @include AllGenerics.R
NULL

#' threshold values of a numeric vector
#'
#' Threshold the values of a numeric vector according to
#' user-specified limits.  Values exceeding the threshold are reset to
#' the threshold and jittered by an amount specified by the user to
#' reduce overplotting.
#'
#' @param x numeric vector
#' @param lim numeric vector of length 2 indicating the values at which to threshold \code{x}
#' @param amount scalar indicating how much to jitter the points at the threshold.
#' @seealso \code{\link{jitter}}
#' @examples
#' x <- rnorm(10)
#' threshold(x, c(-0.5,0.5))
#' @export
threshold <- function(x, lim=c(-Inf,Inf), amount=0){
  notna <- !is.na(x)
  index1 <- which(x <= lim[1] & notna)
  index2 <- which(x >= lim[2] & notna)
  x1 <- runif(length(index1), lim[1], lim[1]+amount)
  x2 <- runif(length(index2), lim[2]-amount, lim[2])
  x[index1] <- x1
  x[index2] <- x2
  return(x)
}

#' Create a GRanges object with bin-level log ratios for visualizing
#' low-level copy number data.
#'
#' @seealso \code{\link[ggbio]{plotGrandLinear}} for complete details
#'   and \code{plotGenomeLogRatios} for the wrapper function
#'   implemented in this package
#'
#' @examples
#' library(svdata)
#' data(pviews)
#' paths(pviews) <- file.path(system.file("extdata", package="svdata"),
#'                            rdsId(pviews))
#' gr <- genomeLogRatios(pviews[, 1])
#' gr
#' 
#' 
#' @export
#' @param view A \code{PreprocessViews2} object 
#' @param select_gr A \code{GRanges} object indicating which intervals to plot
#' @param lim A length-two numeric vector for y-axis limit
#' 
#' @param thin integer indicating how much to thin the data.  A value
#'   of 10 means that every 10th log ratio is returned
genomeLogRatios <- function(view, select_gr, lim=c(-4, 3), thin=10){
  i <- seq(1, nrow(view), thin)
  view <- view[i, ]
  rr <- rowRanges(view)
  rr$log_ratio <- assays(view)[, 1]
  if(!missing(select_gr)){
    rr$select <- overlapsAny(rr, select_gr)
    rr$nselect <- countOverlaps(rr, select_gr)
  }
  rr$log_ratio <- threshold(rr$log_ratio, lim)
  ## so that select values are plotted last
  if(!missing(select_gr)) rr <- rr[order(rr$select)]
  rr
}

#' Visualize log ratios along the genome with plotGrandLinear
#'
#' This is a simple wrapper function for \code{plotGrandLinear}
#' implemented in the \code{ggbio} package.
#'
#' @seealso See \code{\link[ggbio]{plotGrandLinear}} for better
#'   control of plotting genomic data.  See
#'   \code{\link{genomeLogRatios}} for creating a \code{GRanges}
#'   object containing log ratios of the bin-level copy number
#'   estimates.
#'
#' @examples
#'  library(svdata)
#'  data(pviews)
#'  paths(pviews) <- file.path(system.file("extdata", package="svdata"),
#'                             rdsId(pviews))
#'  gr <- genomeLogRatios(pviews[, 1], thin=50)
#'  ##
#'  ## plotGrandLinear currentl ignores size and color
#'  ##
#'  p <- plotGenomeLogRatios(gr, title=colnames(pviews)[1], size=0.1)
#'  print(p)
#' 
#' @export
#' 
#' @param gr A \code{GRanges} object containing 'log_ratio' in
#'   the meta-columns.
#' @param title a character string providing a title for the figure
#' @param ... additional arguments passed to \code{plotGrandLInear}
plotGenomeLogRatios <- function(gr, title, ...){
  log_ratio <- NULL
  p <- plotGrandLinear(gr, geom = "point", coord = "genome",
                       aes(y = log_ratio), ...)

  if(!missing(title)){
    p <- p + ggtitle(title)
  }
  p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
  p
}
