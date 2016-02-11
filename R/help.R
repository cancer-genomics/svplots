#' Tools for visualizing structural variant sequence data
#'
#' Functions and methods for visualizing structural variant sequence data
#'
#' @docType package
#' @name svplots
#' @import methods
#' @import svclasses
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import grid
#' @import ggplot2
#' @importFrom lattice lsegments current.panel.limits xyplot panel.abline panel.xyplot
#' @importFrom ggbio plotGrandLinear
#' @import graph
#' @importMethodsFrom S4Vectors elementLengths
#' @importMethodsFrom GenomicAlignments first last
#' @importFrom IRanges IRanges
#' @importFrom Rgraphviz layoutGraph plot
#' @importMethodsFrom Rgraphviz renderGraph
#' @importFrom RColorBrewer brewer.pal
#' @importFrom biovizBase getIdeoGR
NULL