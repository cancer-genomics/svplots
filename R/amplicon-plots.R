#' @include help.R
NULL

#' Plot amplicons as nodes and the links between amplicons as edges
#'
#' Linked amplicons are represented as a \code{AmpliconGraph}. With
#' amplicons as nodes and edges representing links between amplicons
#' from paired reads, a \code{graphNEL} object encapsulates this data.
#' Plotting utilities in the \code{Rgraphviz} package can be used to
#' visualize these graphs.
#'
#' @export
#'
#' @seealso See \code{\link[svcnvs]{sv_amplicon_exp}} for constructing
#'   an \code{AmpliconGraph} object. See
#'   \code{\link{qualitativeColors}} for a list of qualitative
#'   colors. See \code{\link[graph]{graphNEL-class}} for details on
#'   graph representation.
#'
#' @param ag a \code{AmpliconGraph} object
#' @param col_list a list of color palettes
plot_amplicons <- function(ag, col_list=qualitativeColors()){
  if(length(ag)==0){
    message("No amplicons -- nothing to plot")
    return(invisible())
  }
  dark_colors <- c("#332288", "#661100", "#882255", "#AA4499")
  ##ar <- ampliconRanges(ag)
  ar <- nodes(ag)
  ##ar <- keepSeqlevels(ar, unique(chromosome(ar)))
  chroms <- sapply(strsplit(ar, ":"), "[", 1)
  sl <- unique(chroms)
  ##sl <- seqlevels(ar)
  L <- length(sl)
  L <- ifelse(L > 12, 12, L)
  ##sl <- factor(chromosome(ar), levels=sl)
  sl <- factor(chroms, levels=sl)
  color_nodes <- col_list[[L]][as.integer(sl)]
  g1 <- graph(ag)
  color_nodes <- setNames(color_nodes, nodes(g1))
  nodenames <- setNames(nodes(g1), nodes(g1))
  text_col <- setNames(rep("black", length(nodes(g1))), nodes(g1))
  text_col[color_nodes %in% dark_colors] <- "gray90"
  nodeRenderInfo(g1) <- list(label=nodenames, fill=color_nodes, textCol=text_col)
  ##nodeAttrs <- list(fillcolor=colors)
  nodeAttrs <- list(fillcolor=color_nodes)
  attrs <- list(node=list(shape="rectangle",
                          fixedsize=FALSE),
                graph=list(rankdir="LR"))
  graph_object <- layoutGraph(g1,
                              attrs=attrs,
                              nodeAttrs=nodeAttrs)
  if(numEdges(graph_object) > 0){
    renderGraph(graph_object)
  } else plot(graph_object)
  graph_object
}
