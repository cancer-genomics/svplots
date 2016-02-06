cnTrack <- function(dirs, id, nt=200){
  cbs_file <- file.path(dirs[["0cbs"]], paste0(id, ".rds"))
  gr.cn <- readRDS(cbs_file)
  cn <- gr.cn$seg.mean
  cn[cn < -2.5] <- -2.5
  cn[cn > 4] <- 4
  ## multiply by -1 so that higher copy number goes towards the
  ## center
  cn <- cn * -1
  gr <- granges(gr.cn)
  gr$cn <- cn
  names(gr.cn) <- NULL
  ## microdeletions/duplications
  ##w <- rep(1, length(gr.cn))
  ##w[width(gr.cn) < 20e3] <- 3
  ##gr.cn$lwd <- w
  gr.cn <- sort(gr.cn)
  return(gr)
}

#' Creates three tracks for plotting segment means, rearrangements,
#' and CNVs
#'
#' @return a named list of \code{GRanges} objects that serve as tracks
#'   for making circos plots.
#' 
#' @export
#' @param id a length-one character vector of the sample identifier
#' @param dirs a \code{DataPaths} object
#' @param slstyle one of the available \code{seqlevelStyle}
#' 
#' @param MINSEP the minimum size of a rearrangement ('link'
#'   geometry), where size is measured as the distance between two
#'   intervals.  This parameter is ignored for inter-chromosomal
#'   rearrangements.
#' 
#' @seealso \code{\link[GenomeInfoDb]{seqlevelsStyle}}
circosTracks <- function(id, dirs, slstyle="NCBI", MINSEP=50e3){
  id.rds <- paste0(id, ".rds")
  svs <- readRDS(file.path(dirs["deletions"], id.rds))
  d <- granges(variant(svs))
  d$type <- rep("del", length(d))
  amps <- readRDS(file.path(dirs["amplicons"], id.rds))
  a <- granges(amplicons(amps))
  a$type <- rep("amp", length(a))
  seqlevelsStyle(d) <- slstyle
  seqlevelsStyle(a) <- slstyle
  names(d) <- names(a) <- NULL

  ## ideogram
  hg <- getIdeoGR(d)
  names(hg) <- NULL
  rlist <- readRDS(file.path(dirs["rearrangements/3blat_unmapped"], id.rds))
  r <- linkedBins(rlist)
  if(length(r) > 0){
    seqlevelsStyle(r) <- slstyle
    seqlevelsStyle(r$linked.to) <- slstyle
    r <- keepSeqlevels(r, seqlevels(d))
    seqinfo(r) <- seqinfo(hg)
    lt <- r$linked.to
    seqlevels(lt, force=TRUE) <- seqlevels(r)
    seqinfo(lt) <- seqinfo(r)
    r$linked.to <- lt
    names(r) <- NULL
  } else {
    r$linked.to <- GRanges()
  }
  gr.cn <- cnTrack(dirs, id, 200) ## number of cn estimates per chrom.
  seqlevelsStyle(gr.cn) <- slstyle

  hg <- keepSeqlevels(hg, seqlevels(d))
  gr.cn <- keepSeqlevels(gr.cn, seqlevels(d))
  seqinfo(gr.cn) <- seqinfo(hg)
  if(length(r) > 0){
    types <- c("intra-chrom", "inter-chrom")
    r$rearrangement <- types[as.integer(chromosome(r)!=chromosome(r$linked.to))+1]
    r$rearrangement <- factor(r$rearrangement, levels=c("intra-chrom", "inter-chrom"))
    rintra <- r[r$rearrangement=="intra-chrom"]
    rinter <- r[r$rearrangement=="inter-chrom"]
    sep <- abs(start(rintra$linked.to) - end(rintra)) > MINSEP
    rintra <- rintra[sep]
    r <- c(rintra, rinter)
  }
  cnvs <- c(d, a)
  cnvs$type <- factor(cnvs$type, levels=c("del", "amp"))
  seqinfo(cnvs) <- seqinfo(hg)
  list(cnvs=cnvs, r=r, hg=hg, gr.cn=gr.cn, id=id)
}

#' Create a circos plot with tracks for rearrangements, copy number
#' variants, and segment means.
#'
#' This function is a wrapper for the \code{circle} and \code{ggbio}
#' functions defined the \code{ggbio} package, creating circos plots
#' emphasizing structural variants including rearrangements and DNA
#' copy number alterations.
#'
#' @seealso See \code{\link{ggbio}}
#' 
#' @export
#' 
#' @param tracks a named \code{list} of \code{GRanges} objects as
#'   provided by \code{circosTracks}
#' 
#' @param cbcolors 
circosPlot <- function(tracks, cbcolors){
  r <- tracks[["r"]]
  seqinfo(r) <- seqinfo(tracks[["hg"]])
  seqinfo(r$linked.to) <- seqinfo(r)
  p <- ggbio(buffer=0) + circle(tracks[["hg"]], geom = "text", aes(label = seqnames),
                                    vjust = 0, size = 3, radius=38) +
    circle(r, geom = "link", linked.to = "linked.to",
           radius=25, color=cbcolors[2]) +
    circle(tracks[["cnvs"]], geom = "rect", aes(fill=type, color=type),
           trackWidth=3, radius=28) +  
    circle(tracks[["gr.cn"]], geom = "segment",
           color="black",
           fill="black",
           aes(y=cn),
           grid = FALSE, size=0.5, radius=31, trackwidth=3) +
    circle(tracks[["hg"]], geom = "ideo", fill = "beige",
           color = "gray",
           radius=35,
           trackwidth=1) +
    circle(tracks[["hg"]], geom = "text", aes(label = seqnames),
           vjust = 0, size = 3, radius=38) +
    labs(title=tracks[["id"]])
}
