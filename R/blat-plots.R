prettyLayout <- function(npanels, maxperpage=9, nrow=3, ncol=3){
  if(npanels > 9){
    npages <- ceiling(npanels/9)
    return(c(3, 3, npages))
  }
  if(npanels == 1 ){
    nrow=1
    ncol=2
    return(c(1, 1, 1))
  }
  if(npanels <= 9){
    nrow=3
    ncol=3
  }
  if(npanels <= 6){
    nrow=2
  }
  if(npanels <= 4){
    nrow=2
    ncol=2
  }
  if(npanels <= 2){
    nrow=1
    ncol=2
  }
  return(c(ncol, nrow, 1))
}

#' Lattice-style plot of results from blat analysis
#' 
#' @return a \code{trellis} object
#' @export
#' @param blat data.frame of blat records
#' @param layout layout of lattice plot
#' @param ... additional arguments to panel.xyplot
#'
#' @seealso See \code{\link[lattice]{xyplot}} 
plotBlat <- function(blat, layout=c(5, 5, 1), ...){
  colors <- c("gray60", "blue")[blat$is_overlap + 1]
  R2index <- grep("_R2$", blat$Qname)
  isR2 <- seq_len(nrow(blat)) %in% R2index
  cex <- c(0.3, 0.4)[blat$is_overlap+ 1]
  pch <- rep(20, nrow(blat))
  size_of_target <- abs(blat$Tend-blat$Tstart)
  colors[size_of_target > 120] <- "brown"
  ##pch[isR2] <- 24
  ##uid <- paste(blat$id, blat$fusion, sep="_")
  uid <- paste(blat$id, blat$rearrangement, sep="_")
  blat$uid <- uid
  nfusions <- length(unique(uid))
  layout <- prettyLayout(nfusions)
  p <- xyplot(match~tag_index|uid, blat,
              scales=list(x=list(relation="free", labels=NULL, cex=0.7)),
              pch=pch, sizes=cex,
              col=colors, xlab="tags",
              id=blat$id,
              par.strip.text=list(cex=0.5),
              panel=function(x, y, pch, colors, sizes, id, ..., subscripts){
                vlines <- unique(round(x, 0))
                panel.abline(v=vlines, lty=2, col="gray")
                xx <- jitter(x, amount=0.1)
                colors <- colors[subscripts]
                sizes <- sizes[subscripts]
                pch <- pch[subscripts]
                ##id <- paste(unique(id[subscripts]), collapse=",")
                panel.abline(h=90, lty=2, col="gray")
                panel.xyplot(xx, y, col=colors, cex=sizes, pch=pch, ...)

              },
              ##layout=c(5, 5, pages),
              layout=layout,
              as.table=TRUE, ...)
  return(p)
}
