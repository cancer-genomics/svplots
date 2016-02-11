setGeneric("interstitialLegend", function(object, accent) standardGeneric("interstitialLegend"))
setGeneric("seq_along2", function(along.with) standardGeneric("seq_along2"))

setGeneric("boxIdiogram", function(object) standardGeneric("boxIdiogram"))


setClass("IdiogramParams", representation(seqnames="character",
                                          seqlengths="numeric",
                                          unit="character",
                                          genome="character",
                                          box="list"))

IdiogramParams <- function(seqnames=character(),
                           seqlengths=numeric(), unit="kb", genome="hg19",
                           box=list(xlim=ILimit(), color="blue")){
  new("IdiogramParams", seqnames=seqnames, seqlengths=seqlengths, unit=unit,
      genome=genome, box=box)
}

setClass("ILimit", contains="IRanges")

setMethod("boxIdiogram", "IdiogramParams", function(object) object@box)
setMethod("chromosome", "IdiogramParams", function(object) object@seqnames)
setMethod("genome", "IdiogramParams", function(x) x@genome)
setMethod("seqlengths", "IdiogramParams", function(x) x@seqlengths)
setMethod("seqnames", "IdiogramParams", function(x) x@seqnames)

ILimit <- function(...)  as(IRanges(...), "ILimit")

#' Return the endpoints of an IRanges object as a numeric vector
#'
#' @keywords internal
#' @aliases range,ILimit-method
#' @rdname range-method
setMethod("range", "ILimit", function(x, ...) c(start(x), end(x)))

setMethod("seq_along2", "ILimit", function(along.with){
  seq(start(along.with), end(along.with), 1)
})


.find_xlim_percent <- function(g, percent=0.05){
  wd <- width(g)
  w <- wd/percent
  d <- (w-wd)*1/2
  st <- max(start(g)[1]-d, 1)
  en <- min(end(g)[1]+d, seqlengths(g)[chromosome(g)])

  d1 <- start(g)-st
  d2 <- en-end(g)
  if(d1 & d2 > 0){
    d3 <- min(d1, d2)
    st <- start(g)[1]-d3
    en <- end(g)[1]+d3
  }
  lim <- as.integer(c(st, en))
  ILimit(start=lim[1], end=lim[2])
}

xlimTagDensity <- function(si, xlim, percent=0.1){
  g <- GRanges(seqnames(si), IRanges(xlim[1], xlim[2]))
  seqlevels(g,force=TRUE) <- chromosome(g)
  seqlengths(g) <- seqlengths(si)[seqlevels(g)]
  .find_xlim_percent(g, percent)
}

.ideogramParams <- function(object, params){
  xlim <- params[["xlim"]]
  xlim_tagd <- xlimTagDensity(seqinfo(object)[chromosome(object)[1], ], xlim, 0.1)
  chrom <- chromosome(object)[1]
  iparams <- IdiogramParams(seqnames=chrom,
                            genome="hg19",
                            seqlengths=seqlengths(object)[chrom],
                            box=list(xlim=xlim_tagd, color="blue", lwd=2))
  iparams
}

.subset_pview <- function(object, pview, xlim){
  g <- GRanges(chromosome(object)[1], IRanges(start(xlim), end(xlim)))
  ##seqlevels(segs, force=TRUE) <- chromosome(object)[1]
  seqlevels(g, force=TRUE) <- chromosome(object)[1]
  i <- subjectHits(findOverlaps(g, rowRanges(pview), maxgap=start(g)-start(xlim)))
  pview <- pview[i, ]
  pview
}

tagdensity_filter_viewports <- function(xscale, ylim=c(-3,2)){
  layout <- grid.layout(nrow=2, ncol=1, heights=unit(c(1, 2), c("null", "null")))
  vlayout <- viewport(name="tagd_filters_layout",
                      x=0, y=0,
                      width=unit(1, "npc"),
                      height=unit(1, "npc"),
                      just=c("left", "bottom"),
                      layout=layout)

  vpl1 <- viewport(name="tagd_layout",
                   layout.pos.row=2,
                   layout.pos.col=1)

  tag_density_vp <- viewport(name="tag_density_vp",
                             x=unit(0.1, "npc"),
                             y=unit(0.1, "npc"),
                             width=unit(0.9, "npc"),
                             height=unit(0.9, "npc"),
                             just=c("left", "bottom"),
                             xscale=xscale,
                             yscale=ylim,
                             clip="on")

  axisvp <- viewport(name="axisvp",
                     x=unit(0.1, "npc"),
                     y=unit(0.1, "npc"),
                     width=unit(0.9, "npc"),
                     height=unit(0.9, "npc"),
                     just=c("left", "bottom"),
                     yscale=ylim, xscale=xscale,
                     clip="off")

  axislabel <- viewport(name="axislabel",
                        x=unit(0.05, "npc"),
                        y=unit(0.1, "npc"),
                        width=unit(0.05, "npc"),
                        height=unit(0.9, "npc"),
                        just=c("left", "bottom"),
                        clip="on")

  vpl2 <- viewport(name="filter_layout",
                  layout.pos.row=1,
                  layout.pos.col=1)

  filter_vp <- viewport(name="filter_vp",
                        x=unit(0.1, "npc"),
                        y=unit(0.05, "npc"),
                        width=unit(0.9, "npc"),
                        height=unit(0.95, "npc"),
                        just=c("left", "bottom"),
                        xscale=xscale,
                        yscale=unit(c(0,1), "native"),
                        clip="on")

  filteraxis <- viewport(name="axisvp",
                         x=unit(0.05, "npc"),
                         y=unit(0.05, "npc"),
                         width=unit(0.05, "npc"),
                         height=unit(0.95, "npc"),
                         just=c("left", "bottom"),
                         yscale=c(0,1),
                         clip="off")

  filter_labels <- viewport(name="filter_labels",
                            x=unit(0.1, "npc"),
                            y=unit(0.05, "npc"),
                            width=unit(1, "npc"),
                            height=unit(0.95, "npc"),
                            just=c("left", "bottom"),
                            xscale=xscale,
                            yscale=unit(c(0,1), "native"),
                            clip="off")

  vp <- viewport(name="filter_and_tagd_lay",
                 layout.pos.row=c(1,2),
                 layout.pos.col=1)

  filter_tagd <- viewport(name="filter_and_tagd",
                          x=unit(0.1, "npc"),
                          y=unit(0.1,  "npc"),
                          width=unit(0.9, "npc"),
                          height=unit(0.9, "npc"),
                          just=c("left", "bottom"),
                          xscale=xscale,
                          clip="off")



  vp_list <- vpList(vlayout=vlayout,
                    vpl1=vpl1,
                    tdvp=tag_density_vp,
                    axisvp=axisvp,
                    axislabel=axislabel,
                    vpl2=vpl2,
                    filter_vp=filter_vp,
                    filteraxis=filteraxis,
                    filter_labels=filter_labels,
                    vp=vp,
                    filter_tagd=filter_tagd)
  vp_list
}

.draw_ideogram <- function(object, params, vpdata, vps, idiogram){
  pview <- params[["pview"]]
  filterList <- params[["filterList"]]
  zoom.out <- params[["zoom.out"]]
  xlim <- params[["xlim"]]
  segs <- params[["segs"]]
  cex <- params[["cex"]]
  accent <- params[["accent"]]
  xlim <- xlimTagDensity(seqinfo(object)[chromosome(object)[1], ], xlim, 0.1)
  xscale <- range(c(start(rowRanges(pview)), end(rowRanges(pview))))
  vps2 <- tagdensity_filter_viewports(xscale, ylim=c(-3,2))
  start(xlim) <- xscale[1]
  end(xlim) <- xscale[2]
  
  pushViewport(vpdata)
  ivp <- vps[["idiogram"]]
  pushViewport(ivp)
  grid.draw(idiogram)
  seekViewport("idiogram_vp")
  x <- unit(xscale[1], "native")
  w <- unit(diff(xscale), "native")
  grid.rect(x=x,
            width=w,
            y=unit(1, "npc"),
            height=unit(1, "npc"),
            just=c("left", "top"),
            gp=gpar(fill=addalpha("white", 0.2)))

}

yaxisLabels <- function(ylim, n=6, logscale=TRUE){
  exponent <- pretty(seq(ylim[1], ylim[2], 1), n=n)
  if(!logscale){
    exponent <- exponent[exponent %% 1 == 0]
    is_notpositive <- exponent <= 0
    denom <- 2^(-exponent[is_notpositive])
    ylabels <- paste0("1/", denom)
    if(any(!is_notpositive)){
      ylabels <- c(ylabels, 2^exponent[!is_notpositive])
    }
    ylabels[ylabels=="1/1"] <- 1
    ylabels[1] <- paste0("<", ylabels[1])
  }  else ylabels <- exponent
  at <- exponent
  list(at=at, ylabels=ylabels)
}

subsetFilterList <- function(object, granges){
  object <- lapply(object, function(x, granges) subsetByOverlaps(x, granges), granges=granges)
}

plotTagDensityComplex2 <- function(object,
                                   view,
                                   percent=0.1,
                                   filterList,
                                   zoom.out=1,
                                   xlim,
                                   ylim=c(-3, 2),
                                   logy=TRUE,
                                   segs,
                                   accent,
                                   vps,
                                   params,
                                   ...){
  pch <- params[["pch"]]
  col <- params[["pch_color"]]
  cex <- params[["cex"]]
  seg_col <- params[["seg_color"]]
  g <- GRanges(chromosome(object)[1], IRanges(start(xlim), end(xlim)))
  seqlevels(segs, force=TRUE) <- chromosome(object)[1]
  k <- subjectHits(findOverlaps(g, rowRanges(view), maxgap=start(g)-start(xlim)))
  view <- view[k, ]
  xscale <- range(c(start(rowRanges(view)), end(rowRanges(view))))
  td <- assays(view)[,1]
  td <- threshold(td, lim=ylim)
  segs$seg.mean <- threshold(segs$seg.mean, lim=ylim)
  locs <- pretty(range(xlim), n=10)
  labels <- prettyNum(locs/1000, big.mark=",")
  yaxis <- yaxisLabels(ylim, n=6, logscale=TRUE)
  xlim_g <- GRanges(seqnames(g), IRanges(start(xlim), end(xlim)))
  filterList2 <- subsetFilterList(filterList, xlim_g)
  any_filters <- sum(elementLengths(filterList2)) > 0
  cnv <- variant(object)
  pushViewport(vps[["vlayout"]])
  pushViewport(vps[["vpl1"]])
  pushViewport(vps[["tdvp"]])
  grid.rect(gp=gpar(col="gray"))
  grid.points(x=start(rowRanges(view)),
              y=td, pch=pch,
              default.units="native",
              gp=gpar(col=col, cex=cex))
  grid.segments(x0=start(segs),
                x1=end(segs),
                y0=segs$seg.mean,
                y1=segs$seg.mean,
                default.units="native", gp=gpar(lwd=1.5, col=seg_col))
  upViewport()
  pushViewport(vps[["axisvp"]])
  grid.yaxis(gp=gpar(cex=0.6), main=FALSE)
  at <- pretty(xscale, n=8)
  at <- at[at >= xscale[1] & at <= tail(xscale, 1)]
  labels <- prettyNum(at/1000, big.mark=",")
  grid.xaxis(at=at, label=labels, gp=gpar(cex=0.6))
  upViewport()
  ##tg <- textGrob("log2 tag\ndensity", gp=gpar(cex=0.7))
  pushViewport(vps[["axislabel"]])
  grid.rect(gp=gpar(fill="wheat", col=NA))
  grid.text("log2 tag\ndensity", rot=90, gp=gpar(cex=0.7))
  upViewport(2)
  pushViewport(vps[["vpl2"]])
  pushViewport(vps[["filter_vp"]])
  ##grid.rect(gp=gpar(col="gray"))
  ys <- seq(0, 0.9, length.out=length(filterList))
  h <- 0.1
  grid.segments(y0=unit(ys+h/2, "native"),
                y1=unit(ys+h/2, "native"),
                gp=gpar(lty=2, col="gray"))
  ##ys <- rep(ys, elementLengths(filterList))
  ##h <- diff(ys)[1]*1/2
  ##filters <- unlist(GRangesList(lapply(filterList, reduce)))
  ##seqlevels(filters, force=TRUE) <- chromosome(g)[1]
  ##hits <- subjectHits(findOverlaps(g, filters))
  for(k in seq_along(filterList)){
    f <- filterList[[k]]
    yy <- unit(ys[k], "native")
    hits <- subjectHits(findOverlaps(g, f))
    if(length(hits) > 0){
      st <- start(f)[hits]
      en <- end(f)[hits]
      grid.rect(x=st,
                width=en-st,
                y=yy,
                height=h,
                default.units="native",
                gp=gpar(col=NA, fill="gray40"),
                just=c("left", "bottom"))
    }
  }
  upViewport()
  pushViewport(vps[["filteraxis"]])
  grid.rect(gp=gpar(fill="wheat", col=NA))
  grid.text("filters", rot=90, gp=gpar(cex=0.7))
  yy <- unique(ys)
  grid.yaxis(at=yy,
             label=names(filterList),
             gp=gpar(cex=0.7))
  upViewport(2)
  pushViewport(vps[["vp"]])
  pushViewport(vps[["filter_tagd"]])
  grid.rect(x=start(cnv),
            width=end(cnv)-start(cnv),
            y=unit(-0.01, "npc"),
            height=unit(1, "npc"),
            default.units="native",
            gp=gpar(fill=accent, col=NA),
            just=c("left", "bottom"))
  upViewport(2)
}

.draw_tagdensity <- function(object, params, vpdata, vps){
  pview <- params[["pview"]]
  filterList <- params[["filterList"]]
  zoom.out <- params[["zoom.out"]]
  xlim <- params[["xlim"]]
  segs <- params[["segs"]]
  cex <- params[["cex"]]
  accent <- params[["accent"]]
  xlim <- xlimTagDensity(seqinfo(object)[chromosome(object)[1], ], xlim, 0.1)
  xscale <- range(c(start(rowRanges(pview)), end(rowRanges(pview))))
  vps2 <- tagdensity_filter_viewports(xscale, ylim=c(-3,2))
  start(xlim) <- xscale[1]
  end(xlim) <- xscale[2]  
  upViewport(0)
  pushViewport(vpdata)
  vp <- vps[["tagdens"]]
  pushViewport(vp)
  plotTagDensityComplex2(object,
                         view=pview,
                         filterList=filterList,
                         zoom.out=zoom.out,
                         xlim=xlim,
                         segs=segs, 
                         accent=accent,
                         vps=vps2,
                         params=params)
}

.ideogram_plus_tagdensity <- function(object, params, vpdata, vps, idiogram){
  .draw_ideogram(object, params, vpdata, vps, idiogram)
  .draw_tagdensity(object, params, vpdata, vps)
}

grid.idiogram <- function(params){
  coords <- idiogram_coordinates(params)
  colors <- getOption("biovizBase")$cytobandColor
  colors <- colors[coords$gieStain]

  top <- 0.75
  bot <- 0.25
  h <- top-bot
  deltas <- (coords$xright-coords$xleft)/4
  ##deltas <- with(coords, xright-xleft)/4

  taper_right <- coords$taper_right*h
  taper_left <- coords$taper_left*h
  y <- cbind(bot+coords$taper_left*h, bot, bot, bot+taper_right,
             top-taper_right, top , top, top-taper_left)

  x <- cbind(coords$xleft, coords$xleft+deltas, coords$xright-deltas,
             coords$xright, coords$xright,
             coords$xright-deltas,
             coords$xleft+deltas,
             coords$xleft)
  
##  x <- with(coords, cbind(xleft, xleft+deltas, xright-deltas,
##                          xright, xright,
##                          xright-deltas,
##                          xleft+deltas,
##                          xleft))
  yy <- as.numeric(t(y))
  xx <- as.numeric(t(x))
  id <- rep(1:nrow(y), each=ncol(y))

  polygonGrob(x=unit(xx, "native"), y=unit(yy, "npc"),
              gp=gpar(fill=colors),
              id=id,
              vp=viewport(name="idiogram_vp",
                          xscale=c(0, max(coords$xright))),
              name="idiogram")
}

getCytoband <- function(params){
  path <- system.file("extdata", package="SNPchip", mustWork=TRUE)
  cytoband <- read.table(file.path(path, paste("cytoBand_",
                                               genome(params), ".txt", sep="")),
                         as.is=TRUE, header=FALSE)
  colnames(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
  cytoband <- cytoband[cytoband[, "chrom"] %in% seqnames(params), ]
  rownames(cytoband) <- as.character(cytoband[, "name"])
  return(cytoband)
}

.drawStalk <- function(starts, ends){
  deltas <- (ends-starts)/3
  y0 <- rep(0.2, length(starts))
  y1 <- rep(0.8, length(starts))
  lsegments(starts+deltas, y0, starts+deltas, y1)
  lsegments(ends+deltas, y0, ends+deltas, y1)
}

idiogram_coordinates <- function(params){
  build <- genome(params)
  cytoband <- getCytoband(params)
  cytoband_p <- cytoband[grep("^p", rownames(cytoband), value=TRUE), ]
  cytoband_q <- cytoband[grep("^q", rownames(cytoband), value=TRUE), ]
  N <- nrow(cytoband)
  p.bands <- nrow(cytoband_p)
  left_cuts <- c(1, p.bands+1)
  right_cuts <- c(N, p.bands)
  cut.right <- cut.left <- rep(FALSE, N)
  cut.left[left_cuts] <- TRUE
  cut.right[right_cuts] <- TRUE
  ## this is a "stalk", do not draw box. Draw two vertical lines instead
  is_stalk <- cytoband[, "gieStain"] == "stalk"
  index_stalk <- which(is_stalk)
  cut.right[index_stalk - 1] <- TRUE
  cut.left[index_stalk  + 1] <- TRUE
  colors <- .cytobandColors(cytoband[, "gieStain"])
  xx <- c(0, cytoband[nrow(cytoband), "end"])
##yy <- c(0,1)
  starts <- cytoband[, "start"]
  ends <- cytoband[, "end"]
  if(any(is_stalk)) .drawStalk(starts[is_stalk], ends[is_stalk])
  taper_right <- taper_left <- rep(0, length(starts))
  taper_left[cut.left] <- 0.15
  taper_right[cut.right] <- 0.15

  data.frame(xleft=starts,
             xright=ends,
             taper_right=taper_right,
             taper_left=taper_left,
             gieStain=cytoband$gieStain,
             stringsAsFactors=FALSE)
}

.connect_ideogram <- function(vpdata, vps, vps2, iparams, xscale){
  pushViewport(vpdata)
  pushViewport(vps[["idiogram"]])
  coords <- idiogram_coordinates(iparams)
  vp2 <- viewport(name="idiogram_vp", xscale=c(0, max(coords$xright)))
  pushViewport(vp2)
  x <- unit(xscale[1], "native")
  grid.move.to(x=x - unit(1, "mm"),
               y=unit(0, "npc") - unit(2, "mm"))
  upViewport(0)
  pushViewport(vpdata)
  pushViewport(vps[["tagdens"]])
  pushViewport(vps2[["vlayout"]])
  pushViewport(vps2[["vp"]])
  pushViewport(vps2[["filter_tagd"]])
  ltg <- lineToGrob(x=x+unit(2, "mm"),
                    y=unit(1, "npc") + unit(2, "mm"))
  grid.draw(ltg)
  x2 <- unit(xscale[2], "native")
  grid.move.to(x=x2-unit(2, "mm"),
               y=unit(1, "npc") + unit(2, "mm"))
  upViewport(0)
  pushViewport(vpdata)
  pushViewport(vps[["idiogram"]])
  pushViewport(vp2)
  ltg <- lineToGrob(x=x2+unit(1, "mm"),
                    y=unit(0, "npc") - unit(2, "mm"))
  grid.draw(ltg)  
}

grid_readpairs <- function(rp, del, accent, legend_labels, ...){
  ix <- order(start(first(rp)))
  rp <- rp[ix]
  L <- length(rp)
  i <- which(start(first(rp)) < start(last(rp)))
  j <- seq_along(rp)[-i]
  y <- seq_along(rp)
  ylim <- current.panel.limits()$ylim
  if(L < 10e3){
    if(length(i) > 0){
      grid.segments(x0=unit(end(first(rp))[i], "native"),
                    y0=unit(i, "native"),
                    x1=unit(start(last(rp))[i], "native"),
                    y1=unit(i, "native"),
                    gp=gpar(col="gray"))
    }
    if(length(j) > 0){
      grid.segments(x0=unit(end(last(rp))[j], "native"),
                    y0=unit(j, "native"),
                    x1=unit(start(first(rp))[j], "native"),
                    y1=unit(j, "native"),
                    gp=gpar(col="gray"))
    }
    w <- end(first(rp))-start(first(rp))
    grid.rect(x=unit(start(first(rp)), "native"),
              y=unit(y-0.2, "native"),
              width=unit(w, "native"),
              height=unit(0.4, "native"),
              gp=gpar(col="black", fill="black"))
    w <- end(last(rp))-start(last(rp))
    grid.rect(x=unit(start(last(rp)), "native"),
              y=unit(y-0.2, "native"),
              width=unit(w, "native"),
              height=unit(0.4, "native"),
              gp=gpar(col="royalblue", fill="royalblue"))
  } else {
    d <- start(last(rp)) - end(first(rp))
    j <- which(d < 0)
    if(length(j) > 0)
      d[j] <- end(first(rp))[j]-start(last(rp))[j]
    big <- d > 3e3
    ##panel.points(start(first(rp)), y, pch=20, cex=0.2)
    ##grid.points(start(first(rp)), y, gp=gpar(pch=20, cex=0.2, col="blue"))
    grid.rect(x=unit(start(first(rp)), "native"),
              width=unit(100, "native"),
              y=unit(seq_along(d), units="native")-unit(0.2, "native"),
              height=unit(0.4, "native"),
              gp=gpar(col="blue"), just=c("left", "bottom"))

    ##browser()
    ii <- intersect(which(big), i)
    jj <- intersect(which(big), j)
    if(length(ii) > 0){
      grid.segments(x0=unit(end(first(rp))[ii], "native"),
                    y0=unit(ii, "native"),
                    x1=unit(start(last(rp))[ii], "native"),
                    y1=unit(ii, "native"), gp=gpar(col="gray"))
    }
    if(length(jj) > 0){
      grid.segments(x0=unit(end(last(rp))[jj], "native"),
                    y0=unit(jj, "native"),
                    x1=unit(start(first(rp))[jj], "native"),
                    y1=unit(jj, "native"), gp=gpar(col="gray"))
    }
  }
  ##browser()
  grid.rect(x=unit(start(del), "native"),
            y=unit(0, "npc"),
            height=unit(1, "npc"),
            width=unit(width(del), "native"),
            gp=gpar(fill=accent, col="transparent"),
            just=c("left", "bottom"))
}


.plotReadPairsComplex <- function(object, variant.indices, ##...,
                                  zoom.out=1,
                                  unit="kb",
                                  accent,
                                  xlim,
                                  legend_labels){
  if(missing(variant.indices)) variant.indices <- seq_along(variant(object))
  locs <- pretty(range(xlim), n=10)
  locs <- locs[locs >= xlim[1] & locs <= xlim[2]]
  labels <- prettyNum(locs/1000, big.mark=",")
  ##
  ## Tricky.  Want all read pairs belonging to the group
  ##
  rp <- readPairs(object)
  L <- length(rp)
  if(L > 25e3){
    d <- start(last(rp)) - end(first(rp))
    j <- which(d < 0)
    d[j] <- end(first(rp))[j]-start(last(rp))[j]
    big <- d > 3e3
    thin <- seq(1, length(rp), length.out=25e3)
    thin <- sort(unique(c(thin, which(big))))
    rp <- rp[thin]
  }
  ylim <- c(-1 * length(rp)*0.01, length(rp)+0.5)
  xscale <- range(xlim)
  vprp <- viewport(x=unit(0.09, "npc"),
                   y=unit(0.15, "npc"),
                   width=0.76, height=0.7,
                   xscale=xscale,
                   yscale=ylim, name="readpairs",
                   just=c("left", "bottom"))


  axisvp <- viewport(name="axisvp",
                     x=unit(0.09, "npc"),
                     y=unit(0.15, "npc"),
                     width=unit(0.76, "npc"),
                     height=unit(0.7, "npc"),
                     just=c("left", "bottom"),
                     yscale=ylim, xscale=xscale,
                     clip="off")

  at <- pretty(c(1, length(rp)), n=6)
  axislabel <- viewport(name="axislabel",
                        x=unit(0.06, "npc"),
                        y=unit(0.15, "npc"),
                        width=unit(0.04, "npc"),
                        height=unit(0.7, "npc"),
                        just=c("left", "bottom"),
                        yscale=c(0, max(at)),
                        clip="off")
  ##tg <- textGrob("read pairs", x=-0.15, y=0.5)
  ##grid.draw(editGrob(tg, vp=viewport(angle=90)))
  ##grid.text("read pairs", x=-0.1, y=0.5, rot=90)
  ##grid.points(swf$TEMP, swf$DMC, gp=gpar(col="grey"))
  pushViewport(vprp)
  grid.rect(gp=gpar(col="gray"))
  grid.text("kb", x=0.5, y=-0.2)
  grid_readpairs(rp, del=variant(object)[variant.indices], accent=accent,
                 legend_labels=legend_labels)
  grid.xaxis(at=locs, label=labels, gp=gpar(cex=0.7))
  upViewport()
  pushViewport(axislabel)
  grid.rect(gp=gpar(fill="wheat", col=NA))
  grid.text("read pairs", rot=90, gp=gpar(cex=0.8))
  grid.yaxis(at=at, label=at, gp=gpar(cex=0.7))
  return(vprp)
}


.plot_readpairs <- function(object_all, object, vpdata,
                            vps, xlim, legend_labels, zoom.out,
                            accent, group.index){

  ##
  ## Tricky. If we only passed the object, readpairs in neighboring
  ## regions would not be plotted
  ##
  object2 <- object_all[groupedVariant(object_all)==groupedVariant(object)[1]]
  variant_indices <- match(names(variant(object_all))[group.index], names(variant(object2)))
  .plotReadPairsComplex(object2,
                        variant.indices=variant_indices,
                        zoom.out=zoom.out,
                        accent=accent, xlim=xlim,
                        legend_labels=legend_labels)
}


Viewports <- function(){
  ##idiogram
  v0 <- viewport(name="idiogram_lay",
                 x=unit(0.15, "npc"), y=unit(0.89, "npc"),
                 width=unit(0.65, "npc"),
                 height=unit(0.1, "npc"), just=c("left", "bottom"))

  filterVP <- viewport(x=unit(0.01, "npc"), y=unit(0.8, "npc"),
                       width=unit(0.7, "npc"),
                       height=unit(0.1, "npc"),
                       just=c("left", "bottom"))

  ## tag density
  v1 <- viewport(x=unit(0.01, "npc"), y=unit(0.6, "npc"),
                 width=unit(0.8, "npc"),
                 height=unit(0.2, "npc"), just=c("left", "bottom"))

  ## paired reads
  v2 <- viewport(x=unit(0.005, "npc"), y=unit(0.01, "npc"),
                 width=unit(0.95, "npc"),
                 height=unit(0.55, "npc"), just=c("left", "bottom"),
                 clip="off")

  ## interstitial legend
  v3 <- viewport(x=unit(0.8, "npc"), y=unit(0.7, "npc"),
                 width=unit(0.15, "npc"),
                 height=unit(0.3, "npc"), just=c("left", "bottom"))

  ## intrachromosomal legend
  v4 <- viewport(x=unit(0.8, "npc"), y=unit(0.6, "npc"),
                 width=unit(0.2, "npc"),
                 height=unit(0.1, "npc"), just=c("left", "bottom"))

  ## interchromosomal legend
  v5 <- viewport(x=unit(0.8, "npc"), y=unit(0.2, "npc"),
                 width=unit(0.2, "npc"),
                 height=unit(0.2, "npc"), just=c("left", "bottom"))


  ## inset
  vpR1 <- viewport(x=unit(0.1, "npc"), y=unit(0.55, "npc"),
                   width=unit(0.15, "npc"),
                   height=unit(0.35, "npc"), just=c("left", "bottom"))
  vpR2 <- viewport(x=unit(0.23, "npc"), y=unit(0.55, "npc"),
                   width=unit(0.15, "npc"),
                   height=unit(0.35, "npc"), just=c("left", "bottom"))

  vpR1_2 <- viewport(x=unit(0.65, "npc"), y=unit(0.15, "npc"),
                     width=unit(0.15, "npc"),
                     height=unit(0.35, "npc"), just=c("left", "bottom"))
  vpR2_2 <- viewport(x=unit(0.78, "npc"), y=unit(0.15, "npc"),
                     width=unit(0.15, "npc"),
                     height=unit(0.35, "npc"), just=c("left", "bottom"))

  list(idiogram=v0, tagdens=v1,
       tags=v2,
       interstitial=v3,
       intra=v4,
       inter=v5,
       vpR1, vpR2, vpR1_2, vpR2_2,
       filter=filterVP)
}

.overlapping_segs <- function(object, legend_labels, legvp, label, accent){
  ys <- .legend_locations(object)
  ##h <- diff(ys)[1]*1/2
  h <- abs(diff(ys)[1])
  J <- seq_along(ys)
  ##any_overlapping <- any(calls(object)=="OverlappingHemizygous+")
  ##if(any_overlapping){
  k <- which(calls(object)=="OverlappingHemizygous+")
  calls(object)[k] <- "hemizygous+"
  accent[k] <- "transparent"
  ##}
  ##
  ## Return if no overlapping segs
  ##
  ##if(!any_overlapping) return(invisible())
  J <- J[-k]
  ##
  ## In order to modify the read pairs figure, we need to identify
  ## the viewport that lattice used to draw the tags.
  ##
  st <- start(variant(object))[k]
  en <- end(variant(object))[k]
  ## dashed lines on panel
  grid.segments(x0=unit(c(st, en), "native"),
                x1=unit(c(st, en), "native"),
                y0=unit(rep(0,2), "npc"),
                y1=unit(rep(1,2), "npc"),
                gp=gpar(col="gray", lty=2))
  ##
  ## dendogram-style axis
  ##
  y1 <- unit(rep(1, 2), "npc") + unit(as.integer(factor(k))*4, "mm")
  grid.segments(x0=unit(c(st, en), "native"),
                x1=unit(c(st, en), "native"),
                y0=unit(rep(1, 2), "npc") + unit(1, "mm"),
                y1=y1,
                gp=gpar(col="black"))
  grid.segments(x0=unit(st, "native"),
                x1=unit(en, "native"),
                y0=y1,
                y1=y1,
                gp=gpar(col="black"))
  mid <- (st+en)/2
  ## Add labels to dendrogram
  for(ii in seq_along(k)){
    kk <- k[ii]
    tg <- textGrob(legend_labels[kk], x=unit(mid[ii], "native"),
                   y=y1[ii]+unit(1, "mm"), just="center")
    cg <- circleGrob(x=unit(mid[ii], "native"),
                     y=y1[ii]+unit(1, "mm"),
                     r=0.9*grobWidth(tg))
    circledText <- gTree(children=gList(tg, cg))
    grid.draw(editGrob(circledText, gp=gpar(col="gray", fill="white")))
    grid.draw(editGrob(tg, gp=gpar(col="black")))
  }
  ## Add legends for overlapping regions
  for(m in k){
    upViewport(0)
    ##      pushViewport(rightmargin)
    pushViewport(legvp)
    vp <- viewport(x=0, y=ys[m], width=unit(1, "npc"),
                   height=unit(h, "npc"),
                   just=c("left", "bottom"))
    pushViewport(vp)
    grid.text(legend_labels[m],
              x=unit(0, "npc")+unit(1, "mm"),
              y=unit(1, "npc")-unit(5, "mm"),
              gp=gpar(cex=1))
    interstitialLegend(object[m], accent=accent[m])
  }
  J <- seq_along(ys)[-k]
  if(length(J) > 0){
    for(j in J){
      st <- start(variant(object))[j]
      en <- end(variant(object))[j]
      yy <- unit(1, "npc") - unit(1, "mm")
      mid <- (st+en)/2
      seekViewport("readpairs")
      tg <- textGrob(legend_labels[j],
                     x=unit(mid, "native"),
                     y=yy, just="center")
      cg <- circleGrob(x=unit(mid, "native"),
                       y=yy,
                       r=0.9*grobWidth(tg))
      circledText <- gTree(children=gList(tg, cg))
      grid.draw(editGrob(circledText, gp=gpar(col="gray", fill="white")))
      grid.draw(editGrob(tg, gp=gpar(col="black")))
      upViewport(0)
      pushViewport(legvp)
      vp <- viewport(x=0, y=ys[j], width=unit(1, "npc"),
                     height=unit(h, "npc"),
                     just=c("left", "bottom"))
      pushViewport(vp)
      grid.text(legend_labels[j],
                x=unit(0, "npc")+unit(1, "mm"),
                y=unit(1, "npc")-unit(5, "mm"),
                gp=gpar(cex=1))
      interstitialLegend(object[j], accent=accent[j])
      upViewport()
    }
  }
  upViewport(0)
  grid.text(label, x=unit(0.01, "npc"), y=unit(0.98, "npc"),
            gp=gpar(cex=1.1), just=c("left", "top"))  
}

.legend_locations <- function(object){
  L <- length(object)
  ys <- seq(0, 1, length.out=L+1)
  ys <- ys[-length(ys)]
  ys <- rev(ys)
  h <- abs(diff(ys)[1])
  if(length(ys)==1){
    ys <- 0.5
    h <- 0.35
  }
  ys
}

setMethod("interstitialLegend", "StructuralVariant", function(object, accent){
  jun <- variant(object)
  fold_change <- round(2^copynumber(object), 2)
  size <- prettyNum(round(width(jun)/1e3, 2), big.mark=",")
  st <- prettyNum(round(start(jun)/1e3, 2), big.mark=",")
  en <- prettyNum(round(end(jun)/1e3, 2), big.mark=",")
  nRPs <- length(improper(object))
  labels <- paste0("size (kb)     :", size, "\n",
                   "fold change   :", fold_change, "\n",
                   "rearranged RPs:", nRPs, "\n",
                   "start (kb)    :", st, "\n",
                   "end   (kb)    :", en, "\n")
                   ##"id            :", names(jun), "\n")
  grid.rect(gp=gpar(fill=accent, col="transparent"))
  ## title of legend
  grid.text(calls(object), x=unit(0.5, "npc"),
            y=unit(1, "npc"),
            just=c("center", "top"),
            gp=gpar(font=2))
  grid.text(labels, x=unit(0.05, "npc"),
            y=unit(0.75, "npc"), just=c("left", "top"),
            gp=gpar(cex=0.6, fontfamily="mono"))
})

.legend_deletions <- function(object, legend_labels, legvp, label, accent){
  ys <- .legend_locations(object)
  J <- seq_along(ys)
  h <- abs(diff(ys)[1])
  if(length(ys)==1){
    ys <- 0.5
    h <- 0.35
  }
  if(length(J) > 0){
    for(j in J){
      st <- start(variant(object))[j]
      en <- end(variant(object))[j]
      yy <- unit(1, "npc") - unit(1, "mm")
      mid <- (st+en)/2
      seekViewport("readpairs")
      tg <- textGrob(legend_labels[j],
                     x=unit(mid, "native"),
                     y=yy, just="center")
      cg <- circleGrob(x=unit(mid, "native"),
                       y=yy,
                       r=0.9*grobWidth(tg))
      circledText <- gTree(children=gList(tg, cg))
      grid.draw(editGrob(circledText, gp=gpar(col="gray", fill="white")))
      grid.draw(editGrob(tg, gp=gpar(col="black")))
      upViewport(0)
      pushViewport(legvp)
      vp <- viewport(x=0, y=ys[j], width=unit(1, "npc"),
                     height=unit(h, "npc"),
                     just=c("left", "bottom"))
      pushViewport(vp)
      grid.text(legend_labels[j],
                x=unit(0, "npc")+unit(1, "mm"),
                y=unit(1, "npc")-unit(5, "mm"),
                gp=gpar(cex=1))
      interstitialLegend(object[j], accent=accent[j])
      upViewport()
    }
  }
  upViewport(0)
  grid.text(label, x=unit(0.01, "npc"), y=unit(0.98, "npc"),
            gp=gpar(cex=1.1), just=c("left", "top"))
}

gridComplex <- function(object, group.index, legend.index, params){
  filterList <- params[["filterList"]]
  zoom.out <- params[["zoom.out"]]
  segs <- params[["segs"]]
  pview <- params[["pview"]]
  aview <- params[["aview"]]
  accent <- params[["accent"]][legend.index]
  xlim <- params[["xlim"]]
  legend_labels <- params[["legend_labels"]][legend.index]
  
  vps <- Viewports()
  object_all <- object
  object <- object[group.index]
  iparams <- .ideogramParams(object, params)
  idiogram <- grid.idiogram(iparams)
  vpdata <- viewport(x=0, y=0, width=unit(0.8, "npc"), height=unit(1, "npc"),
                     just=c("left", "bottom"))
  legvp <- viewport(x=0.75, y=0, width=unit(0.25, "npc"), height=unit(1, "npc"),
                    just=c("left", "bottom"))
  label <- paste(names(aview), chromosome(object)[1], sep="\n")
  vplist <- list()
  ##
  ## Begin drawing figure
  ##
  xlim_tagd <- xlimTagDensity(seqinfo(object)[chromosome(object)[1], ], xlim, 0.1)
  pview <- params[["pview"]] <- .subset_pview(object, pview, xlim_tagd)
  xscale <- range(c(start(rowRanges(pview)), end(rowRanges(pview))))
  vps2 <- tagdensity_filter_viewports(xscale, ylim=c(-3,2))
  grid.newpage()
  .ideogram_plus_tagdensity(object, params, vpdata, vps, idiogram)
  upViewport(0)
  .connect_ideogram(vpdata, vps, vps2, iparams, xscale)
  ##print(idiogram, vp=vp, newpage=FALSE)
  upViewport(0)
  pushViewport(vpdata)
  vp <- vps[["tags"]]
  pushViewport(vp)
  vprp <- .plot_readpairs(object_all, object, vpdata, vps,
                          xlim, legend_labels,
                          zoom.out, accent, group.index)
  upViewport()
  pushViewport(vprp)
  any_overlapping <- any(calls(object)=="OverlappingHemizygous+")
  if(any_overlapping){
    .overlapping_segs(object, legend_labels, legvp, label, accent)
  } else{
    .legend_deletions(object, legend_labels, legvp, label, accent)
  }
  upViewport(0)  
}

loadFilters <- function(){
  data(lowMappabilityBins_hg19, package="svfilters", envir=environment())
  data(tx_hg19, package="svfilters", envir=environment())
  data(binAssemblyGaps_hg19, package="svfilters", envir=environment())
  data(lymphoblast_filters_hg19, package="svfilters",
       envir=environment())
  filters <- as.list(lymphoblast_filters_hg19)
  filters$map <- lowMappabilityBins_hg19
  filters$gc <- binAssemblyGaps_hg19
  filters
}

.xlimit_deletion <- function(object, group, zoom.out=1){
  groups <- groupedVariant(object)
  k <- which (groups == group)
  g2 <- expandGRanges(variant(object)[k], width(variant(object)[k])*zoom.out)
  g2 <- GRanges(seqlevels(g2)[1], IRanges(min(start(g2)), max(end(g2))))
  seqlengths(g2) <- seqlengths(variant(object)[k])[1]
  xlim <- c(start(g2), end(g2))
  xlim
}

.colors_deletion <- function(object, group, accent){
  groups <- groupedVariant(object)
  k <- which (groups == group)
  L <- length(k)
  accent <- rep(accent, length.out=L)
  accent
}

.legend_labels <- function(object, group){
  legend_labels <- NULL
  groups <- groupedVariant(object)
  k <- which (groups == group)
  L <- length(k)
  if(L > 1) legend_labels <- LETTERS[seq_len(L)] 
  legend_labels
}

.group_indices <- function(object, group){
  k <- which(groupedVariant(object) == group)
  klist <- split(k, factor(rep(seq_len(length(k)), each=4, length.out=length(k))))
  klist
}

.legend_indices <- function(object, group){
  k <- which (groupedVariant(object) == group)
  L <- length(k)  
  ilist <- split(seq_len(L), factor(rep(seq_len(length(k)), each=4, length.out=length(k))))
}


gridDeletionParams <- function(object,
                               ##pview,
                               ##aview,
                               dirs,
                               id,
                               group,
                               filterList,
                               zoom.out=1,
                               cex=0.3,
                               accent,
                               pch=20,
                               pch_color="gray",
                               gaps_gr=gaps_gr,
                               seg_color="black"){
  segs <- readRDS(file.path(dirs[["0cbs"]], paste0(id, ".rds")))
  frac <- intOverWidth(segs, gaps_gr)
  segs <- segs[ frac < 0.8 ]
  xlim <- .xlimit_deletion(object, group, zoom.out)
  accent <- .colors_deletion(object, group, accent)
  legend_labels <- .legend_labels(object, group)
  klist <- .group_indices(object, group)
  ilist <- .legend_indices(object, group)
  list(filterList=filterList,
       segs=segs,
       zoom.out=zoom.out,
       accent=accent,
       xlim=xlim,
       group.index=klist,
       legend.index=ilist,
       cex=cex,
       legend_labels=legend_labels,
       pch=pch,
       pch_color=pch_color,
       seg_color=seg_color)
}

#' Parameters for plotting deletions
#'
#' 
#' @param sv a \code{StructuralVariant} object
#' @param dirs a \code{DataPaths} object
#' @param id  sample id character string
#' @param group grouping factor for the deletions
#' @param gaps a \code{GRanges} object containing centromeres,
#'   heterochromatin regions, and telomeres
#' 
#' @export
grid_params <- function(sv, dirs, id, group=1, gaps){
  accent <- addalpha(brewer.pal(12, "Paired"), 0.2)
  gray <- addalpha("gray10", alpha=0.2)
  lightblue <- addalpha("lightblue", alpha=0.9)
  grayLight <- addalpha("gray", alpha=0.1)
  gray <- addalpha("gray10", alpha=0.2)
  accent <- addalpha(brewer.pal(12, "Paired"), 0.2)
  ##gaps_gr <- readRDS(file.path(dirs[["extdata"]], "gaps_gr_hg19.rds"))
  ##views <- loadViews(id, dirs)
  ##filters <- loadFilters()
  filters <- listGenomeFilters("hg19")
  ## make names short
  names(filters) <- c("centr", "gaps", "germ", "out", "tx")
  params <- gridDeletionParams(sv,
                               ##pview=views[["lr"]],
                               ##aview=views[["aln"]],
                               id=id,
                               accent=accent,
                               filterList=filters,
                               dirs=dirs,
                               group=group,
                               pch=".",
                               pch_color=gray,
                               gaps_gr=gaps,
                               seg_color=lightblue)
}


#' Plot data supporting deletions
#' 
#' @export 
#' @param object a \code{StructuralVariant} object
#' @param params a list of parameters as created by \code{grid_params}
#' @param pviews a \code{PreprocessViews2} object
#' @seealso \code{\link{grid_params}}
gridPlot2 <- function(object, params, pview){
  klist <- params[["group.index"]]
  ilist <- params[["legend.index"]]
  params$pview <- pview
  for(m in seq_along(klist)){
    p <- gridComplex(object, 
                     group.index=klist[[m]],
                     legend.index=ilist[[m]],
                     params=params)
  }
  p
}


draw_ideogram2 <- function(object, vps, params, pview, subset=TRUE){
  if(missing(vps)) vps <- Viewports()
  group.index <- params[["group.index"]][[1]]
  ##aview <- params[["aview"]]
  object_all <- object
  object <- object[group.index]
  iparams <- .ideogramParams(object, params)
  idiogram <- grid.idiogram(iparams)
  vpdata <- viewport(x=0, y=0, width=unit(0.8, "npc"), height=unit(1, "npc"),
                     just=c("left", "bottom"))
  label <- paste(colnames(pview), chromosome(object)[1], sep="\n")
  vplist <- list()
  xlim_tagd <- xlimTagDensity(seqinfo(object)[chromosome(object)[1], ],
                              params[["xlim"]], 0.1)
  if(subset){
    params$pview <- .subset_pview(object, pview, xlim_tagd)
  } else{
    params$pview <- keepSeqlevels(pview, chromosome(object)[1])
  }
  pview <- params$pview
  xscale <- range(c(start(rowRanges(pview)), end(rowRanges(pview))))
  vps2 <- tagdensity_filter_viewports(xscale, ylim=c(-3,2))  
  grid.newpage()
  .draw_ideogram(object, params, vpdata, vps, idiogram)
}


.draw_filter_tracks <- function(object,
                                view,
                                percent=0.1,
                                filterList,
                                zoom.out=1,
                                xlim,
                                ylim=c(-3, 2),
                                logy=TRUE,
                                segs,
                                accent,
                                vps,
                                params,
                                ...){
  pch <- params[["pch"]]
  col <- params[["pch_color"]]
  cex <- params[["cex"]]
  seg_col <- params[["seg_color"]]
  g <- GRanges(chromosome(object)[1], IRanges(start(xlim), end(xlim)))
  seqlevels(segs, force=TRUE) <- chromosome(object)[1]
  k <- subjectHits(findOverlaps(g, rowRanges(view), maxgap=start(g)-start(xlim)))
  view <- view[k, ]
  ##rowRanges(view) <- rowRanges(view)[k]
  xscale <- range(c(start(rowRanges(view)), end(rowRanges(view))))
  td <- assays(view)[,1]
  td <- threshold(td, lim=ylim)
  segs$seg.mean <- threshold(segs$seg.mean, lim=ylim)
  locs <- pretty(range(xlim), n=10)
  labels <- prettyNum(locs/1000, big.mark=",")
  yaxis <- yaxisLabels(ylim, n=6, logscale=TRUE)
  xlim_g <- GRanges(seqnames(g), IRanges(start(xlim), end(xlim)))
  filterList2 <- subsetFilterList(filterList, xlim_g)
  any_filters <- sum(elementLengths(filterList2)) > 0
  cnv <- variant(object)
  pushViewport(vps[["vlayout"]])
  pushViewport(vps[["vpl2"]])
  pushViewport(vps[["filter_vp"]])
  ys <- seq(0, 0.9, length.out=length(filterList))
  h <- 0.1
  grid.segments(y0=unit(ys+h/2, "native"),
                y1=unit(ys+h/2, "native"),
                gp=gpar(lty=2, col="gray"))
  for(k in seq_along(filterList)){
    f <- filterList[[k]]
    yy <- unit(ys[k], "native")
    hits <- subjectHits(findOverlaps(g, f))
    if(length(hits) > 0){
      st <- start(f)[hits]
      en <- end(f)[hits]
      grid.rect(x=st,
                width=en-st,
                y=yy,
                height=h,
                default.units="native",
                gp=gpar(col=NA, fill="gray40"),
                just=c("left", "bottom"))
    }
  }
  upViewport()
  pushViewport(vps[["filteraxis"]])
  grid.rect(gp=gpar(fill="wheat", col=NA))
  grid.text("filters", rot=90, gp=gpar(cex=0.7))
  yy <- unique(ys)
  grid.yaxis(at=yy,
             label=names(filterList),
             gp=gpar(cex=0.7))
  upViewport(2)
  pushViewport(vps[["vp"]])
  pushViewport(vps[["filter_tagd"]])
  ## this highlights the rearrangement region
  grid.rect(x=start(cnv),
            width=end(cnv)-start(cnv),
            y=unit(-0.01, "npc"),
            height=unit(1, "npc"),
            default.units="native",
            gp=gpar(fill=accent, col=NA),
            just=c("left", "bottom"))
  upViewport(2)
}

.draw_filter_tracks2 <- function(object, params, vpdata, vps){
  pview <- params[["pview"]]
  filterList <- params[["filterList"]]
  zoom.out <- params[["zoom.out"]]
  xlim <- params[["xlim"]]
  segs <- params[["segs"]]
  cex <- params[["cex"]]
  accent <- params[["accent"]]
  xlim <- xlimTagDensity(seqinfo(object)[chromosome(object)[1], ], xlim, 0.1)
  xscale <- range(c(start(rowRanges(pview)), end(rowRanges(pview))))
  vps2 <- tagdensity_filter_viewports(xscale, ylim=c(-3,2))
  start(xlim) <- xscale[1]
  end(xlim) <- xscale[2]  
  upViewport(0)
  pushViewport(vpdata)
  vp <- vps[["tagdens"]]
  pushViewport(vp)
  .draw_filter_tracks(object,
                      view=pview,
                      filterList=filterList,
                      zoom.out=zoom.out,
                      xlim=xlim,
                      segs=segs, 
                      accent=accent,
                      vps=vps2,
                      params=params)  
}

draw_filters <- function(object, vps, params, pview, subset=TRUE){
  if(missing(vps)) vps <- Viewports()
  group.index <- params[["group.index"]][[1]]
  ##aview <- params[["aview"]]
  object_all <- object
  object <- object[group.index]
  iparams <- .ideogramParams(object, params)
  idiogram <- grid.idiogram(iparams)
  vpdata <- viewport(x=0, y=0, width=unit(0.8, "npc"), height=unit(1, "npc"),
                     just=c("left", "bottom"))
  label <- paste(colnames(pview), chromosome(object)[1], sep="\n")
  vplist <- list()
  ##
  ## Begin drawing figure
  ##
  xlim_tagd <- xlimTagDensity(seqinfo(object)[chromosome(object)[1], ],
                              params[["xlim"]], 0.1)
  if(subset){
    params$pview <- .subset_pview(object, pview, xlim_tagd)
  } else{
    params$pview <- keepSeqlevels(pview, chromosome(object)[1])
  }
  pview <- params$pview
  xscale <- range(c(start(rowRanges(pview)), end(rowRanges(pview))))
  vps2 <- tagdensity_filter_viewports(xscale, ylim=c(-3,2))  
  .draw_filter_tracks2(object, params, vpdata, vps)
}

.draw_tagdensity3 <- function(object,
                              view,
                              percent=0.1,
                              filterList,
                              zoom.out=1,
                              xlim,
                              ylim=c(-3, 2),
                              logy=TRUE,
                              segs,
                              accent,
                              vps,
                              params,
                              ...){
  pch <- params[["pch"]]
  col <- params[["pch_color"]]
  cex <- params[["cex"]]
  seg_col <- params[["seg_color"]]
  g <- GRanges(chromosome(object)[1], IRanges(start(xlim), end(xlim)))
  seqlevels(segs, force=TRUE) <- chromosome(object)[1]
  k <- subjectHits(findOverlaps(g, rowRanges(view), maxgap=start(g)-start(xlim)))
  view <- view[k, ]
  xscale <- range(c(start(rowRanges(view)), end(rowRanges(view))))
  td <- assays(view)[,1]
  td <- threshold(td, lim=ylim)
  segs$seg.mean <- threshold(segs$seg.mean, lim=ylim)
  locs <- pretty(range(xlim), n=10)
  labels <- prettyNum(locs/1000, big.mark=",")
  yaxis <- yaxisLabels(ylim, n=6, logscale=TRUE)
  xlim_g <- GRanges(seqnames(g), IRanges(start(xlim), end(xlim)))
  filterList2 <- subsetFilterList(filterList, xlim_g)
  any_filters <- sum(elementLengths(filterList2)) > 0
  cnv <- variant(object)
  pushViewport(vps[["vlayout"]])
  pushViewport(vps[["vpl1"]])
  pushViewport(vps[["tdvp"]])
  grid.rect(gp=gpar(col="gray"))
  grid.points(x=start(rowRanges(view)),
              y=td, pch=pch,
              default.units="native",
              gp=gpar(col=col, cex=cex))
  grid.segments(x0=start(segs),
                x1=end(segs),
                y0=segs$seg.mean,
                y1=segs$seg.mean,
                default.units="native", gp=gpar(lwd=1.5, col=seg_col))
  upViewport()
  pushViewport(vps[["axisvp"]])
  grid.yaxis(gp=gpar(cex=0.6), main=FALSE)
  at <- pretty(xscale, n=8)
  at <- at[at >= xscale[1] & at <= tail(xscale, 1)]
  labels <- prettyNum(at/1000, big.mark=",")
  grid.xaxis(at=at, label=labels, gp=gpar(cex=0.6))
  upViewport()
  ##tg <- textGrob("log2 tag\ndensity", gp=gpar(cex=0.7))
  pushViewport(vps[["axislabel"]])
  grid.rect(gp=gpar(fill="wheat", col=NA))
  grid.text("log2 tag\ndensity", rot=90, gp=gpar(cex=0.7))
  upViewport(2)
}

.draw_tagdensity2 <- function(object, params, vpdata, vps){
  pview <- params[["pview"]]
  filterList <- params[["filterList"]]
  zoom.out <- params[["zoom.out"]]
  xlim <- params[["xlim"]]
  segs <- params[["segs"]]
  cex <- params[["cex"]]
  accent <- params[["accent"]]
  xlim <- xlimTagDensity(seqinfo(object)[chromosome(object)[1], ], xlim, 0.1)
  xscale <- range(c(start(rowRanges(pview)), end(rowRanges(pview))))
  vps2 <- tagdensity_filter_viewports(xscale, ylim=c(-3,2))
  start(xlim) <- xscale[1]
  end(xlim) <- xscale[2]  
  upViewport(0)
  pushViewport(vpdata)
  vp <- vps[["tagdens"]]
  pushViewport(vp)
  .draw_tagdensity3(object,
                    view=pview,
                    filterList=filterList,
                    zoom.out=zoom.out,
                    xlim=xlim,
                    segs=segs, 
                    accent=accent,
                    vps=vps2,
                    params=params)
}

draw_logratios2 <- function(object, vps, params, pview, subset=TRUE){
  if(missing(vps)) vps <- Viewports()
  group.index <- params[["group.index"]][[1]]
  object_all <- object
  object <- object[group.index]
  iparams <- .ideogramParams(object, params)
  idiogram <- grid.idiogram(iparams)
  vpdata <- viewport(x=0, y=0, width=unit(0.8, "npc"), height=unit(1, "npc"),
                     just=c("left", "bottom"))
  label <- paste(colnames(pview), chromosome(object)[1], sep="\n")
  vplist <- list()
  ##
  ## Begin drawing figure
  ##
  xlim_tagd <- xlimTagDensity(seqinfo(object)[chromosome(object)[1], ],
                              params[["xlim"]], 0.1)
  if(subset){
    params$pview <- .subset_pview(object, pview, xlim_tagd)
  } else{
    params$pview <- keepSeqlevels(pview, chromosome(object)[1])
  }
  pview <- params$pview
  xscale <- range(c(start(rowRanges(pview)), end(rowRanges(pview))))
  vps2 <- tagdensity_filter_viewports(xscale, ylim=c(-3,2))  
  .draw_tagdensity2(object, params, vpdata, vps)
}

draw_connection <- function(object, vps, params, pview, subset=TRUE){
  if(missing(vps)) vps <- Viewports()
  group.index <- params[["group.index"]][[1]]
  object_all <- object
  object <- object[group.index]
  iparams <- .ideogramParams(object, params)
  idiogram <- grid.idiogram(iparams)
  vpdata <- viewport(x=0, y=0, width=unit(0.8, "npc"), height=unit(1, "npc"),
                     just=c("left", "bottom"))
  label <- paste(colnames(pview), chromosome(object)[1], sep="\n")
  vplist <- list()
  xlim_tagd <- xlimTagDensity(seqinfo(object)[chromosome(object)[1], ],
                              params[["xlim"]], 0.1)
  if(subset){
    params$pview <- .subset_pview(object, pview, xlim_tagd)
  } else{
    params$pview <- keepSeqlevels(pview, chromosome(object)[1])
  }  
  pview <- params[["pview"]]
  xscale <- range(c(start(rowRanges(pview)), end(rowRanges(pview))))
  vps2 <- tagdensity_filter_viewports(xscale, ylim=c(-3,2))  
  upViewport(0)
  .connect_ideogram(vpdata, vps, vps2, iparams, xscale)
}
