#' @include AllGenerics.R
NULL

IdiogramParams <- function(seqnames=character(),
                           seqlengths=numeric(), unit="kb", genome="hg19",
                           box=list(xlim=ILimit(), color="blue")){
  new("IdiogramParams", seqnames=seqnames, seqlengths=seqlengths, unit=unit,
      genome=genome, box=box)
}

setMethod("boxIdiogram", "IdiogramParams", function(object) object@box)
setMethod("chromosome", "IdiogramParams", function(object) object@seqnames)
setMethod("genome", "IdiogramParams", function(x) x@genome)
setMethod("seqlengths", "IdiogramParams", function(x) x@seqlengths)
setMethod("seqnames", "IdiogramParams", function(x) x@seqnames)

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
  seqlevels(g, force=TRUE) <- chromosome(g)
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

#' @export
ggIdeogram <- function(chroms, fusion.dat, g.params){
  ip1 <- IdiogramParams(chroms[1])
  coords1 <- idiogram_coordinates(ip1)
  if(length(unique(chroms)) > 1){
    ip2 <- IdiogramParams(chroms[2])
    coords2 <- idiogram_coordinates(ip2)
    coords1$seqnames <- chroms[1]
    coords2$seqnames <- chroms[2]
    coords <- rbind(coords1, coords2)
    coords$seqnames <- factor(coords$seqnames, chroms)
    ##load_all("integration/integration")
    ideo <- .gg_ideogram_twochrom(coords, fusion.dat, g.params[[1]],
                        g.params[[2]])
    g.ideo <- ggplotGrob(ideo)
    ##g.ideo$heights[1] <- unit(1, "pt")
  } else {
    ideo <- .gg_ideogram(coords1, fusion.dat, g.params[[1]],
                       g.params[[2]])
    g.ideo <- ggplotGrob(ideo)
  }
  g.ideo
}

##ggIdeogram <- function(chroms, fusion.dat, g.params){
##  ip1 <- svplots:::IdiogramParams(chroms[1])
##  coords1 <- svplots:::idiogram_coordinates(ip1)
##  if(length(unique(chroms)) > 1){
##    ip2 <- svplots:::IdiogramParams(chroms[2])
##    coords2 <- svplots:::idiogram_coordinates(ip2)
##    coords1$seqnames <- chroms[1]
##    coords2$seqnames <- chroms[2]
##    coords <- rbind(coords1, coords2)
##    coords$seqnames <- factor(coords$seqnames, chroms)
##    ##load_all("integration/integration")
##    ideo <- .gg_ideogram_twochrom(coords, fusion.dat, g.params[[1]],
##                        g.params[[2]])
##    g.ideo <- ggplotGrob(ideo)
##    ##g.ideo$heights[1] <- unit(1, "pt")
##  } else {
##    ideo <- .gg_ideogram(coords1, fusion.dat, g.params[[1]],
##                       g.params[[2]])
##    g.ideo <- ggplotGrob(ideo)
##  }
##  g.ideo
##}

## gg_ideogram <- function(coords, fusion.dat, gene1.params, gene2.params){
##  colors <- getOption("biovizBase")$cytobandColor
##  gene1 <- gene1.params$gene.name
##  gene2 <- gene2.params$gene.name
## 
##  ustains <- unique(coords$gieStain)
##  colors <- colors[ustains]
##  coords$gieStain <- factor(coords$gieStain, levels=ustains)
##  coords$ymin <- -0.5
##  coords$ymax <- 0.5
## 
##  taper.coords <- taperIdeogram(coords)
##  coords2 <- coords[-c(1, nrow(coords)), ]
##  coords3 <- coords2[coords2$gieStain != "acen", ]
## 
##  genes <- strsplit(fusion.dat$fusion, "-")[[1]]
##  tx1 <- fusion.dat$exons[[gene1]]
##  tx2 <- fusion.dat$exons[[gene2]]
##  chroms <- as.character(c(tx1$seqnames[1], tx2$seqnames[2]))
##  stopifnot(identical(chroms[1], chroms[2]))
##  gene.coords <- data.frame(chrom=chroms,
##                            start=c(min(tx1$start),
##                                    min(tx2$start)),
##                            end=c(max(tx1$end),
##                                  max(tx2$end)),
##                            gene=genes)
##  gene.coords$gene <- factor(gene.coords$gene, levels=genes)
##  coords3$chrom <- as.character(seqnames(fusion.dat[[1]]))[1]
##  gene.colors <- c(gene1.params$exon.color, gene2.params$exon.color)
## 
##  ggplot(coords3, aes(xmin=xleft, xmax=xright, ymin=ymin, ymax=ymax,
##                           fill=gieStain, color=gieStain)) +
##    geom_polygon(data=subset(taper.coords, id=="p.telomere"),
##                 aes(x, y, group=id),
##                 fill=colors["gneg"],
##                 inherit.aes=FALSE, color="black") +
##    geom_polygon(data=subset(taper.coords, id %in% c("p.centromere", "q.centromere")),
##                 aes(x, y, group=id),
##                 fill=colors["acen"],
##                 inherit.aes=FALSE, color="black") +
##    geom_polygon(data=subset(taper.coords, id=="q.telomere"),
##                 aes(x, y, group=id),
##                 fill=colors["gneg"],
##                 inherit.aes=FALSE, color="black") +
##    geom_rect(color="black") +
##    geom_rect(data=subset(gene.coords, gene==genes[2]),
##              aes(xmin=start, xmax=end,
##                  ymin=-Inf, ymax=+Inf),
##              color=gene.colors[2],
##              inherit.aes=FALSE, fill=gene.colors[2], alpha=0.3) +
##    geom_rect(data=subset(gene.coords, gene==genes[1]),
##              aes(xmin=start, xmax=end,
##                  ymin=-Inf, ymax=+Inf),
##              color=gene.colors[1],
##              inherit.aes=FALSE,
##              fill=gene.colors[1], alpha=0.3) +
##    scale_fill_manual(values=colors) +
##    scale_color_manual(values=colors) +
##    scale_x_continuous(expand=c(0, 0)) +
##    scale_y_continuous(breaks=0, labels=coords3$chrom[1], expand=c(0.5, 0)) + 
##    theme(legend.position="none",
##          axis.text.x=element_blank(),
##          axis.text.y=element_text(angle=0, size=15), 
##          axis.line=element_blank(),
##          axis.ticks=element_blank(),
##          axis.title=element_blank(),
##          panel.background=element_rect(fill="transparent"))
## }

#' @export
ideogramData <- function(chrom, param.list, build="hg19"){
  gene1.params <- param.list[[1]]
  gene2.params <- param.list[[2]]
  ip <- IdiogramParams(chrom)
  coords <- idiogram_coordinates(ip)
  colors <- getOption("biovizBase")$cytobandColor
  gene1 <- gene1.params$gene.name
  gene2 <- gene2.params$gene.name

  ustains <- unique(coords$gieStain)
  colors <- colors[ustains]
  coords$gieStain <- factor(coords$gieStain, levels=ustains)
  coords$ymin <- -0.5
  coords$ymax <- 0.5

  taper.coords <- taperIdeogram(coords) %>%
    as.tibble
  coords2 <- coords[-c(1, nrow(coords)), ]
  coords3 <- coords2[coords2$gieStain != "acen", ] %>%
    as.tibble
  genes <- c(gene1.params$gene.name,
             gene2.params$gene.name)
  tx <- loadTx(build)
  tx <- tx[tx$gene_name %in% genes]
  tx <- tx[!duplicated(tx$gene_name)]
  ##  genes <- strsplit(fusion.dat$fusion, "-")[[1]] %>%
  ##    as.data.frame %>%
  ##    mutate(chrom=chromosome(tx))
  gene.coords <- as.data.frame(tx) %>%
    mutate(chrom=chromosome(tx),
           ymin=-0.75,
           ymax=0.75)
  genes <- gene.coords$gene_name
  coords3$chrom <- gene.coords$chrom
  ##gene.coords$gene <- factor(gene.coords$gene, levels=genes)
  ##coords3$chrom <- as.character(seqnames(fusion.dat[[1]]))[1]
  gene.coords$gene.colors <- c(gene1.params$exon.color,
                               gene2.params$exon.color)
  list(gene.coords=gene.coords,
       taper.coords=taper.coords,
       ideogram.coords=coords3,
       colors=colors)
}

#' @export
ggIdeogram2 <- function(dat.list){
  coords3 <- dat.list[["ideogram.coords"]]
  gene.coords <- dat.list[["gene.coords"]]
  taper.coords <- dat.list[["taper.coords"]]
  colors <- dat.list[["colors"]]
  gene.colors <- gene.coords$gene.colors
  genes <- gene.coords$gene_name
  ##ggplot(coords3, aes(xmin=xleft, xmax=xright, ymin=ymin, ymax=ymax,
  ggplot(coords3,
         aes(xmin=xleft, xmax=xright, ymin=ymin, ymax=ymax,
             fill=gieStain, color=gieStain)) +
    geom_polygon(data=subset(taper.coords, id=="p.telomere"),
                 aes(x, y, group=id),
                 fill=colors["gneg"],
                 inherit.aes=FALSE, color="black") +
    geom_polygon(data=subset(taper.coords,
                             id %in% c("p.centromere", "q.centromere")),
                 aes(x, y, group=id),
                 fill=colors["acen"],
                 inherit.aes=FALSE, color="black") +
    geom_polygon(data=subset(taper.coords, id=="q.telomere"),
                 aes(x, y, group=id),
                 fill=colors["gneg"],
                 inherit.aes=FALSE, color="black") +
    geom_rect(color="black") +
    geom_rect(data=subset(gene.coords, gene_name==genes[2]),
              aes(xmin=start, xmax=end,
                  ymin=ymin, ymax=ymax),
              color=gene.colors[2],
              inherit.aes=FALSE, fill=gene.colors[2], alpha=0.3) +
    geom_rect(data=subset(gene.coords, gene_name==genes[1]),
              aes(xmin=start, xmax=end,
                  ymin=ymin, ymax=ymax),
              color=gene.colors[1],
              inherit.aes=FALSE,
              fill=gene.colors[1], alpha=0.3) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(breaks=0, labels=coords3$chrom[1],
                       expand=c(0.5, 0)) +
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.text.y=element_text(angle=0, size=15), 
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          panel.background=element_rect(fill="transparent"))
}


.gg_ideogram <- function(coords, fusion.dat,
                         gene1.params, gene2.params){
  colors <- getOption("biovizBase")$cytobandColor
  gene1 <- gene1.params$gene.name
  gene2 <- gene2.params$gene.name

  ustains <- unique(coords$gieStain)
  colors <- colors[ustains]
  coords$gieStain <- factor(coords$gieStain, levels=ustains)
  coords$ymin <- -0.5
  coords$ymax <- 0.5

  taper.coords <- taperIdeogram(coords)
  coords2 <- coords[-c(1, nrow(coords)), ]
  coords3 <- coords2[coords2$gieStain != "acen", ]

  genes <- strsplit(fusion.dat$fusion, "-")[[1]]
  tx1 <- fusion.dat$exons[[gene1]]
  tx2 <- fusion.dat$exons[[gene2]]
  chroms <- as.character(c(tx1$seqnames[1], tx2$seqnames[1]))
  stopifnot(identical(chroms[1], chroms[2]))
  gene.coords <- data.frame(chrom=chroms,
                            start=c(min(tx1$start),
                                    min(tx2$start)),
                            end=c(max(tx1$end),
                                  max(tx2$end)),
                            gene=genes)
  gene.coords$gene <- factor(gene.coords$gene, levels=genes)
  coords3$chrom <- as.character(seqnames(fusion.dat[[1]]))[1]
  gene.colors <- c(gene1.params$exon.color, gene2.params$exon.color)

  ggplot(coords3, aes(xmin=xleft, xmax=xright, ymin=ymin, ymax=ymax,
                           fill=gieStain, color=gieStain)) +
    geom_polygon(data=subset(taper.coords, id=="p.telomere"),
                 aes(x, y, group=id),
                 fill=colors["gneg"],
                 inherit.aes=FALSE, color="black") +
    geom_polygon(data=subset(taper.coords, id %in% c("p.centromere", "q.centromere")),
                 aes(x, y, group=id),
                 fill=colors["acen"],
                 inherit.aes=FALSE, color="black") +
    geom_polygon(data=subset(taper.coords, id=="q.telomere"),
                 aes(x, y, group=id),
                 fill=colors["gneg"],
                 inherit.aes=FALSE, color="black") +
    geom_rect(color="black") +
    geom_rect(data=subset(gene.coords, gene==genes[2]),
              aes(xmin=start, xmax=end,
                  ymin=-Inf, ymax=+Inf),
              color=gene.colors[2],
              inherit.aes=FALSE, fill=gene.colors[2], alpha=0.3) +
    geom_rect(data=subset(gene.coords, gene==genes[1]),
              aes(xmin=start, xmax=end,
                  ymin=-Inf, ymax=+Inf),
              color=gene.colors[1],
              inherit.aes=FALSE,
              fill=gene.colors[1], alpha=0.3) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(breaks=0, labels=coords3$chrom[1], expand=c(0.5, 0)) + 
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.text.y=element_text(angle=0, size=15), 
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          panel.background=element_rect(fill="transparent"))
}

## .gg_ideogram_twochrom <- function(coords,
##                                   fusion.dat, gene1.params,
##                                   gene2.params){
##   colors <- getOption("biovizBase")$cytobandColor
##   gene1 <- gene1.params$gene.name
##   gene2 <- gene2.params$gene.name
## 
##   ustains <- unique(coords$gieStain)
##   colors <- colors[ustains]
##   coords$gieStain <- factor(coords$gieStain, levels=ustains)
##   coords$ymin <- -0.5
##   coords$ymax <- 0.5
## 
##   
##   genes <- strsplit(fusion.dat$fusion, "-")[[1]]
##   tx1 <- fusion.dat$exons[[gene1]]
##   tx2 <- fusion.dat$exons[[gene2]]
## 
##   chroms <- levels(coords$seqnames)
##   coords1 <- coords[coords$seqnames==chroms[1], ]
##   coords2 <- coords[coords$seqnames==chroms[2], ]
##   taper.coords1 <- taperIdeogram(coords1)
##   taper.coords2 <- taperIdeogram(coords2)
##   taper.coords1$seqnames <- chroms[1]
##   taper.coords2$seqnames <- chroms[2]
##   taper.coords <- rbind(taper.coords1, taper.coords2)
##   taper.coords$seqnames <- factor(taper.coords$seqnames, levels=chroms)
##   coords2 <- coords[-c(1, nrow(coords)), ]
##   coords3 <- coords2[coords2$gieStain != "acen", ]
## 
##   chroms <- as.character(c(tx1$seqnames[1], tx2$seqnames[2]))
##   gene.coords <- data.frame(chrom=chroms,
##                             start=c(min(tx1$start),
##                                     min(tx2$start)),
##                             end=c(max(tx1$end),
##                                   max(tx2$end)),
##                             gene=genes)
##   gene.coords$gene <- factor(gene.coords$gene, levels=genes)
##   gene.coords$seqnames <- factor(levels(coords$seqnames),
##                                  levels=levels(coords$seqnames))
##   coords3$chrom <- as.character(seqnames(fusion.dat[[1]]))[1]
##   gene.colors <- c(gene1.params$exon.color, gene2.params$exon.color)
## 
##   ggplot(coords3, aes(xmin=xleft, xmax=xright, ymin=ymin, ymax=ymax,
##                            fill=gieStain, color=gieStain)) +
##     geom_polygon(data=subset(taper.coords, id=="p.telomere"),
##                  aes(x, y, group=id),
##                  fill=colors["gneg"],
##                  inherit.aes=FALSE, color="black") +
##     geom_polygon(data=subset(taper.coords, id %in% c("p.centromere", "q.centromere")),
##                  aes(x, y, group=id),
##                  fill=colors["acen"],
##                  inherit.aes=FALSE, color="black") +
##     geom_polygon(data=subset(taper.coords, id=="q.telomere"),
##                  aes(x, y, group=id),
##                  fill=colors["gneg"],
##                  inherit.aes=FALSE, color="black") +
##     geom_rect(color="black") +
##     geom_rect(data=subset(gene.coords, gene==genes[2]),
##               aes(xmin=start, xmax=end,
##                   ymin=-Inf, ymax=+Inf),
##               color=gene.colors[2],
##               inherit.aes=FALSE, fill=gene.colors[2], alpha=0.3) +
##     geom_rect(data=subset(gene.coords, gene==genes[1]),
##               aes(xmin=start, xmax=end,
##                   ymin=-Inf, ymax=+Inf),
##               color=gene.colors[1],
##               inherit.aes=FALSE,
##               fill=gene.colors[1], alpha=0.3) +
##     scale_fill_manual(values=colors) +
##     scale_color_manual(values=colors) +
##     scale_x_continuous(expand=c(0, 0)) +
##     theme(legend.position="none",
##           axis.text.x=element_blank(),
##           axis.text.y=element_blank(), 
##           axis.line=element_blank(),
##           axis.ticks=element_blank(),
##           axis.title=element_blank(),
##           panel.background=element_rect(fill="transparent"),
##           strip.text=element_text(lineheight=0.05)) +
##     facet_wrap(~seqnames, scales="free_x")
## }

taperIdeogram <- function(coords, b=0.5){
  t <- seq(0, 2*pi, length.out=100)
  a <- coords$xright[1]
  x <- a*cos(t)
  y <- b*sin(t)
  x <- x - min(x)
  p.telomere <- data.frame(x=x, y=y, id="p.telomere",
                           gieStain=coords$gieStain[1])

  a <- coords$xright[nrow(coords)] - coords$xleft[nrow(coords)]
  x <- a*cos(t)
  y <- b*sin(t)
  x <- x + coords$xleft[nrow(coords)]
  q.telomere <- data.frame(x=x, y=y, id="q.telomere",
                           gieStain=coords$gieStain[nrow(coords)])

  centrom <- coords[coords$gieStain=="acen", ]
  a <- centrom$xright[1] - centrom$xleft[1] 
  x <- a*cos(t)
  y <- b*sin(t)
  x <- x + centrom$xleft[1]
  p.centrom <- data.frame(x=x, y=y, id="p.centromere",
                          gieStain="acen")
  a <- centrom$xright[2] - centrom$xleft[2] 
  x <- a*cos(t)
  y <- b*sin(t)
  x <- x + centrom$xright[2]
  q.centrom <- data.frame(x=x, y=y, id="q.centromere",
                          gieStain="acen")

  rbind(p.telomere, p.centrom, q.centrom, q.telomere)
}

.gg_ideogram_twochrom <- function(coords,
                                  fusion.dat, gene1.params,
                                  gene2.params){
  colors <- getOption("biovizBase")$cytobandColor
  gene1 <- gene1.params$gene.name
  gene2 <- gene2.params$gene.name

  ustains <- unique(coords$gieStain)
  colors <- colors[ustains]
  coords$gieStain <- factor(coords$gieStain, levels=ustains)
  coords$ymin <- -0.5
  coords$ymax <- 0.5

  
  genes <- strsplit(fusion.dat$fusion, "-")[[1]]
  tx1 <- fusion.dat$exons[[gene1]]
  tx2 <- fusion.dat$exons[[gene2]]

  chroms <- levels(coords$seqnames)
  coords1 <- coords[coords$seqnames==chroms[1], ]
  coords2 <- coords[coords$seqnames==chroms[2], ]
  taper.coords1 <- taperIdeogram(coords1)
  taper.coords2 <- taperIdeogram(coords2)
  taper.coords1$seqnames <- chroms[1]
  taper.coords2$seqnames <- chroms[2]
  taper.coords <- rbind(taper.coords1, taper.coords2)
  taper.coords$seqnames <- factor(taper.coords$seqnames, levels=chroms)
  coords2 <- coords[-c(1, nrow(coords)), ]
  coords3 <- coords2[coords2$gieStain != "acen", ]

  chroms <- as.character(c(tx1$seqnames[1], tx2$seqnames[2]))
  gene.coords <- data.frame(chrom=chroms,
                            start=c(min(tx1$start),
                                    min(tx2$start)),
                            end=c(max(tx1$end),
                                  max(tx2$end)),
                            gene=genes)
  gene.coords$gene <- factor(gene.coords$gene, levels=genes)
  gene.coords$seqnames <- factor(levels(coords$seqnames),
                                 levels=levels(coords$seqnames))
  coords3$chrom <- as.character(seqnames(fusion.dat[[1]]))[1]
  gene.colors <- c(gene1.params$exon.color, gene2.params$exon.color)

  ggplot(coords3, aes(xmin=xleft, xmax=xright, ymin=ymin, ymax=ymax,
                           fill=gieStain, color=gieStain)) +
    geom_polygon(data=subset(taper.coords, id=="p.telomere"),
                 aes(x, y, group=id),
                 fill=colors["gneg"],
                 inherit.aes=FALSE, color="black") +
    geom_polygon(data=subset(taper.coords, id %in% c("p.centromere", "q.centromere")),
                 aes(x, y, group=id),
                 fill=colors["acen"],
                 inherit.aes=FALSE, color="black") +
    geom_polygon(data=subset(taper.coords, id=="q.telomere"),
                 aes(x, y, group=id),
                 fill=colors["gneg"],
                 inherit.aes=FALSE, color="black") +
    geom_rect(color="black") +
    geom_rect(data=subset(gene.coords, gene==genes[2]),
              aes(xmin=start, xmax=end,
                  ymin=-Inf, ymax=+Inf),
              color=gene.colors[2],
              inherit.aes=FALSE, fill=gene.colors[2], alpha=0.3) +
    geom_rect(data=subset(gene.coords, gene==genes[1]),
              aes(xmin=start, xmax=end,
                  ymin=-Inf, ymax=+Inf),
              color=gene.colors[1],
              inherit.aes=FALSE,
              fill=gene.colors[1], alpha=0.3) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    scale_x_continuous(expand=c(0, 0)) +
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.text.y=element_blank(), 
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          panel.background=element_rect(fill="transparent"),
          strip.text=element_text(lineheight=0.05)) +
    facet_wrap(~seqnames, scales="free_x")
}
