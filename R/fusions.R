.gg_clipped_tx1 <- function(data.list, params, roi, clip=FALSE){
  manual.colors <- c(params[["exon.color"]], params[["clipped.color"]])
  ##rreads <- data.list[["rearranged.reads"]]
  gene <- params[["gene.name"]]
  ##gr <- rreads[gene]
  gr <- roi[gene]
  ##y.lab <- params[["y.axis.label"]]
  exons <- data.list$exons[[gene]]
  strand <- params$rearrangedStrand
  is.minus <- strand=="-"
  exons$is_clipped <- factor(exons$is_clipped, levels=c("FALSE", "TRUE"))
  if(clip){
    exons <- exons[exons$is_clipped != TRUE]
  }
  manual.colors <- unique(manual.colors[sort(as.integer(exons$is_clipped))])
  if(!"rank" %in% colnames(exons)){
    exons$rank <- exons$exon_rank
  }
  ranks <- exons$rank
  ranks <- ranks[ranks != ""]
  if(is.minus) xmin <- +Inf else xmin <- -Inf
  clip.at <- gr$bp.jxn
  p <- ggplot(subset(exons, rank %in% ranks),
         aes(xmin=start, xmax=end, ymin=-0.5, ymax=0.5,
             label=rank,
             color=is_clipped,
             fill=is_clipped)) +
    geom_rect(fill="transparent", color="transparent") + ## set axes
    geom_rect(data=exons, aes(xmin=xmin,
                              xmax=clip.at,
                              ymin=-Inf, ymax=+Inf),
              fill=params[["background.color"]],
              color=params[["background.color"]],
              alpha=0.2, inherit.aes=FALSE) +
    geom_rect(data=exons,
              aes(xmin=clip.at,
                  xmax=-xmin,
                  ymin=-Inf, ymax=+Inf),
              fill="gray",
              color="transparent", alpha=0.9, inherit.aes=FALSE) +
    geom_hline(yintercept=0, color="black") +
    geom_rect(size=1) +
    scale_y_continuous(breaks=0, labels="5'", expand=c(0, 0), limits=c(-0.5, 0.5)) +
    scale_x_continuous(expand=c(0,0)) +
    guides(fill="none", color="none") +
    geom_vline(xintercept=c(clip.at-20, clip.at + 20), linetype="dashed") +
    scale_color_manual(values=manual.colors) +
    scale_fill_manual(values=manual.colors) +
    theme(axis.text.x=element_blank(),
          axis.title=element_blank(),
          axis.text.y=element_text(size=12),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          panel.background=element_rect(fill=params[["background.color"]]),
          panel.grid.minor.x=element_line(color="transparent"),
          panel.border=element_rect(color="transparent", fill=NA),
          ## top, right, bottom, left
          ##plot.margin=unit(c(1, -0.1, 1, 1), "lines")) +
          plot.margin=unit(c(1, 0, 1, 0), "lines")) +
    ggtitle(params[["gene.name"]])
  if(is.minus) p <- p + scale_x_reverse(expand=c(0,0))
  p
}


.gg_clipped_tx2 <- function(data.list, params, roi, clip=FALSE){
  manual.colors <- c(params[["exon.color"]], params[["clipped.color"]])
  rreads <- data.list[["rearranged.reads"]]
  gene <- params[["gene.name"]]
  ##gr <- rreads[gene]
  gr <- roi[gene]
  exons <- data.list$exons[[gene]]
  strand <- params$strand
  is.minus <- strand=="-"
  exons$is_clipped <- factor(exons$is_clipped, levels=c("FALSE", "TRUE"))
  if(clip){
    exons <- exons[exons$is_clipped != TRUE]
  }
  manual.colors <- unique(manual.colors[sort(as.integer(exons$is_clipped))])
  if(!"rank" %in% colnames(exons)){
    exons$rank <- exons$exon_rank
  }
  ranks <- exons$rank
  ranks <- ranks[ranks != ""]
  if(is.minus) xmin <- -Inf else xmin <- +Inf
  clip.at <- gr$bp.jxn
  ## since we're flipping the axis if on minus strand, the y-axis label is always 5'
  p <- ggplot(subset(exons, rank %in% ranks),
         aes(xmin=start, xmax=end, ymin=-0.5, ymax=0.5,
             label=rank,
             color=is_clipped,
             fill=is_clipped)) +
    geom_rect(fill="transparent", color="transparent") + ## set axes
    geom_rect(data=exons, aes(xmin=xmin,
                              xmax=clip.at,
                              ymin=-Inf, ymax=+Inf),
              fill=params[["background.color"]],
              color=params[["background.color"]],
              alpha=0.2, inherit.aes=FALSE) +
    geom_rect(data=exons,
              aes(xmin=clip.at,
                  xmax=-xmin,
                  ymin=-Inf, ymax=+Inf),
              fill="gray",
              color="transparent", alpha=0.9, inherit.aes=FALSE) +
    geom_hline(yintercept=0, color="black") +
    geom_rect(size=1) +
    scale_y_continuous(breaks=0, labels="5'", expand=c(0, 0), limits=c(-0.5, 0.5)) +
    scale_x_continuous(expand=c(0,0)) +
    guides(fill="none", color="none") +
    geom_vline(xintercept=clip.at, linetype="dashed") +
    scale_color_manual(values=manual.colors) +
    scale_fill_manual(values=manual.colors) +
    theme(axis.text.x=element_blank(),
          axis.title=element_blank(),
          axis.text.y=element_text(size=12),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          panel.background=element_rect(fill=params[["background.color"]]),
          panel.grid.minor.x=element_line(color="transparent"),
          panel.border=element_rect(color="transparent", fill=NA),
          ## top, right, bottom, left
          plot.margin=unit(c(1, 0, 1, 0), "lines")) +
    ggtitle(params[["gene.name"]])
  if(is.minus) p <- p + scale_x_reverse(expand=c(0,0))
  p
}

#' @param data.list a list.  see readFusionData
ggClippedExons <- function(data.list, params, roi, clip=FALSE){
  is.first <- params$is.first
  if(is.first){
    p <- .gg_clipped_tx1(data.list, params, roi, clip=clip)
  } else{
    p <- .gg_clipped_tx2(data.list, params, roi, clip=clip)
  }
  p
}

#'@include ideogram.R
NULL

ggAxisLabel <- function(data.list, params, label){
  gene <- params$gene.name
  exons <- data.list$exons[[gene]]
  ranks <- exons$exon_rank
  ranks <- ranks[ranks!=""]
  ggplot(subset(exons, exon_rank %in% ranks),
         aes(xmin=start, xmax=end, ymin=-0.5, ymax=0.5, label=exon_rank)) +
    geom_rect(size=1, color="transparent", fill="transparent") +
    geom_text(aes(x=midx, y=0.6), color="transparent", size=3) +
    scale_y_continuous(breaks=0.05, labels=label) +
    guides(fill="none", color="none") +
    theme(axis.text.x=element_blank(),
          axis.title=element_blank(),
          axis.text.y=element_text(size=10),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          panel.background=element_rect(fill="transparent"),
          panel.grid.minor.x=element_line(color="transparent"),
          panel.border=element_rect(color="transparent", fill=NA),
          ## top, right, bottom, left
          plot.margin=unit(c(1, 0, 1, 0), "lines")) +
    ggtitle("")
}

ggTxFusion <- function(dat, gene1.params, gene2.params){
  strand <- as.character(dat$exons[["strand"]][1])
  if(strand=="+"){
    p <- gg_fusion_plus(dat, gene1.params, gene2.params)
  }
  if(strand=="-"){
    p <- gg_fusion_minus(dat, gene1.params, gene2.params)
  }
  ##annotate("text", max(exons$end2)+5000, 1, label="3'", size=4) ##+
  ##facet_wrap(~tumor, switch="y")
  p
}


gg_fusion_minus <- function(dat, gene1.params, gene2.params){
  exons <- dat$exons
  rreads <- dat$rearranged.reads
  gene.colors <- c(gene1.params$exon.color, gene2.params$exon.color)
  back.colors <- c(gene1.params$background.color,
                   gene2.params$background.color)
  gene1 <- gene1.params$gene.name
  gene2 <- gene2.params$gene.name
  sequence_jxn <- exons$sequence_junction[1]
  jxn.label <- paste(unique(exons$sequence_junction), collapse="-")
  exons$gene <- factor(exons$gene, levels=c(gene1, gene2))
  exons$ymin <- -0.5
  exons$ymax <- 0.5
  ggplot(exons,
         aes(xmin=start, xmax=end, ymin=ymin,
             ymax=ymax, color=gene,
             fill=gene, label=exon_rank)) +
    geom_rect(fill="transparent", color="transparent") + ## just to set axes
    ## background for first gene
    geom_rect(data=exons, aes(xmin=sequence_jxn,
                              xmax=+Inf,
                              ymin=-Inf, ymax=+Inf),
              fill=back.colors[1],
              color=back.colors[1], alpha=0.2, inherit.aes=FALSE) +
    ## background for second gene
    geom_rect(aes(xmin=-Inf, xmax=sequence_jxn,
                  ymin=-Inf, ymax=+Inf),
              fill=back.colors[2],
              color=back.colors[2],
              alpha=0.2, inherit.aes=FALSE) +
    geom_vline(xintercept=sequence_jxn, color="gray20", linetype="dashed") +
    geom_segment(data=exons,
                 aes(x=min(start), xend=max(end), y=0, yend=0),
                 color="black") +
    geom_rect(size=1) +
    scale_y_continuous(expand=c(0, 0), limits=c(-0.5, 0.5), breaks=0, label="5'") +
    scale_x_reverse(expand=c(0, 0), breaks=sequence_jxn, labels=jxn.label) +
    ##geom_text(aes(x=midx, y=1.3), color="black", size=3) +
    ylab("") + xlab("") +
    guides(fill="none", color="none") +
    scale_color_manual(values=gene.colors) +
    scale_fill_manual(values=gene.colors) +
    theme(axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=12),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          panel.background=element_rect(fill="transparent"),
          panel.grid.minor.x=element_line(color="transparent"),
          panel.border=element_rect(color="transparent", fill=NA),
          ## top, right, bottom, left
          plot.margin=unit(c(0.5, 1, 0, 1), "mm")) #+
}

gg_fusion_plus <- function(dat, gene1.params, gene2.params){
  exons <- dat$exons
  rreads <- dat$rearranged.reads
  gene.colors <- c(gene1.params$exon.color, gene2.params$exon.color)
  back.colors <- c(gene1.params$background.color,
                   gene2.params$background.color)
  gene1 <- gene1.params$gene.name
  gene2 <- gene2.params$gene.name
  sequence_jxn <- exons$sequence_junction[1]
  jxn.label <- paste(unique(exons$sequence_junction), collapse="-")
  exons$ymin <- -0.5
  exons$ymax <- 0.5
  exons$gene <- factor(exons$gene, levels=c(gene1, gene2))
  ggplot(exons,
         aes(xmin=start, xmax=end, ymin=ymin,
             ymax=ymax, color=gene,
             fill=gene, label=exon_rank)) +
    geom_rect(fill="transparent", color="transparent") + ## just to set axes
    ## background for second gene
    geom_rect(aes(xmin=sequence_jxn, xmax=+Inf,
                  ymin=-Inf, ymax=+Inf),
              fill=back.colors[2],
              color=back.colors[2],
              alpha=0.2, inherit.aes=FALSE) +
    ## background for first gene
    geom_rect(data=exons, aes(xmin=-Inf,
                            xmax=sequence_jxn-1,
                            ymin=-Inf, ymax=+Inf),
              fill=back.colors[1],
              color=back.colors[1], alpha=0.2, inherit.aes=FALSE) +
    geom_vline(xintercept=sequence_jxn, color="gray20", linetype="dashed") +
    geom_segment(data=exons,
                 aes(x=min(start), xend=max(end), y=0, yend=0),
                 color="black") +
    geom_rect(size=1) +
    scale_y_continuous(limits=c(-0.5, 0.5),
                       expand=c(0, 0), breaks=0, label="5'") +
    scale_x_continuous(breaks=sequence_jxn, labels=jxn.label, expand=c(0,0)) +
    ylab("") + xlab("") +
    guides(fill="none", color="none") +
    scale_color_manual(values=gene.colors) +
    scale_fill_manual(values=gene.colors) +
    theme(axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=12),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          panel.background=element_rect(fill="transparent"),
          panel.grid.minor.x=element_line(color="transparent"),
          panel.border=element_rect(color="transparent", fill=NA),
          ## top, right, bottom, left
          plot.margin=unit(c(0.5, 1, 0, 1), "mm")) #+
}

getLegend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


ggRearrangedReads2 <- function(reads.df, g.params, roi, reverse.roi){
  if(reverse.roi) {
    roi <- roi[2:1]
    g.params <- g.params[2:1]
  }
  gene1.params <- g.params[[1]]
  gene2.params <- g.params[[2]]
  genes <- c(gene1.params$gene.name, gene2.params$gene.name)
  gene1 <- genes[1]
  gene2 <- genes[2]
  ##reads.df$type <- paste0(reads.df$read, reads.df$strand)
  ##reads.df$type[reads.df$is.split] <- "split read"
  ##reads.df$type <- factor(reads.df$type, levels=unique(reads.df$type))
  back.colors <- c(gene1.params$background.color, gene2.params$background.color)
  reads.gene1 <- reads.df[reads.df$transcript==gene1, ]
  reads.gene2 <- reads.df[reads.df$transcript==gene2, ]
  reads.gene2 <- reads.gene2[order(reads.gene2$pair.id, decreasing=FALSE), ]
  reads.gene1 <- reads.gene1[order(reads.gene1$pair.id, decreasing=FALSE), ]
  reads.gene1$type <- paste0(reads.gene1$read,
                             reads.gene1$strand, " : ",
                             reads.gene2$read, reads.gene2$strand)
  reads.gene1$type[reads.gene1$is.split] <- "split read"
  reads.gene1$type <- factor(reads.gene1$type, levels=unique(reads.gene1$type))
  reads.gene2$type <- reads.gene1$type
  stopifnot(identical(reads.gene1$pair.id, reads.gene2$pair.id))
  strands <- as.character(strand(roi))
  strand1 <- strands[1]
  strand2 <- strands[2]
  improp <- reads.gene1[!reads.gene1$is.split, ]
  improp <- improp[order(improp$start),  ]
  reads.gene1[!reads.gene1$is.split, ] <- improp
  ##xlim <- c(min(improp$start), roi$bp.jxn[1])
  if(FALSE){
    ## this dosn't work because sometimes we need to reverse the coordinates of only one panel
	  x <- rbind(reads.gene1, reads.gene2)
	  x$transcript <- factor(x$transcript, levels=c(reads.gene1$transcript[1],
	                                                reads.gene2$transcript[2]))
	  ggplot(x,
	         aes(xmin=start, xmax=end, ymin=pair.id-0.3,
	             ymax=pair.id+0.3)) +
	    geom_rect(##data=subset(reads.df, transcript==gene1)[1, ],
        aes(xmin=-Inf, xmax=+Inf, ymin=-Inf, ymax=+Inf),
        fill=back.colors[1], alpha=0.5, inherit.aes=FALSE) +
	    geom_rect(alpha=0.9, aes(fill=type)) +
	    scale_y_continuous(breaks=pretty(unique(x$pair.id), n=4)) +
	    scale_x_continuous(expand=c(0, 0)) +
	    theme(axis.text.x=element_blank(),
	          axis.ticks=element_blank(),
	          axis.text.y=element_blank(),
	          panel.background=element_rect(fill="transparent"),
	          panel.border=element_rect(color="transparent", fill=NA),
	          strip.text=element_blank(),
	          strip.background=element_blank(),
	          legend.text=element_text(size=7),
	          legend.key=element_rect(size=0.5))+
      ## top, right, bottom, left
      ##plot.margin=unit(c(0.1, 0, 0.01, 0.5), "npc")) +
	    ##ylab("rearranged\nread pairs") +
      xlab("") +
	    guides(fill="none") +
	    facet_wrap(~transcript, scale="free_x")
  }
  g1 <- ggplot(reads.gene1,
         aes(xmin=start, xmax=end, ymin=pair.id-0.3,
             ymax=pair.id+0.3)) +
    geom_rect(##data=subset(reads.df, transcript==gene1)[1, ],
              aes(xmin=-Inf, xmax=+Inf, ymin=-Inf, ymax=+Inf),
              fill=back.colors[1], alpha=0.5, inherit.aes=FALSE) +
    geom_rect(alpha=0.9, aes(fill=type)) +
    scale_y_continuous(breaks=pretty(unique(reads.gene1$pair.id), n=4)) +
    scale_x_continuous(expand=c(0, 0)) +
    theme(axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.text.y=element_blank(),
          panel.background=element_rect(fill="transparent"),
          panel.border=element_rect(color="transparent", fill=NA),
          strip.text=element_blank(),
          strip.background=element_blank(),
          legend.text=element_text(size=7),
          legend.key=element_rect(size=0.5)) +
          ## top, right, bottom, left
          ##plot.margin=unit(c(0.1, 0, 0.01, 0), "npc")) +
    ##ylab("rearranged\nread pairs") + xlab("") +
    guides(fill="none")
  if(strand1 == "-") g1 <- g1 + scale_x_reverse(expand=c(0, 0))
  g2 <- ggplot(reads.gene2,
         aes(xmin=start, xmax=end, ymin=pair.id-0.3,
             ymax=pair.id+0.3)) +
    geom_rect(aes(xmin=-Inf, xmax=+Inf, ymin=-Inf, ymax=+Inf),
              fill=back.colors[2], alpha=0.5, inherit.aes=FALSE) +
    geom_rect(alpha=0.9, aes(fill=type)) +
    scale_x_continuous(expand=c(0, 0)) +
    theme(axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.text.y=element_blank(),
          panel.background=element_rect(fill="transparent"),
          panel.border=element_rect(color="transparent", fill=NA),

          strip.text=element_blank(),
          strip.background=element_blank(),
          legend.text=element_text(size=7),
          legend.key=element_rect(size=0.5),
          ## top, right, bottom, left
          ## plot.margin=unit(c(0.1, 0, 0.01, 0), "npc"),  +
          ## plot.margin=unit(c(0.1, 0.45, 0.01, 0), "npc"),
          legend.key.size=unit(0.5, "cm")) +
    ylab("") + xlab("") +
    guides(fill=guide_legend(title=""))
  if(strand2=="-") g2 <- g2 + scale_x_reverse(expand=c(0, 0))
  legend <- getLegend(g2)
  g2 <- g2 + theme(legend.position="none")
  gg1 <- ggplotGrob(g1)
  gg2 <- ggplotGrob(g2)
  gg2$widths <- gg1$widths
  gg2$heights <- gg1$heights
  list(gg1, gg2, legend)
}

ggProtein <- function(domains, params){
  protein <- params$protein
  xlimit <- c(1, domains$aa_len[1])
  domains$description <- factor(domains$description)
  domain.color <- params$domain.color
  back.color <- params$background.color
  description.size <- params$description.size
  ggplot(domains, aes(xmin=start, xmax=end,
                      ymin=-0.5, ymax=0.5,
                      label=short.desc)) +
    geom_rect(data=domains, aes(xmin=xlimit[1], xmax=xlimit[2],
                                ymin=-0.5, ymax=0.5),
              color="black", fill=back.color, inherit.aes=FALSE) +
    geom_rect(fill=domain.color) +
    geom_text(aes(x=midx, y=0), size=description.size, angle=90) +
    ##scale_fill_manual(values=c(params[["background.color"]],
    ##rep(params[["domain.color"]], 3))) +
    scale_x_continuous(expand=c(0, 1), limits=xlimit) +
    ##ylab(protein) +
    theme(axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.y=element_blank(),
          ##axis.ticks=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_rect(fill="transparent"))
}


ggClippedProtein1 <- function(p.dat, params){
  manual.colors <- c(params[["clipped.color"]], params[["domain.color"]])
  protein <- params[["protein"]]
  xlimit <- c(1, p.dat$aa_len[1])
  ##i <- as.integer(factor(p.dat$is.clipped, levels=c(TRUE, FALSE)))
  ##manual.colors <- manual.colors[i]
  back.color <- params$background.color
  description.size <- params$description.size
  jxn <- p.dat$aa.jxn[1]
  jxn.in.dom <- jxn > p.dat$start & jxn < p.dat$end
  if(any(jxn.in.dom)){
    index <- which(jxn.in.dom)
    end <- p.dat$end[index]
    p.dat$end[index] <- jxn
    p.dat$is.clipped[index] <- FALSE
    clipped.dom <- p.dat[1, ]
    clipped.dom$start <- jxn
    clipped.dom$end <- end
    clipped.dom$is.clipped <- TRUE
    p.dat <- rbind(p.dat, clipped.dom)
  }
  p.dat$is.clipped <- factor(p.dat$is.clipped, levels=c("TRUE", "FALSE"))
  ## if all TRUE or all FALSE, ggplot will only use the first color
  manual.colors <- unique(manual.colors[sort(as.integer(p.dat$is.clipped))])
  ggplot(data=p.dat,
         aes(xmin=start, xmax=end, ymin=-0.5, ymax=0.5)) +
    geom_rect(data=p.dat, aes(xmin=xlimit[1], xmax=xlimit[2],
                                ymin=-0.5, ymax=0.5),
              color="black", fill=back.color, inherit.aes=FALSE) +
    geom_rect(data=p.dat,
              aes(xmin=p.dat$clip.start[1], ##start(rreads)[1]+20,
                  xmax=p.dat$clip.end[1],
                  ymin=-0.5, ymax=+0.5),
              fill="gray",
              color="transparent", alpha=0.9, inherit.aes=FALSE) +
    geom_rect(aes(fill=is.clipped, color=is.clipped), size=1) +
    geom_text(aes(x=midx, y=0, label=short.desc),
              size=description.size, angle=90)  +
    scale_x_continuous(expand=c(0, 1),
                       breaks=p.dat$clip.start[1],
                       labels=p.dat$clip.start[1]) +
    guides(fill="none", color="none") +
    geom_vline(xintercept=p.dat$aa.jxn[1], linetype="dashed") +
    scale_color_manual(values=manual.colors) +
    scale_fill_manual(values=manual.colors) +
    theme(axis.title=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          panel.background=element_rect(fill="transparent"),
          panel.grid.minor.x=element_line(color="transparent"),
          panel.border=element_rect(color="transparent", fill=NA))
  ## top, right, bottom, left
  ##plot.margin=unit(c(1, -0.1, 1, 1), "lines"))
  ##ggtitle(protein)
}

ggClippedProtein2 <- function(p.dat, params){
  manual.colors <- c(params[["clipped.color"]], params[["domain.color"]])
  protein <- params[["protein"]]
  xlimit <- c(1, p.dat$aa_len[1])
  ##i <- as.integer(factor(p.dat$is.clipped, levels=c(TRUE, FALSE)))
  ##manual.colors <- manual.colors[i]
  back.color <- params$background.color
  description.size <- params$description.size
  jxn <- p.dat$aa.jxn[1]
  jxn.in.dom <- jxn > p.dat$start & jxn < p.dat$end
  if(any(jxn.in.dom)){
    index <- which(jxn.in.dom)
    end <- p.dat$end[index]
    p.dat$end[index] <- jxn
    p.dat$is.clipped[index] <- FALSE
    clipped.dom <- p.dat[1, ]
    clipped.dom$start <- jxn
    clipped.dom$end <- end
    clipped.dom$is.clipped <- TRUE
    p.dat <- rbind(p.dat, clipped.dom)
  }
  p.dat$is.clipped <- factor(p.dat$is.clipped, levels=c("TRUE", "FALSE"))
  ## if all TRUE or all FALSE, ggplot will only use the first color
  manual.colors <- unique(manual.colors[sort(as.integer(p.dat$is.clipped))])
  ggplot(data=p.dat,
         aes(xmin=start, xmax=end, ymin=-0.5, ymax=0.5)) +
    geom_rect(data=p.dat, aes(xmin=xlimit[1], xmax=xlimit[2],
                                ymin=-0.5, ymax=0.5),
              color="black", fill=back.color, inherit.aes=FALSE) +
    geom_rect(data=p.dat,
              aes(xmin=p.dat$clip.start[1], ##start(rreads)[1]+20,
                  xmax=p.dat$clip.end[1],
                  ymin=-0.5, ymax=+0.5),
              fill="gray",
              color="transparent", alpha=0.9, inherit.aes=FALSE) +
    geom_rect(aes(fill=is.clipped, color=is.clipped), size=1) +
    geom_text(data=p.dat, aes(x=midx, y=0, label=short.desc),
              size=description.size, angle=90, inherit.aes=FALSE)  +
    scale_x_continuous(expand=c(0, 1), breaks=p.dat$clip.end[1], labels=p.dat$clip.end[1]) +
    guides(fill="none", color="none") +
    geom_vline(xintercept=p.dat$aa.jxn[1], linetype="dashed") +
    scale_color_manual(values=manual.colors) +
    scale_fill_manual(values=manual.colors) +
    theme(axis.title=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          panel.background=element_rect(fill="transparent"),
          panel.grid.minor.x=element_line(color="transparent"),
          panel.border=element_rect(color="transparent", fill=NA))
  ## top, right, bottom, left
  ##plot.margin=unit(c(1, -0.1, 1, 1), "lines"))
  ##ggtitle(protein)
}

ggProteinFusion <- function(data.list, p1.params, p2.params){
  coords <- data.list[["coords"]]
  d <- data.list[["domains"]]
  coords$ymin <- -0.5
  coords$ymax <- 0.5
  proteins <- c(p1.params[["protein"]], p2.params[["protein"]])
  domain.colors <- c(p1.params[["domain.color"]],
                      p2.params[["domain.color"]])
  names(domain.colors) <- proteins
  back.color1 <- p1.params[["background.color"]]
  back.color2 <- p2.params[["background.color"]]
  if(nrow(d)==0){
    dtext <- data.frame(short.desc="", midx=1)
    d$ymin <- rep(-0.5, nrow(d))
    d$ymax <- rep(0.5, nrow(d))
  } else{
    dtext <- d
    d$ymin <- -0.5
    d$ymax <- 0.5
  }
  jxn <- coords$end[1]
  breaks <- c(1, jxn, coords$end[2])
  description.size <- p1.params$description.size
  ggplot(coords,
         aes(xmin=start, xmax=end, ymin=ymin,
             ymax=ymax, color=hugo,
             fill=hugo)) +
    geom_rect(fill="transparent", color="black") +
    ## background for first gene
    geom_rect(data=coords[1, ],
              aes(xmin=start, xmax=end, ymin=ymin,
                  ymax=ymax),
              color="transparent",
              fill=back.color1, inherit.aes=FALSE) +
    ## background for second gene
    geom_rect(data=coords[2, ],
              aes(xmin=start, xmax=end, ymin=ymin,
                  ymax=ymax),
              color="transparent",
              fill=back.color2, inherit.aes=FALSE) +
    ## domains
    geom_rect(data=d, aes(xmin=start, xmax=end,
                           ymin=ymin, ymax=ymax,
                           color=hugo,
                           fill=hugo),
              inherit.aes=FALSE, alpha=0.9) +
    geom_text(data=dtext, aes(x=midx, y=0, label=short.desc),
              size=description.size, angle=90, inherit.aes=FALSE) +
    ## background for first gene
    geom_vline(xintercept=d$aa.jxn[1],
               color="gray20", linetype="dashed") +
    ylab("") + xlab("") +
    guides(fill="none", color="none") +
    scale_color_manual(values=domain.colors) +
    scale_fill_manual(values=domain.colors) +
    scale_x_continuous(expand=c(0, 1), breaks=breaks, labels=breaks) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background=element_rect(fill="transparent"),
          panel.grid.minor.x=element_line(color="transparent"),
          panel.border=element_rect(color="transparent", fill=NA))
  ## top, right, bottom, left
  ##plot.margin=unit(c(0.5, 1, 0, 1), "mm"))
}

compositeFusionParams <- function(){
  widths <- c(0.95, 0.05, 0.95, 0.05)
  widths <- widths/sum(widths)
  widths <- unit(widths, "npc")

  heights1 <- c(0.3, rep(0.7, 3))
  heights2 <- heights1[2:4]
  heights3 <- c(heights1, heights2)
  heights4 <- heights3/sum(heights3)
  heights <- unit(heights4, "npc")

  layout1 <- rbind(rep(1, 4), 2:5,
                   c(6, 6, 7, 7),
                   c(8, 8, 8, 9))
  layout2 <- rbind(c(1, 1, 2, 2),
                   c(3, 3, 4, 4),
                   c(5, 6, 6, 5)) + max(layout1)
  layout <- rbind(layout1, layout2)
  list(widths=widths, heights=heights, layout=layout)
}

gg_blank <- function(){
  df <- data.frame(start=1:10, end=1:10)
  ggplot(df, aes(start, end)) + geom_blank() +
    theme(text=element_blank(),
          axis.ticks=element_blank(),
          panel.background=element_rect(fill="transparent"))
}

#' @export
transcriptGrobs <- function(dat){
  g1.clip <- ggClippedExons(dat[["fusion.dat"]],
                            dat[["g.params"]][[1]],
                            dat[["roi"]])
  g.axis <- ggAxisLabel(dat[["fusion.dat"]],
                        dat[["g.params"]][[1]], "3'")
  gg1.clip <- ggplotGrob(g1.clip)
  gg.axis <- ggplotGrob(g.axis)
  gg.axis$heights <- gg1.clip$heights
  g2.clip <- ggClippedExons(dat[["fusion.dat"]],
                            dat[["g.params"]][[2]],
                            dat[["roi"]])
  gg2.clip <- ggplotGrob(g2.clip)

  g.fused <- ggTxFusion(dat[["fused.transcripts"]],
                        dat[["g.params"]][[1]],
                        dat[["g.params"]][[2]])
  g.fused <- g.fused + ggtitle(dat[["fusion.dat"]]$fusion)
##  fused1 <- ggClippedExons(dat[["fusion.dat"]],
##                           dat[["g.params"]][[1]],
##                           dat[["roi"]],
##                           clip=TRUE)
##  fused2 <- ggClippedExons(dat[["fusion.dat"]],
##                           dat[["g.params"]][[2]],
##                           dat[["roi"]],
##                           clip=TRUE)
##
  gg.fused <- ggplotGrob(g.fused)
  gg.fused$heights <- gg1.clip$heights
  gg.fused$widths <- gg1.clip$widths

##  gg.fused <- list(ggplotGrob(fused1), ggplotGrob(fused2))
  gg.readlist <- ggRearrangedReads2(dat[["rreads"]], dat[["g.params"]],
                                    dat[["roi"]], dat[["reverse.roi"]])

  list(gg1.clip=gg1.clip,
       gg2.clip=gg2.clip,
       gg.axis=gg.axis,
       gg.readlist=gg.readlist,
       gg.fused=gg.fused)
}


#' @export
proteinGrobs <- function(dat){
  ## plot protein
  g.p1 <- ggProtein(dat[["protein1"]],
                    dat[["protein1.params"]])
  g.p2 <- ggProtein(dat[["protein2"]],
                    dat[["protein2.params"]])
  gp1.clip <- ggClippedProtein1(dat[["protein1.clipped"]],
                                dat[["protein1.params"]])
  gp2.clip <- ggClippedProtein2(dat[["protein2.clipped"]],
                                dat[["protein2.params"]])
  ##load_all("integration/integration")
  gp.fuse <- ggProteinFusion(dat[["protein.fusion"]],
                             dat[["protein1.params"]],
                             dat[["protein2.params"]])
  list(g.p1=g.p1, g.p2=g.p2,
       gp.fuse=gp.fuse,
       gp1.clip=gp1.clip,
       gp2.clip=gp2.clip)
}

arrangeFusionGrobs <- function(g.ideo, tx, pr, layout){
  if(missing(layout)){
    layout <- compositeFusionParams()
  }
  blank <- gg_blank()
  grid.arrange(g.ideo,
               tx[["gg1.clip"]],
               tx[["gg.axis"]],
               tx[["gg2.clip"]],
               tx[["gg.axis"]],
               tx[["gg.readlist"]][[1]],
               tx[["gg.readlist"]][[2]],
               tx[["gg.fused"]],
               tx[["gg.axis"]],
               pr[["g.p1"]],
               pr[["g.p2"]],
               pr[["gp1.clip"]],
               pr[["gp2.clip"]],
               blank,
               pr[["gp.fuse"]],
               widths=layout[["widths"]],
               heights=layout[["heights"]],
               layout_matrix=layout[["layout"]],
               newpage=FALSE)
}

arrangeFusionGrobs2 <- function(tx.grobs, pt.grobs){
  pt.widths <- c(tx.protein$protein1$aa_len[1],
                 tx.protein$protein2$aa_len[1])
  pt.widths1 <- pt.widths/(sum(pt.widths))

  p1 <- tx.protein$protein1.clipped
  p2 <- tx.protein$protein2.clipped
  pt.widths2 <- c(p1$clip.start[1],
                  p2$aa_len[1]-p2$clip.end[1])
  pt.width2 <- sum(pt.widths2)/sum(pt.widths)

  roi <- tx.protein[["roi"]]
  ##
  ## Draw the transcripts and proteins to scale
  ##
  chroms <- as.character(seqnames(roi))
  g1.clipped <- GRanges(chroms[1], IRanges(start(roi)[1], roi$bp.jxn[1]))
  g2.clipped <- GRanges(chroms[2], IRanges(start(roi)[2], roi$bp.jxn[2]))
  widths <- width(roi)/sum(width(roi))
  width2 <- sum(width(c(g1.clipped, g2.clipped)))/sum(width(roi))
  ##fused.width <- sum(width(roi2))/sum(width(roi))
  just <- c("left", "bottom")
  vp <- viewport(x=unit(0, "npc"), y=unit(0.3, "npc"),
                 width=unit(1, "npc"), height=unit(0.7, "npc"), just=c("left", "bottom"))
  pushViewport(vp)
  txLayout1(g.ideo, 1)
  txLayout2(tx.grobs[["gg.readlist"]][[1]], nrow=2, x=0.2, w=0.2)
  txLayout2(tx.grobs[["gg.readlist"]][[2]], nrow=2, x=0.4, w=0.2)
  txLayout2(tx.grobs[["gg.readlist"]][[3]], nrow=2, x=0.6, w=0.05)
  txLayout2(tx.grobs[["gg1.clip"]], nrow=3, x=0, w=widths[1])
  txLayout2(tx.grobs[["gg2.clip"]], nrow=3, x=widths[1], w=widths[2])
  txLayout2(tx.grobs[["gg.fused"]], nrow=4, x=0, w=width2-.014)
  upViewport()
  vp <- viewport(x=unit(0, "npc"), y=unit(0, "npc"),
                 width=unit(1, "npc"), height=unit(0.25, "npc"), just=c("left", "bottom"))
  pushViewport(vp)
  ptLayout2(mygrob=ggplotGrob(pr.grobs[["gp1.clip"]]), nrow=1, x=0, w=pt.widths1[1])
  ptLayout2(mygrob=ggplotGrob(pr.grobs[["gp2.clip"]]), nrow=1, x=pt.widths1[1], w=pt.widths1[2])
  ptLayout2(mygrob=ggplotGrob(pr.grobs[["gp.fuse"]]), nrow=2, x=0-0.02, w=pt.width2-0.01)
}

txLayout1 <- function(mygrob, nrow, just=c("left", "center")){
  vp <- viewport(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                 width=unit(0.95, "npc"), height=unit(0.95, "npc"),
                 just=rep("center", 2))
  pushViewport(vp)
  gl <- grid.layout(nrow=4, ncol=1,
                    heights=unit(c(2, 4, 2, 2), "null"))
  vp0 <- viewport(layout=gl, x=unit(0.5, "npc"),
                  y=unit(0.5, "npc"))
  pushViewport(vp0)
  vp2 <- viewport(layout.pos.row=nrow, just=just)
  pushViewport(vp2)
  ##grid.rect()
  grid.draw(mygrob)
  popViewport(3)
}

ptLayout1 <- function(mygrob, nrow, just=c("left", "center")){
  vp <- viewport(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                 width=unit(0.95, "npc"), height=unit(0.95, "npc"),
                 just=rep("center", 2))
  pushViewport(vp)
  gl <- grid.layout(nrow=3, ncol=1,
                    heights=unit(c(1, 1, 1), "null"))
  vp0 <- viewport(layout=gl, x=unit(0.5, "npc"),
                  y=unit(0.5, "npc"))
  pushViewport(vp0)
  vp2 <- viewport(layout.pos.row=nrow, just=just)
  pushViewport(vp2)
  grid.draw(mygrob)
  popViewport(3)
}

txLayout2 <- function(mygrob, nrow, x, w, just=c("left", "center")){
  vp <- viewport(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                 width=unit(0.95, "npc"), height=unit(0.95, "npc"),
                 just=rep("center", 2))
  pushViewport(vp)
  gl <- grid.layout(nrow=4, ncol=1,
                    heights=unit(c(2, 4, 2, 2), "null"))
  vp0 <- viewport(layout=gl, x=unit(0.5, "npc"),
                  y=unit(0.5, "npc"))
  pushViewport(vp0)
  vp1 <- viewport(layout.pos.row=nrow)
  pushViewport(vp1)
  ##grid.rect()
  vp2 <- viewport(x=unit(x, "npc"),
                  y=unit(0.5, "npc"),
                  width=unit(w, "npc"),
                  just=just)
  pushViewport(vp2)
  ##grid.rect(gp=gpar(col="blue"))
  grid.draw(mygrob)
  popViewport(4)
}

ptLayout2 <- function(mygrob, nrow, x, w, just=c("left", "center")){
  vp <- viewport(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                 width=unit(0.95, "npc"), height=unit(0.95, "npc"),
                 just=rep("center", 2))
  pushViewport(vp)
  gl <- grid.layout(nrow=2, ncol=1,
                    heights=unit(c(1, 1), "null"))
  vp0 <- viewport(layout=gl, x=unit(0.5, "npc"),
                  y=unit(0.5, "npc"))
  pushViewport(vp0)
  vp1 <- viewport(layout.pos.row=nrow)
  pushViewport(vp1)
  ##grid.rect()
  vp2 <- viewport(x=unit(x, "npc"),
                  y=unit(0.5, "npc"),
                  width=unit(w, "npc"),
                  just=just)
  pushViewport(vp2)
  ##grid.rect(gp=gpar(col="blue"))
  grid.draw(mygrob)
  popViewport(4)
}


.pt_viewports <- function(roi){
  gl2 <- grid.layout(nrow=4, ncol=1,
                     widths=unit(1, "npc"),
                     heights=unit(c(2, 4, 2, 2), "null"))
  gl.pt <- grid.layout(nrow=1, ncol=2,
                       widths=unit(pt.widths, "npc"))

}
