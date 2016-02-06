##
## https://github.com/mylesmharrison/colorRampPaletteAlpha/blob/master/colorRampPaletteAlpha.R
##
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

colorRampPaletteAlpha <- function(colors, n=32, interpolate='linear') {
  # Create the color ramp normally
  cr <- colorRampPalette(colors, interpolate=interpolate)(n)
  # Find the alpha channel
  a <- col2rgb(colors, alpha=T)[4,]
  # Interpolate
        if (interpolate=='linear') {
                l <- approx(a, n=n)
        } else {
                l <- spline(a, n=n)
        }
        l$y[l$y > 255] <- 255 # Clamp if spline is > 255
        cr <- addalpha(cr, l$y/255.0)
        return(cr)
}

.cytobandColors <- function(stains){
  colors <- rep("white", length(stains))
  colors[stains == "gneg" | stains == "gvar"] <- "gray100"
  colors[stains=="gpos25"] <- "gray90"
  colors[stains=="gpos50"] <- "gray70"
  colors[stains=="gpos75"] <- "gray40"
  colors[stains=="gpos100"] <- "gray0"
  colors[stains=="stalk"] <- "brown3"
  colors[stains=="acen"] <- "brown4"
  colors
}

#' Create qualitative color palettes of different lengths
#'
#' TODO: Reference URL where these colors were obtained
#'
#' @seealso \code{\link{plot_amplicons}}
#' @return a list
#' @examples
#' pal <- qualitativeColors()
#' pal[[3]]
#' @export
qualitativeColors <- function(){
  tol1qualitative=c("#4477AA")
  tol2qualitative=c("#4477AA", "#CC6677")
  tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
  tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
  tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
  tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677", "#AA4499")
  tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77",
                    "#CC6677","#AA4499")
  tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
                    "#DDCC77", "#CC6677","#AA4499")
  tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
                    "#DDCC77", "#CC6677", "#882255", "#AA4499")
  tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
                     "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
  tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733",
                     "#999933", "#DDCC77", "#661100", "#CC6677", "#882255",
                     "#AA4499")
  tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733",
                     "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466",
                     "#882255", "#AA4499")
  tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA",
                  "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744",
                  "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77",
                  "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
  col_list <- list(tol1qualitative,
                   tol2qualitative,
                   tol3qualitative,
                   tol4qualitative,
                   tol5qualitative,
                   tol6qualitative,
                   tol7qualitative,
                   tol8qualitative,
                   tol9qualitative,
                   tol10qualitative,
                   tol11qualitative,
                   tol12qualitative,
                   tol21rainbow)
}


## From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
## To use for fills, add
##  scale_fill_manual(values=cbPalette)
##
## To use for line and point colors, add
##  scale_colour_manual(values=cbPalette)
## The palette with grey:
cbPalette <- function(){
  c("#999999", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}

## The palette with black:
cbbPalette <- function(){
  c("#000000", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}

