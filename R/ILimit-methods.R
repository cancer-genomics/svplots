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
