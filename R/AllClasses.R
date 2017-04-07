#' @include help.R
NULL

setClass("IdiogramParams", representation(seqnames="character",
                                          seqlengths="numeric",
                                          unit="character",
                                          genome="character",
                                          box="list"))

setClass("ILimit", contains="IRanges")
