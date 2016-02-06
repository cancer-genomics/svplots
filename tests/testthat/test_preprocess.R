context("Plot pviews")
test_that("genomeLogRatios", {
  library(svdata)
  data(pviews)
  paths(pviews) <- file.path(system.file("extdata", package="svdata"),
                             rdsId(pviews))
  logratios <- genomeLogRatios(pviews, thin=20)
  expect_is(logratios, "GRanges")
  expect_is(logratios$log_ratio, "numeric")
})

test_that("genomeLogRatios", {
  library(svdata)
  data(pviews)
  paths(pviews) <- file.path(system.file("extdata", package="svdata"),
                             rdsId(pviews))
  logratios <- genomeLogRatios(pviews, thin=50)
  fig <- plotGenomeLogRatios(logratios)
  expect_is(fig, "GGbio")

  colors <- rep("black", length(logratios))
  index <- chromosome(logratios) %in% paste0("chr", seq(1, 22, 2))
  colors[index] <- "gray"
  plot(logratios$log_ratio, pch=".", col=colors)
})
