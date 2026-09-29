# plot(): one point per site, fixed against added, top sites labelled.

plot_fixture <- function() {
  gam_influence(reef_models$poisson_off, cluster = "site", which = as.character(1:8),
                progress = FALSE)
}

test_that("plot draws every site and labels the top n", {
  i <- plot_fixture()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  out <- plot(i, n = 3)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("id", "n_obs", "fixed", "added", "total", "status", "labelled"))
  expect_equal(nrow(out), 8)
  expect_equal(sum(out$labelled), 3)
  expect_setequal(out$id[out$labelled], top_units(i, by = "total", n = 3)$id)
  expect_setequal(plot(i, n = 2, by = "added")$id[plot(i, n = 2, by = "added")$labelled],
                  top_units(i, by = "added", n = 2)$id)
  expect_equal(sum(plot(i, n = 0)$labelled), 0)
  expect_equal(sum(plot(i, n = 50)$labelled), 8)
  expect_error(plot(i, by = "nope"), "must be")
})

test_that("failed sites are left out with a message; flagged sites are drawn", {
  i <- plot_fixture()
  i$sites[1, c("fixed", "added", "total")] <- NA
  i$sites$status[2] <- "two optima"
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_message(out <- plot(i), "1 site\\(s\\) not shown")
  expect_equal(nrow(out), 7)
  expect_true(i$sites$id[2] %in% out$id)
  i$sites[, c("fixed", "added", "total")] <- NA
  expect_error(suppressMessages(plot(i)), "no site")
})

test_that("placed labels overlap neither each other nor the legend", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  graphics::plot.new()
  graphics::plot.window(c(0, 1), c(0, 1), asp = 1, xaxs = "i", yaxs = "i")
  set.seed(3)
  px <- 0.5 + stats::runif(10, -0.05, 0.05)       # a tight cluster of ten sites
  py <- 0.5 + stats::runif(10, -0.05, 0.05)
  leg <- c(0.75, 1, 0, 0.12)                      # a legend box in the corner
  b <- gamfluence:::place_labels(px, py, paste0("site", 1:10), px, py, list(leg),
                                 graphics::par("usr"), 0.8, 1.3, "black", "grey40")
  expect_length(b, 10)
  hit <- function(a, c) !(a[2] < c[1] || c[2] < a[1] || a[4] < c[3] || c[4] < a[3])
  pairs <- utils::combn(10, 2)
  expect_false(any(apply(pairs, 2, function(ij) hit(b[[ij[1]]], b[[ij[2]]]))))
  expect_false(any(vapply(b, hit, logical(1), c = leg)))
  inside <- vapply(b, function(z) z[1] >= 0 && z[2] <= 1 && z[3] >= 0 && z[4] <= 1, logical(1))
  expect_true(all(inside))
})
