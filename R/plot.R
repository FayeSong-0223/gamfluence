#' Plot fixed against added deletion changes
#'
#' One point per site. The x-axis is `fixed`, the change in the other sites'
#' fitted values when the site is removed with the hyperparameters held; the
#' y-axis is `added`, the extra change when they are re-estimated. Both are
#' weighted RMS changes on the link scale (see [gam_influence()]), so the axes
#' share one scale. Sites above the dashed diagonal change the fit more through
#' the smoothing parameters (and theta) than through the held-parameter refit.
#'
#' The `n` sites with the largest `by` are labelled. Sites whose re-estimated
#' refit had two REML optima (status `"two optima"`) are drawn as orange
#' triangles; read their `added` with care. Sites whose refits failed have no
#' changes to plot and are left out, with a message.
#'
#' @param x A `gam_influence` object.
#' @param n Number of sites to label.
#' @param by Which change ranks the sites to label: `"total"`, `"added"` or
#'   `"fixed"`.
#' @param main Plot title.
#' @param ... Unused.
#' @return Invisibly, a data frame of the plotted sites (`id`, `n_obs`,
#'   `fixed`, `added`, `total`, `status`) with a logical column `labelled`.
#' @examples
#' library(mgcv)
#' set.seed(1)
#' d <- data.frame(site = factor(rep(1:12, each = 6)), x = runif(72))
#' d$y <- sin(2 * pi * d$x) + rnorm(12, 0, 0.4)[d$site] + rnorm(72, 0, 0.5)
#' m <- gam(y ~ s(x) + s(site, bs = "re"), data = d, method = "REML")
#' infl <- gam_influence(m, cluster = "site", progress = FALSE)
#' plot(infl, n = 3)
#' @export
plot.gam_influence <- function(x, n = 5, by = "total",
                               main = "Site deletion: fixed vs added change", ...) {
  if (!(by %in% c("fixed", "added", "total")))
    stop("`by` must be \"fixed\", \"added\" or \"total\".", call. = FALSE)
  s  <- x$sites
  ok <- is.finite(s$fixed) & is.finite(s$added)
  if (any(!ok))
    message(sum(!ok), " site(s) not shown: their refits failed or did not converge.")
  s <- s[ok, c("id", "n_obs", "fixed", "added", "total", "status"), drop = FALSE]
  rownames(s) <- NULL
  if (!nrow(s)) stop("no site has changes to plot.", call. = FALSE)

  two <- grepl("two optima", s$status, fixed = TRUE)
  lab <- rep(FALSE, nrow(s))
  lab[utils::head(order(s[[by]], decreasing = TRUE), max(0, n))] <- TRUE

  col <- c(point = "#2a78d6", flag = "#eb6834", ink = "#0b0b0b",
           ink2 = "#52514e", grid = "#e6e5e1", ref = "#8c8b86")
  pt_cex <- 1.3
  lab_cex <- 0.8
  top <- max(c(s$fixed, s$added))
  lim <- c(0, if (top > 0) top * 1.1 else 1)

  op <- graphics::par(mar = c(4.4, 4.6, 4.2, 1.2), mgp = c(2.7, 0.7, 0), las = 1,
                      tcl = -0.3, col.axis = col[["ink2"]], col.lab = col[["ink"]],
                      cex.axis = 0.8)
  on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(lim, lim, asp = 1, xaxs = "i", yaxs = "i")
  usr <- graphics::par("usr")
  graphics::abline(v = graphics::axTicks(1), h = graphics::axTicks(2),
                   col = col[["grid"]], lwd = 0.8)
  graphics::abline(0, 1, col = col[["ref"]], lty = 2, lwd = 1.2)
  graphics::axis(1, col = NA, col.ticks = col[["ref"]])
  graphics::axis(2, col = NA, col.ticks = col[["ref"]])
  graphics::box(bty = "l", col = col[["ref"]])
  graphics::title(xlab = "fixed: change with hyperparameters held",
                  ylab = "added: extra change when re-estimated")
  graphics::title(main = main, adj = 0, line = 2.6, font.main = 2, cex.main = 1,
                  col.main = col[["ink"]])
  graphics::mtext(paste("Weighted RMS change in the other sites' fitted values (link scale).",
                        "Above the line, re-estimation adds more.", sep = "\n"),
                  side = 3, line = 0.4, adj = 0, cex = 0.72, col = col[["ink2"]])

  # the diagonal's label, along the line near its upper end
  dx <- usr[2] - usr[1]
  graphics::text(usr[1] + 0.88 * dx, usr[1] + 0.88 * dx, "added = fixed", srt = 45,
                 adj = c(0.5, -0.5), cex = 0.7, col = col[["ref"]])

  graphics::points(s$fixed[!two], s$added[!two], pch = 21, bg = col[["point"]],
                   col = "white", cex = pt_cex, lwd = 1)
  graphics::points(s$fixed[two], s$added[two], pch = 24, bg = col[["flag"]],
                   col = "white", cex = pt_cex, lwd = 1)

  obstacles <- list()
  if (any(two)) {
    corner <- emptiest_corner(s$fixed, s$added, usr)
    lg <- graphics::legend(corner, legend = c("site", "two REML optima"),
                           pch = c(21, 24), pt.bg = c(col[["point"]], col[["flag"]]),
                           col = "white", pt.cex = pt_cex, bty = "n", cex = 0.75,
                           text.col = col[["ink2"]], inset = 0.01)
    r <- lg$rect
    obstacles <- list(c(r$left, r$left + r$w, r$top - r$h, r$top))
  }
  if (any(lab)) {
    i <- which(lab)[order(s[[by]][lab], decreasing = TRUE)]
    place_labels(s$fixed[i], s$added[i], s$id[i], s$fixed, s$added, obstacles,
                 usr, lab_cex, pt_cex, col[["ink"]], col[["ink2"]])
  }
  s$labelled <- lab
  invisible(s)
}

# Corner of the plot region with the fewest points (for the legend).
emptiest_corner <- function(px, py, usr) {
  w <- usr[2] - usr[1]; h <- usr[4] - usr[3]
  box <- list(bottomright = c(usr[2] - 0.45 * w, usr[2], usr[3], usr[3] + 0.18 * h),
              topleft     = c(usr[1], usr[1] + 0.45 * w, usr[4] - 0.18 * h, usr[4]),
              bottomleft  = c(usr[1], usr[1] + 0.45 * w, usr[3], usr[3] + 0.18 * h),
              topright    = c(usr[2] - 0.45 * w, usr[2], usr[4] - 0.18 * h, usr[4]))
  n_in <- vapply(box, function(b) sum(px >= b[1] & px <= b[2] & py >= b[3] & py <= b[4]),
                 numeric(1))
  names(box)[which.min(n_in)]
}

# Greedy label placement: for each label in turn, try positions around its
# point at increasing distances and keep the first that stays inside the plot
# and clear of points, the legend and labels already placed. Labels placed
# away from their point, or above or below it, get a thin leader line.
place_labels <- function(px, py, labels, all_x, all_y, obstacles, usr,
                         cex, pt_cex, ink, ink2) {
  ch <- graphics::strheight("M", cex = cex)
  pr <- 0.5 * graphics::strheight("M", cex = pt_cex)
  pad <- 0.15 * ch
  hits <- function(a, b) !(a[2] + pad < b[1] || b[2] + pad < a[1] ||
                           a[4] + pad < b[3] || b[4] + pad < a[3])
  dirs <- list(c(1, 0), c(-1, 0), c(1, 1), c(-1, 1), c(1, -1), c(-1, -1),
               c(0, 1), c(0, -1))
  radii <- c(1.1, 2.2, 3.4, 4.8) * ch
  boxes <- obstacles
  for (k in seq_along(labels)) {
    w <- graphics::strwidth(labels[k], cex = cex)
    chosen <- NULL
    for (r in radii) {
      for (d in dirs) {
        u  <- d / sqrt(sum(d^2))
        ax <- px[k] + u[1] * r
        ay <- py[k] + u[2] * r
        x0 <- if (u[1] > 0.1) ax else if (u[1] < -0.1) ax - w else ax - w / 2
        y0 <- if (u[2] > 0.1) ay else if (u[2] < -0.1) ay - ch else ay - ch / 2
        b  <- c(x0, x0 + w, y0, y0 + ch)
        if (b[1] < usr[1] || b[2] > usr[2] || b[3] < usr[3] || b[4] > usr[4]) next
        if (any(vapply(boxes, hits, logical(1), b = b))) next
        if (any(all_x > b[1] - pr & all_x < b[2] + pr &
                all_y > b[3] - pr & all_y < b[4] + pr)) next
        chosen <- list(box = b, r = r, side = u[2] == 0)
        break
      }
      if (!is.null(chosen)) break
    }
    if (is.null(chosen))   # nowhere free: put it to the right and accept overlap
      chosen <- list(box = c(px[k] + radii[1], px[k] + radii[1] + w,
                             py[k] - ch / 2, py[k] + ch / 2), r = radii[1], side = TRUE)
    b <- chosen$box
    boxes <- c(boxes, list(b))
    if (chosen$r > radii[1] || !chosen$side) {   # leader unless right beside the point
      cx <- min(max(px[k], b[1]), b[2])
      cy <- min(max(py[k], b[3]), b[4])
      len <- sqrt((cx - px[k])^2 + (cy - py[k])^2)
      gap <- 0.3 * ch                                 # keep the line off the text
      if (len - pr - gap > 0.2 * ch)
        graphics::segments(px[k] + (cx - px[k]) * pr / len, py[k] + (cy - py[k]) * pr / len,
                           cx - (cx - px[k]) * gap / len, cy - (cy - py[k]) * gap / len,
                           col = ink2, lwd = 0.7)
    }
    graphics::text(mean(b[1:2]), mean(b[3:4]), labels[k], cex = cex, col = ink)
  }
  invisible(utils::tail(boxes, length(labels)))       # the label boxes, for tests
}
