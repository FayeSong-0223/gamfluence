## Resubmission

This is a resubmission of a new package. In response to the review, I have:

* Added `\value` to `gam_influence-methods.Rd`. It now describes the class
  and the meaning of what each method returns: `print()`, `summary()`, the
  `print()` method for summaries, `as.data.frame()` and `top_units()`. The
  print methods are documented as returning their input invisibly, called for
  their side effect of printing. The `\value` of `plot.gam_influence.Rd` now
  also says that it is called for its side effect of drawing the plot.
* Removed the `<<-` in `gam_influence()`. Warnings from each site's refits are
  now collected in a local environment created inside the function, so the
  package does not write to `.GlobalEnv`. A new test checks that
  `gam_influence()` leaves the global environment unchanged.

## Test environments

* GitHub Actions: macOS (R release), Windows (R release), Ubuntu (R devel,
  R release, R oldrel-1)
* Local: Ubuntu 24.04, R 4.3.3

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for the reviewer

* The dataset `coral_trout` is derived from Australian Institute of Marine
  Science data published under CC BY 3.0 Australia. Its citation, attribution
  and the changes made are given in `LICENSE.note` and in `?coral_trout`.
* The vignette loads a saved result (`inst/extdata/coral_trout_influence.rds`)
  for one step that takes about three minutes; the script that makes it is in
  `data-raw/` in the source repository.
