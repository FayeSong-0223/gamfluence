## Submission

This is a new submission.

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
