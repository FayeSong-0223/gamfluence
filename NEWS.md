# gamfluence 0.0.1

* `gam_influence(model, cluster)`: for each site, refits the model without it
  with the estimated hyperparameters fixed and re-estimated, on the frozen
  full-data representation, and reports the fixed, added and total changes in
  the other sites' fitted values, plus a table of hyperparameter changes.
* Scope guard and setup self-check; row alignment to the rows mgcv used;
  rank-deficiency handling after deletion.
* Tests include independent references (closed form, ordinary refits, an
  independently coded REML criterion) and the identity fixed + added = total
  as a regression test. `validation/` holds the scripts behind the numbers.
* The re-estimated refit is searched from two starting points (mgcv's defaults
  and the full-fit values) and the fit with the better REML score is kept.
  `status` says `"two optima"` when the two end at clearly different fits.
* Random effects must be `s(g, bs = "re")` with a single factor `g`; random
  slopes are rejected.
