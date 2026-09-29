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
* `plot()` method: each site's `fixed` change against its `added` change, with
  the diagonal where they are equal, the top sites labelled, and sites with two
  REML optima marked.
* Works with mgcv 1.9-4, where `nb()` and `negbin()` read the link with
  `substitute()`: the refits now pass the link as a literal string. Before
  this, every negative binomial model with estimated theta failed at setup.
* New dataset `coral_trout`: 467 surveys of coral trout at 71 inshore reef
  sites (Palm and Whitsunday island groups, 2007-2018), derived from AIMS data
  under CC BY 3.0 AU (see `LICENSE.note`). `data-raw/coral_trout.R` rebuilds it.
