# gamfluence 0.1.1

First release.

* `gam_influence()` deletes each site (cluster) of a `gam()` fit in turn and
  refits the model twice, with the estimated smoothing parameters (and
  negative binomial theta) held and re-estimated. It reports the fixed, added
  and total changes in the other sites' fitted values, and a table of
  hyperparameter changes.
* The re-estimated refit is searched from two starting points and the better
  REML fit is kept; `status` flags sites with two clearly different optima.
* `plot()` shows each site's fixed change against its added change.
* A scope guard rejects models outside the supported class, and a self-check
  confirms that the refits reproduce the original fit.
* Dataset `coral_trout`: 467 surveys at 71 inshore reef sites of the Great
  Barrier Reef, derived from Australian Institute of Marine Science data under
  CC BY 3.0 Australia (see `LICENSE.note`).
* Vignette: "Which sites change a GAM, and why? A coral trout example".
