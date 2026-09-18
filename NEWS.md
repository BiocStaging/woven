# woven 0.99.3

* Removed a stale DESCRIPTION claim about sparse projection matrices (PMD);
  no such code path exists (`woven()` recovers dense projection matrices by
  design, optionally ridge-regularized and feature-screened via `ridge_w`
  and `screen_top`). DESCRIPTION now describes that mechanism instead.
* `.scale_fit()` (internal) rewritten without `<<-` (BiocCheck NOTE).
* `summary()`'s printed report now goes through a dedicated internal
  `wovenSummary` S4 class with its own `show` method, so `cat()` output
  only ever appears inside a `show` method (BiocCheck NOTE). `summary()`
  still returns the metrics vector invisibly; printed output is unchanged.
* Escaped literal braces in a `solver_mcca_dual.R` roxygen comment that
  `R CMD check` was parsing as Rd markup ("Lost braces" NOTE).
* `.Rbuildignore` now excludes `*.Rcheck`/`*.BiocCheck` directories, so a
  leftover local check-output directory can no longer ship inside the
  built source tarball (one had shipped in the 0.99.2 tarball).
* Added `tests/testthat/test-precompute-plot-summary.R`, covering
  `woven_precompute()` reuse/equivalence, a malformed-precomp error path,
  `summary()` with and without labels, `plot()` with/without labels and
  anchor-highlighting (including its out-of-range-dims error), and
  `woven_nystrom_error()` as a regression test for the 0.99.2 `Za_list`
  field-name bug fix. These close out the last open items from the
  agentic Bioconductor review's testing-gap finding.
* Verified clean from a rebuilt tarball: `R CMD check` 0 errors/0
  warnings/0 notes; `BiocCheck` 0 errors/0 warnings/4 notes (down from 6,
  all cosmetic); 88/88 tests pass.

# woven 0.99.2

* The object returned by `woven()` is now an S4 class (`woven`, defined in
  `AllClasses.R`) instead of a plain S3 list. `show`, `summary`, and `plot`
  are now real S4 methods. New accessor generics `Z()`, `WList()`,
  `anchorIdx()`, and `singularValues()` are the recommended way to pull the
  most commonly used slots out of a fit; see `?\`woven-class\`` for the
  complete slot list, or use `@` directly for anything an accessor does not
  cover.
* Fixed a real bug in `woven_nystrom_error()`: it referenced a field
  (`fit$Za_list`) that did not exist on the fit object (the correct name is
  `Z_anchors`), so the function always silently returned `NaN` rather than
  ever computing a real projection error. It also passed a too-short `Y`
  vector to the internal leave-one-out refit; both are fixed, and the
  function now returns a real, non-degenerate value.
* `plot()` on a woven fit now returns the ggplot object directly (visibly)
  instead of forcing a `print()` call internally and returning it invisibly;
  this is the standard R plot-method convention and avoids printing outside
  a `show` method.
* `woven()` now accepts a `MultiAssayExperiment` directly as `X_list`. Each
  experiment becomes one modality; subjects are aligned to `colData` via
  `sampleMap`, and a subject absent from a given experiment gets an all-NA
  row for it automatically, matching the existing block-missing convention.
  `Y` may be a colData column name (a single string) instead of an explicit
  vector. See `?woven` for a worked example and `.mae_to_xlist()` for the
  conversion internals.
* Laplacian construction in `woven_precompute()` now parallelizes via
  `BiocParallel` (`MulticoreParam`/`SnowParam`/`SerialParam` chosen by
  platform) instead of `parallel::mclapply`, which does not fork on Windows.
* Fixed a stale `?woven` description of solver routing (it described a
  since-removed V=2-only closed-form / V>=3 ALS split; the package has used
  a single closed-form dual SUMCOR MCCA solver for all V since 0.99.0).
* Fixed "Nystrm" typo (-> "Nystrom") in source comments and documentation.
* Vignette introduction rewritten: defines "assay block", explains the
  Nystrom method and its use as an out-of-sample extension (with a link to
  further reading), and states WOVEN's comparison to DIABLO concretely
  instead of via an unqualified "standard methods" claim.
* (Retroactive note, first available in 0.99.1) `woven()` exposes `ridge_w`
  and `screen_top`, ridge-regularized / feature-screened recovery of the
  projection matrices, useful when a modality has many more features than
  anchor subjects and the true signal is diffuse. Select jointly via
  cross-validation on held-out classification error.

# woven 0.99.1

* `woven()` now standardizes each modality to zero mean and unit variance by
  default (`scale = TRUE`), matching the internal scaling of comparator methods
  so that no modality dominates by measurement scale. Stored centers and scales
  are reapplied to new subjects in `woven_scores()` and `woven_predict()`. Set
  `scale = FALSE` to retain the previous raw-feature behavior.

# woven 0.99.0

* Initial Bioconductor submission.
* Unified dual SUMCOR MCCA solver (`woven_mcca_dual`) for all V >= 2 modalities.
* Nyström out-of-sample extension for block-missing subjects.
* Graph Laplacian regularization built on observed rows only (no imputation).
* Label-augmented cross-covariance for supervised integration.
* Full benchmark against DIABLO, MOFA2, IntegrAO, and Impute+DIABLO across four simulation arms and ADNI real-data validation.
