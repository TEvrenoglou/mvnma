## mvnma, version 0.4-0 (2026-mm-dd)

### User-visible changes

* Metadata for `mvnma`, `mvrank`, and `netsplit.mvnma` objects is now stored in
  named list elements rather than attributes.

* `mvrank` objects now contain separate `ranks` and `ranks.shared` components,
  together with `trts` and `trts.shared`.

* print.mvrank():
  - uses clearer labels for ranking summaries: `P(best)`, `Median rank`,
    `Mean rank`, and `95% CrI for rank`
  - combines the lower and upper credible limits into a single displayed
    credible-interval column.

### Internal changes

* Added migration support for objects created with previous package versions.

* Refactored downstream functions to use the new result-object structures.


## mvnma, version 0.3-0 (2026-09-17)

### Major changes

* mvnma():
  - accepts one or more outcomes and studies with two or more arms
  - a common-effects model can be fitted
  - between-study heterogeneity parameter(s) can be prespecified

* New functions netsplit.mvnma(), print.netsplit.mvnma(), and
  forest.netsplit.mvnma() support consistency checks using the node-splitting
  method.

### User-visible changes

* mvnma():
  - new argument 'pooled' to fit a common-effects or random-effects model
  - new argument 'psi.preset' to prespecify the between-study heterogeneity
    parameter(s)


## mvnma, version 0.2-0 (2026-07-09)

### Major changes

* mvnma():
  - default for argument 'reference.group' extracted from pairwise() objects
    (if it is identical across all objects)

### User-visible changes

* mvnma():
  - new argument 'varTE.missing' to specify the variance for outcomes not
    reported in a study

### Bug fixes

* mvnma():
  - fix within-study variance-covariance matrix for missing outcomes
  - fix the error "Error: connections left open" due to not closing the
    connection with the model code

## mvnma, version 0.1-0 (2026-05-15)

* initial release on CRAN
