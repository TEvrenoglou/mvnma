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
