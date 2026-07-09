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
