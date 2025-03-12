#
# (1) Make R packages available
#

library("devtools")
library("roxygen2")


#
# (2) Create documentation file(s)
#

document()


#
# (3) Build R package and PDF file with help pages
#

build()
#build(args = "--compact-vignettes=gs+qpdf")
build_manual()


#
# (4) Install R package
#

install()
#install(build_vignettes = TRUE)

check()

