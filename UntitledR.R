library(devtools)
library(knitr)
library(roxygen2)
library(testthat)
 devtools::has_devel()
usethis::use_mit_license("Ziyue Li")
usethis::use_package("mclust", type = "Imports")
usethis::use_readme_md()
usethis::use_readme_rmd()

devtools::load_all()
devtools::document()


usethis::use_git()

library(gitcreds)

# link local repository with Github
# slide 56-57 of L4













