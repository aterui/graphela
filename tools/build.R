
# build package -----------------------------------------------------------

usethis::use_mit_license(copyright_holder = "Akira Terui")
usethis::use_roxygen_md()
devtools::document()
devtools::load_all()
devtools::check(vignettes=FALSE)


# check syntax ------------------------------------------------------------

lintr::lint_package()


# test call ---------------------------------------------------------------

file.copy(
  "src/graphela.cpp",
  file.path(tempdir(), "graphela.cpp")
)

Rcpp::sourceCpp(file.path(tempdir(), "graphela.cpp"))

beta <- matrix(1, 3, 3)
diag(beta) <- 0
alpha <- c(1, 1, 1)
state <- c(1, 0, 1)
energy(state = state, alpha = alpha, beta = beta)
