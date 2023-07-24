# update dependencies
devtools::install_deps()

# run unit tests
devtools::test() 

# run a local R CMD check
devtools::check()

# check for CRAN specific requirements
results = rhub::check_for_cran()
results$cran_summary()

# integrate comments into cran-comments.md
usethis::use_cran_comments()

# devtools::check_rhub()

# check win-builder
devtools::check_win_devel()

# spell check
devtools::spell_check()

# run goodpractice check
goodpractice::gp()

# submit to cran
# devtools::release()

# build pdf manual
devtools::build_manual()
