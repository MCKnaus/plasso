
# plasso <img src='docs/figures/plasso.png' align="right" height="139" />

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/plasso)](https://CRAN.R-project.org/package=plasso)
[![License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Downloads_Total](https://cranlogs.r-pkg.org/badges/grand-total/plasso)](https://CRAN.R-project.org/package=plasso)
[![Downloads_Monthly](https://cranlogs.r-pkg.org/badges/plasso)](https://CRAN.R-project.org/package=plasso)
[![Project_Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

Built on top of the `glmnet` library by Friedman, Hastie, and Tibshirani
(2010), the `plasso` package follows Knaus (2022) and comes up with two
functions that estimate least squares Lasso and Post-Lasso models. The
`plasso()` function adds coefficient paths for a Post-Lasso model to the
standard `glmnet()` output. On top of that `cv.plasso()` cross-validates
the coefficient paths for both the Lasso and Post-Lasso model and
provides optimal hyperparameter values for the penalty term lambda.

### Bug reports & support

For reporting a bug, simply [open an
issue](https://github.com/stefan-1997/plasso/issues/new) on GitHub. For
personal contact, you can write an email to
michael.knaus@uni-tuebingen.de.

### References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-glmnet" class="csl-entry">

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software* 33 (1): 1–22.
<https://doi.org/10.18637/jss.v033.i01>.

</div>

<div id="ref-knaus" class="csl-entry">

Knaus, Michael C. 2022. “<span class="nocase">Double machine
learning-based programme evaluation under unconfoundedness</span>.” *The
Econometrics Journal* 25 (3): 602–27.
<https://doi.org/10.1093/ectj/utac015>.

</div>

</div>
