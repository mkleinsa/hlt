# R package: hlt <img src="man/figures/logo.png" align="right" width="120" />

<!-- badges: start -->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rstanarm?color=blue)](http://cran.r-project.org/package=rstanarm)
[![Downloads](http://cranlogs.r-pkg.org/badges/rstanarm?color=blue)](http://cran.rstudio.com/package=rstanarm)
[![R-CMD-check](https://github.com/stan-dev/rstanarm/workflows/R-CMD-check/badge.svg)](https://github.com/stan-dev/rstanarm/actions)
<!-- badges: end -->

### Flexible Item Response Theory (IRT) for Complex Survey Design

This is an R package for estimation of item response theory models (IRT) when
the user needs (1) flexable specifiaction of univariate and multivariate (in a variety of forms) IRT models
for (2) dichotemous and polytonomous items responses using a wide array of likelihoods,
and/or (3) complex survey desing adjustments (cluster, strata, and weights).

The package was designed with a simple user interface in mind for the applied researcher. 
It is flexible enough to accomodate any IRT model with univariate and multivariate latent
trait structures. If the survey sampling accounted for the population of interest, the 
user can do population level inference by accounting for complex survey design. If
you are interesting in a convenient interface for fitting flexible multivariate item 
response theory models without complex survey adjustments, then this is also the package
for you.

Click the arrows for more details:

<details><summary>Quick start guide</summary>

Hello.
</details>

<details><summary>More details</summary>

Hello.
</details>

<details><summary>Modeling</summary>

* [__`stan_jm`__](https://mc-stan.org/rstanarm/reference/stan_jm.html)

   Estimates shared parameter joint models for longitudinal and time-to-event
   (i.e. survival) data. The joint model can be univariate (i.e. one longitudinal
   outcome) or multivariate (i.e. more than one longitudinal outcome). A variety
   of parameterisations are available for linking the longitudinal and event
   processes (i.e. a variety of association structures).

</details>

<details><summary>Estimation</summary>
Hi again
</details>

---

### Resources

* [mc-stan.org/rstanarm](https://mc-stan.org/rstanarm) (online documentation, vignettes)
* [Ask a question](https://discourse.mc-stan.org) (Stan Forums on Discourse)
* [Open an issue](https://github.com/stan-dev/rstanarm/issues) (GitHub issues for bug reports, feature requests)

### Installation

#### Latest Release

The most recent **rstanarm** release can be installed from CRAN via

```r
install.packages("rstanarm")
```

#### Development Version

To install from GitHub, first make sure that you can install the **rstan**
package and C++ toolchain by following these
[instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).
Once **rstan** is successfully installed, you can install **rstanarm** from
GitHub using the **remotes** package by executing the following in R:

```r
Sys.setenv(MAKEFLAGS = "-j4") # change 4 to however many cores you can/want to use to parallelize install 
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github("stan-dev/rstanarm", INSTALL_opts = "--no-multiarch", force = TRUE)
```

You can switch `build_vignettes` to `TRUE` but it takes a lot longer to install and the 
vignettes are already separately available from the 
[Stan website](http://mc-stan.org/rstanarm/articles/index.html) 
and 
[CRAN](https://cran.r-project.org/package=rstanarm/vignettes). 
If installation fails, please let us know by [filing an issue](https://github.com/stan-dev/rstanarm/issues).

#### Survival Analysis Version

The `feature/survival` branch on GitHub contains a development version of **rstanarm** that includes survival analysis functionality (via the `stan_surv` modelling function). Until this functionality is available in the CRAN release of **rstanarm**, users who wish to use the survival analysis functionality can install a binary version of the survival branch of **rstanarm** from the Stan R packages repository with:

```
install.packages("rstanarm", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

Note that this binary is static (i.e. it is not automatically updated) and is only hosted so that users can access the (experimental) survival analysis functionality without needing to go through the time consuming (and sometimes painful) task of installing the development version of **rstanarm** from source.

### Contributing 

If you are interested in contributing to the development of **rstanarm** please 
see the [developer notes](http://mc-stan.org/rstanarm/dev-notes/index.html) page.
