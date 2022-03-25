# R package: hlt <img src="man/figures/logo.png" align="right" width="120" />

### Higher-Order Item Response Theory (Latent Trait Theory) with Regression

Higher-order latent trait theory (item response theory). We implement the generalized partial credit model with a second-order latent trait structure. Latent regression can be done on the second-order latent trait.

---

### Resources

* [Learn about MCMC](https://m-clark.github.io/docs/ld_mcmc/index_onepage.html#preface) (MCMC book)
* [Ask a question/ Open an issue](https://github.com/mkleinsa/hlt/issues) (GitHub issues for bug reports, feature requests)

### Installation

#### Latest Release

The most recent **hlt** release can be installed from CRAN via

```r
install.packages("hlt")
```

#### Development Version

The most recent **hlt** development release can be installed from Github via devtools. If you do not have devtools installed on your system, go here to install it.

```r
install.packages("devtools")
```

```r
devtools::install_github("mkleinsa/hlt")
```

### Help/Getting started

Once installed, load the package with 

```r
library("hlt")
```

#### Examples from `?hlt`

Type `?hlt` for examples and main function documentation.

Here we analyze the Adult Self-Transcendence Inventory data from the R package `MPsychoR`.

Description of the ASTI data:

Adult Self-Transcendence Inventory scale which measures wisdom using five dimensions. The five dimensions are: self-knowledge and integration (SI), peace of mind (PM), non-attachment (NA), self-transcendence (ST), and presence in the here-and-now and growth (PG).

We are interested in extracting a single general latent measure of wisdom and doing regression on the wisdom dimension.

Levenson, M. R., Jennings, P. A., Aldwin, C. M., & Shiraishi, R. W. (2005). Self-transcendence: conceptualization and measurement. The International Journal of Aging and Human Development, 60, 127-143.

```r
# load the package
library("hlt")

# set seed
set.seed(153)

# load the asti data set
data("asti")

# shift responses to range from 0 instead of 1
x = as.matrix(asti[, 1:25]) - 1

# subset and transform predictor data
z = asti[, 26:27]
z[, 1] = (z[, 1] == "students") * 1
z[, 2] = (z[, 2] == "male") * 1
z = as.matrix(z)

# specify which items from x belong to each domain
id = c(0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4)

# fit the model
mod = hlt(x, z = z, id = id, iter = 2.2e6, burn = 2e6, delta = 0.013)

mod$accept.rate # ideally 0.234

# inspect the results
plot(mod, "lambda1", type = "trace")
plot(mod, "lambda2", type = "trace")
plot(mod, "lambda3", type = "trace")
plot(mod, "a1", type = "trace")
plot(mod, "d2", type = "trace")
plot(mod, "beta1", type = "trace")

plot(mod, 7, type = "icc", min = -10, max = 10)

summary(mod, param = "beta")
summary(mod, param = "lambda")
summary(mod, param = "alpha")
summary(mod, param = "delta")
summary(mod, param = "theta", dimension = 1)
summary(mod, param = "theta", dimension = 2)
summary(mod, param = "theta", dimension = 3)
summary(mod, param = "theta", dimension = 4)
```

### Contributing 

If you are interested in contributing to the development of **hlt** please open an issue to request.

### References and other literature

Paper under review.
