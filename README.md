# R package: hlt <img src="man/figures/logo.png" align="right" width="120" />

### Latent Regression in Higher-Order Item Response Theory with the R Package hlt

Item response theory (IRT) has become a standard method in the analysis of survey assessment data. In the IRT paradigm, the performance of test takers (i.e., their capacity to answer a question correctly) is explained by individual ability and the properties or characteristics of the test. Classically, ability is assumed to be 1) univariate and 2) homogenous across person characteristics. First, ability, which is a placement of each test taker on a continuum, can often be explained by heterogeneity among study participants. Second, recent applications involve characterizing a general, higher-order, ability which explain multiple, first-order, domains of ability. The natural progression, then, is to explain the heterogeneity of general ability. The open-source R package hlt implements the random-walk Metropolis-Hastings algorithm to estimate the higher-order IRT model by implementing a flexible Bayesian framework. We implement a higher-order generalize partial credit model and its extension of latent regression with the goal of explaining the relationship between the general latent construct and a set of explanatory variables.

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
summary(mod, param = "cor.theta")
summary(mod, param = "alpha")
summary(mod, param = "delta")
summary(mod, param = "theta", dimension = 1)
summary(mod, param = "theta", dimension = 2)
summary(mod, param = "theta", dimension = 3)
summary(mod, param = "theta", dimension = 4)
```

### Extensive examples
<details>
  <summary>Extra examples from the old package documentation</summary>
  
```r

set.seed(7794)
xdat = hltsim(n = 150, ntheta = ntheta, lambda = lambda, 
              tJ = id, dL = dL, mua = 1, mud = 1.4, siga = 0.8, sigd = 1,
              beta = 1)
x = xdat$x
mod = hlt(x, ntheta = 3, id = id, iter = 40000, burn = 30000, delta = 0.05)
mod$accept.rate
post = mod$post

Regression example 
set.seed(2034)
nB = 1
n = 300
z = matrix(sample(0:1, n, replace = TRUE), nrow = n, ncol = nB)
xdat = hltsim(n = n, ntheta = 3, lambda = c(0.8, 0.4, 0.3), 
              tJ = id, dL = 4, mua = 1, mud = 1.4, siga = 0.8, sigd = 1,
              regression = TRUE, z = z, nB = nB, beta = 1)
apply(xdat$x, 2, table)
lm(xdat$theta[,4] ~ xdat$z)
xdat$s.beta
x = xdat$x
x[,1] = ifelse(x[,1] == 3, 2, x[,1])
x[,1] = ifelse(x[,1] == 2, 1, x[,1])
x[,4] = ifelse(x[,4] == 3, 2, x[,4])
x[,5] = ifelse(x[,5] == 3, 2, x[,5])
mod = hlt(x, z = z, id = id, iter = 1e5, delta = 0.017)
mod = hlt(x, z = z, id = id, iter = 1e6, burn = 9e5, delta = 0.017)

mod = hlt(x, z = z, id = id, iter = 2.2e6, burn = 2e6, delta = 0.013)
mod$accept.rate
post = mod$post
apply(post, 2, mean)

smy = function(x) {c(mean = mean(x), se = sd(x), quantile(x, probs = c(0.025, 0.5, 0.7, 0.975)))}
apply(post[, "beta1", drop = FALSE], 2, smy)

summary(mod, param = "all")
summary(mod, param = "beta")
summary(mod, param = "lambda")
summary(mod, param = "alpha")
summary(mod, param = "delta")
summary(mod, param = "theta", dimension = 1)
summary(mod, param = "theta", dimension = 2)
summary(mod, param = "theta", dimension = 3)
summary(mod, param = "theta", dimension = 4)

cor(as.vector(xdat$s.delta[-1,]), 
    as.vector(matrix(summary(mod, param = "delta")[1,], nrow = 3)))
    
cor(summary(mod, param = "alpha")[1,], xdat$s.alpha)

summary(mod, param = "cor.theta", cor.theta = c(1,2))
xdat$s.lam.cor

th = summary(mod, param = "theta", dimension = 1)
cor(th, xdat$theta[,1])

plot(mod, "lambda1", type = "trace")
plot(mod, "lambda2", type = "trace")
plot(mod, "lambda3", type = "trace")
plot(mod, "a1", type = "trace")
plot(mod, "d2", type = "trace")
plot(mod, "beta1", type = "trace")

plot(mod, 1, type = "icc")
plot(mod, 2, type = "icc")
plot(mod, 3, type = "icc")
plot(mod, 4, type = "icc")
plot(mod, 5, type = "icc")
plot(mod, 6, type = "icc")
plot(mod, 7, type = "icc", min = -10, max = 10)

example 
data("asti")
x = as.matrix(asti[, 1:25]) - 1
z = asti[, 26:27]
z[, 1] = (z[, 1] == "students") * 1
z[, 2] = (z[, 2] == "male") * 1
z = as.matrix(z)
id = c(0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4)
mod = hlt(x, z = z, id = id, iter = 2.2e6, burn = 2e6, delta = 0.013)


##Test each combination of input values ###
n = 200
nB = 5
z = matrix(sample(0:1, n, replace = TRUE), nrow = n, ncol = nB)

PCM, No regression, 2 dimensions
xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
              tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
              dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
              regression = FALSE, z = NULL)
apply(xdat$x, 2, table)
mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
          iter = 200, burn = 100, delta = 0.013)

PCM, Regression, 2 dimensions
xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
              tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
              dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
              regression = FALSE, z = z)
apply(xdat$x, 2, table)
mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
          iter = 200, burn = 100, delta = 0.013)

PCM, No regression, 4 dimensions
xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
              tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
                     3,3,3,3,3,3),
              dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
              regression = FALSE, z = NULL)
apply(xdat$x, 2, table)
mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
          3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)

PCM, Regression, 4 dimensions
xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
              tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
                     3,3,3,3,3,3),
              dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
              regression = FALSE, z = z)
apply(xdat$x, 2, table)
mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
          3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)






GPCM, No regression, 2 dimensions
xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
              tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
              dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
              regression = FALSE, z = NULL)
apply(xdat$x, 2, table)
mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
          iter = 200, burn = 100, delta = 0.013)

GPCM, Regression, 2 dimensions
xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
              tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
              dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
              regression = FALSE, z = z)
apply(xdat$x, 2, table)
mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
          iter = 200, burn = 100, delta = 0.013)

GPCM, No regression, 4 dimensions
xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
              tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
                     3,3,3,3,3,3),
              dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
              regression = FALSE, z = NULL)
apply(xdat$x, 2, table)
mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
          3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)

GPCM, Regression, 4 dimensions
xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
              tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
                     3,3,3,3,3,3),
              dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
              regression = FALSE, z = z)
apply(xdat$x, 2, table)
mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
          3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)

}

```

</details>

### Contributing 

If you are interested in contributing to the development of **hlt** please open an issue to request.

### References and other literature

Paper under review.
