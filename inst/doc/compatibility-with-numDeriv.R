## ----include = FALSE----------------------------------------------------------
knitr::knit_hooks$set(pngquant = knitr::hook_pngquant)
knitr::opts_chunk$set(
  dev = "png",
  dev.args = list(type = if (Sys.info()["sysname"] == "Darwin") "quartz" else "cairo-png"),
  fig.width = 512 / 72,
  fig.height = 320 / 72,
  out.width="10cm",
  dpi = 72,
  fig.retina = 1,
  collapse = TRUE,
  comment = "#>"
)
if (.Platform$OS.type == "unix") knitr::opts_chunk$set(pngquant = "--speed=1 --quality=50-60")

## ----setup--------------------------------------------------------------------
library(pnd)
library(numDeriv)

## -----------------------------------------------------------------------------
grad(sin, pi)  # Old
Grad(sin, pi)  # New

## -----------------------------------------------------------------------------
grad(sin, (0:10)*2*pi/10)  # Old
Grad(sin, (0:10)*2*pi/10)  # New

## ----message=FALSE------------------------------------------------------------
func0 <- function(x) sum(sin(x))
grad(func0, (0:10)*2*pi/10)  # Old
Grad(func0, (0:10)*2*pi/10)  # New

## ----message=FALSE, R.options = list(digits = 4)------------------------------
func1 <- function(x) sin(10*x) - exp(-x)
curve(func1, from = 0, to = 5)
x <- 2.04
numd1 <- grad(func1, x)  # Old
numd2 <- Grad(func1, x)  # New
numd3 <- Grad(func1, x, h = gradstep(func1, x)$par)  # New auto-selection
exact <- 10*cos(10*x) + exp(-x)
c(Exact = exact, Old = numd1, New = numd2, NewAuto = numd2,
  OldErr = (numd1-exact)/exact, NewErr = (numd2-exact)/exact,
  NewAutoErr = (numd3-exact)/exact)

## ----message=FALSE, R.options = list(digits = 4)------------------------------
x <- c(1:10)
numd1 <- grad(func1, x)  # Old
numd2 <- Grad(func1, x)  # New
numd3 <- Grad(func1, x, h = sapply(x, function(y) gradstep(func1, y)$par))
exact <- 10*cos(10*x) + exp(-x)
cbind(Exact = exact, Old = numd1, New = numd2, NewAuto = numd2,
 OldErr = (numd1-exact)/exact, NewErr = (numd2-exact)/exact,
  NewAutoErr = (numd3-exact)/exact)

## -----------------------------------------------------------------------------
func2 <- function(x) c(sin(x), cos(x))
x <- (0:1)*2*pi
jacobian(func2, x)
Jacobian(func2, x)

## -----------------------------------------------------------------------------
x <- 0.25 * pi
hessian(sin, x)
fun1e <- function(x) sum(exp(2*x))
x <- c(1, 3, 5)
hessian(fun1e, x, method.args=list(d=0.01))

## -----------------------------------------------------------------------------
system.time(replicate(1e5, sin(rnorm(10))))
system.time(sin(rnorm(1e6)))

## ----examplevec, error = TRUE-------------------------------------------------
try({
f <- function(x) quantile(x, 1:3/4)
grad(f, x = 1:2)
grad(f, x = 1:4)
grad(f, x = 1:3)
})

## -----------------------------------------------------------------------------
jacobian(f, x = 1:3)
jacobian(f, x = 1:4)

## ----sincos-------------------------------------------------------------------
f <- function(x) c(sin(x), cos(x))
f(1)
jacobian(f, 1:2)
jacobian(f, 1:3)

## ----error = TRUE-------------------------------------------------------------
try({
f2 <- function(x) c(sin(x), cos(x))  # Vector output -> gradient is unsupported
grad(f2, x = 1:4)
hessian(f2, x = 1:4)

Grad(f2, x = 1:4)
})

## ----error = TRUE-------------------------------------------------------------
try({
f2 <- function(x) c(sum(sin(x)), sum(cos(x)))
grad(f2, x = 1:4)
hessian(f2, x = 1:4)
jacobian(f2, x = 1:4)

Grad(f2, x = 1:4)
Jacobian(f2, x = 1:4)
})

## ----error = TRUE-------------------------------------------------------------
try({
f2 <- function(x) c(sin(x), cos(x))
grad(f2, x = 1:4)
jacobian(f2, x = 1:4)
})

## -----------------------------------------------------------------------------
f <- function(x) {cat(x, " "); sin(x)}
grad(f, 1.7e-5, method.args = list(r = 2))  # step 1.0e-4
grad(f, 1.8e-5, method.args = list(r = 2))  # step 1.8e-9

## -----------------------------------------------------------------------------
xseq <- 10^seq(-10, 2, 0.25)
sseq2 <- stepx(xseq, acc.order = 2)
sseq4 <- stepx(xseq, acc.order = 4)
matplot(xseq, cbind(sseq2, sseq4), lty = 1:2, col = 1, type = "l", bty = "n",
        log = "xy", main = "Default step size", xlab = "Point", ylab = "Step")
legend("topleft", c("2nd-order accurate", "4th-order accurate"), lty = 1:2)

## ----richardsonprint, message=FALSE-------------------------------------------
f <- function(x) {cat(x, " "); sin(x)}
x0 <- 1
g1 <- numDeriv::grad(f, x0)
print(g1)

cat("Auto-detected step:", step.SW(sin, x0)$par, "\n")
hgrid <- 10^seq(-10, -4, 1/32)
errors <- sapply(hgrid, function(h) Grad(sin, x0, h = h, cores = 1,
            elementwise = TRUE, vectorised = TRUE, multivalued = FALSE)) - cos(x0)
plot(hgrid, abs(errors), log = "xy", cex = 0.6)

## -----------------------------------------------------------------------------
b <- fdCoef(stencil = c(-(2^(3:0)), 2^(0:3)))
print(b)

## ----stencil------------------------------------------------------------------
fd <- sin(x0 + b$stencil / 8 * 1e-4) * b$weights
abs(fd[1:4]) / sum(abs(fd[1:4]))

## ----richardson---------------------------------------------------------------
g2 <- Grad(f, x0, h = 1.25e-05, acc.order = 4,
           elementwise = TRUE, vectorised = TRUE, multivalued = FALSE)
print(g2)

c(diff = g1 - g2, Error8 = cos(x0) - g1, Error4 = cos(x0) - g2)

## -----------------------------------------------------------------------------
f <- function(x) {print(x, digits = 16); x^9} 
fp1 <- numDeriv::grad(f, x = 1, method.args = list(r = 4, d = 2^-10, show.details = TRUE))
print(fp1, digits = 16)

## -----------------------------------------------------------------------------
# f <- function(x) x^9
# fp2 <- pnd::Grad(f, x = 1, h = "SW", acc.order = 8, vectorised1d = TRUE)
# print(attributes(fp2)$step.search$iterations, digits = 16)

## -----------------------------------------------------------------------------
g <- function(x) sum(sin(x))
hessian(g, 1:3, method.args = list(show.details = TRUE))

