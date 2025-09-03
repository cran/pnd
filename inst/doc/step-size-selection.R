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

## -----------------------------------------------------------------------------
h0 <- 1e-4
x0 <- 1
fun <- function(x) sin(x)
getRatio <- function(FUN, x, h) {
  f  <- FUN(x)
  fplus  <- FUN(x + h)
  fminus <- FUN(x - h)
  fd <- (fplus - f) / h
  bd <- (f - fminus) / h
  cd <- (fplus - fminus) / h / 2
  et <- abs(cd - fd)
  er <- 0.5 * abs(f) * .Machine$double.eps / h
  ret <- et / er
  attr(ret, "vals") <- c(`h` = h,
                         bd = bd, cd = cd, fd = fd,
                         e_trunc = et, e_round = er)
  return(ret)
}
print(u0 <- getRatio(FUN = fun, x = x0, h = h0))  # 45035996, too high

## -----------------------------------------------------------------------------
h1 <- h0 * sqrt(100 / max(u0, 1))
print(u1 <- getRatio(FUN = fun, x = x0, h = h1))

## -----------------------------------------------------------------------------
uopt <- getRatio(FUN = fun, x = x0, h = (1.5 * tan(1) * .Machine$double.eps)^(1/3))
nd <- c(step0 = attr(u0, "vals")["cd"], step1 = attr(u1, "vals")["cd"],
        optimal = attr(uopt, "vals")["cd"])
print(total.err <- cos(x0) - nd)

## -----------------------------------------------------------------------------
h0 <- 1e-5
x0 <- 0.1
fun <- function(x) pi*x + exp(1)
print(u0 <- getRatio(FUN = fun, x = x0, h = h0))
h1 <- h0 * sqrt(100 / max(u0, 1))
print(u1 <- getRatio(FUN = fun, x = x0, h = h1))
h2 <- h1 * sqrt(100 / max(u0, 1))
print(u2 <- getRatio(FUN = fun, x = x0, h = h2))

## -----------------------------------------------------------------------------
h0 <- 2^-16
x0 <- sqrt(2)
fun  <- function(x) x^6 - 2*x^4 - 4*x^2
fun1 <- function(x) 6*x^5 - 8*x^3 - 8*x  # f'
fun3 <- function(x) 120*x^3 - 48*x       # f'''
print(u0 <- getRatio(FUN = fun, x = x0, h = h0))
h1 <- h0 * sqrt(100 / max(u0, 1))
print(u1 <- getRatio(FUN = fun, x = x0, h = h1))
hopt <- abs(1.5 * fun(x0) / fun3(x0) * .Machine$double.eps)^(1/3)
uopt <- getRatio(FUN = fun, x = x0, h = hopt)
fp.est <- c(step0 = attr(u0, "vals")["cd"], step1 = attr(u1, "vals")["cd"],
            optimal = attr(uopt, "vals")["cd"])
print(total.err <- fun1(x0) - fp.est)

## -----------------------------------------------------------------------------
dsin <- function(x, h) (sin(x+h) - sin(x-h)) / h / 2
totErr <- function(x, h, t = 0.5) (dsin(x, t*h) - dsin(x, h)) / (1 - t^2)
hgrid <- 2^seq(-52, 12, 0.5)
suppressWarnings(plot(hgrid, totErr(x = pi/4, hgrid), log = "xy",
     main = "Truncation + round-off error of d/dx sin(x) at x = pi/4",
     bty = "n", xlab = "Step size", ylab = "Sum of errors"))

## -----------------------------------------------------------------------------
h   <- c(2^-8, 2^-9)
te  <- totErr(x = pi/4, h = h)
print(diff(log(te)) / diff(log(h)))

## -----------------------------------------------------------------------------
h <- c(1e-4, 1.490116e-07)
b <- c(-h, 0, rev(h)) / min(h)
fc <- fdCoef(deriv.order = 1, stencil = b)
print(fc, 3)

## -----------------------------------------------------------------------------
cos(1) - sum(sin(1 + fc$stencil * h[2]) * fc$weights) / h[2]

## -----------------------------------------------------------------------------
hOpt <- function(x) (1.5 * abs(tan(x)) * .Machine$double.eps)^(1/3)

set.seed(1)
xgrid <- sort(runif(10000, max = 2*pi))
hgrid <- hOpt(xgrid)
df1 <- (sin(xgrid + hgrid) - sin(xgrid - hgrid)) / hgrid / 2
df2 <- (sin(xgrid + 7e-6) - sin(xgrid - 7e-6)) / 7e-6 / 2
fc  <- fdCoef(deriv.order = 1, stencil = c(-500, -1, 1, 500))
h4grid <- outer(hgrid, c(-500, -1, 1, 500))
df4 <- rowSums(sweep(sin(xgrid + h4grid), 2, fc$weights, "*")) / hgrid
df.true <- cos(xgrid)
err <- df.true - data.frame(df1, df2, df4)
abserr <- abs(err)

print(summary(err))
print(summary(abserr))

squish <- function(x, pow = 1/6, shift = 1e-12) ((abs(x) + shift)^pow - shift^pow) * sign(x)
xnice <- seq(0, 2*pi, length.out = 401)
plot(xnice, hOpt(xnice), type = "l", log = "y", main = "Optimal step size for sin(x)", bty = "n")

par(mar = c(2, 4, 0, 0) + .1)
matplot(xgrid, squish(abserr[, 2:3] - abserr[, 1]), type = "p", bty = "n",
        xlab = "", ylab = "", ylim = squish(c(-5e-10, 5e-11)), yaxt = "n",
        col = c("#0000AA88", "#AA000088"), pch = 16, cex = 0.3)
abline(h = 0, lwd = 3, col = "white")
abline(h = 0, lty = 2)
legend("bottomleft", c("Fixed vs. optimal", "4th-order vs. optimal"), bty = "n", pch = 16, col = c("#00008888", "#88000088"))
yax.vals <- c(-3e-10, -1e-10, -3e-11, -1e-11, -1e-12, 0, 1e-12, 1e-11, 3e-11)
axis(2, squish(yax.vals), yax.vals, las = 1)

## -----------------------------------------------------------------------------
S <- function(r) {
  if (is.unsorted(r)) r <- sort(r)
  n_2 <- length(r)
  n <- length(r) * 2
  s <- numeric(n_2)
  for (k in 1:n_2) {
    j <- setdiff(1:n_2, k)
    s[k] <- prod(abs(r[j]^2 / (r[k]^2 - r[j]^2))) / r[k]
  }
  prod(r)^(1/n_2) * sum(s)
} 

## -----------------------------------------------------------------------------
S(c(1, 2))
S(c(1, 3))    # Better
S(c(1, 2.6))  # Even better

## -----------------------------------------------------------------------------
ff <- Vectorize(function(x) S(c(1, x)))
oldpar <- par(mar = c(4, 4, 0, 0) + .1)
curve(ff, 1.2, 5, bty = "n", xlab = expression(c[1]), ylab = "S")

## -----------------------------------------------------------------------------
optim(par = 2, fn = ff, method = "BFGS")$par
# THE optimal length-4 grid is x0 + (-2.62, -1, 1, 2.62)h

# Finding the optimum for multiple lengths
nn <- 9
res <- vector("list", nn)
for (i in seq_along(res)) {
  set.seed(i)
  n <- 2 + 2*i
  # init.val <- (1:i)*2
  ff <- function(x) S(c(1, x))
  lower <- 1:i
  upper <- (1+(1:i)*3)^0.9
  init.pop <- sapply(1:i, function(j) runif(n*20, lower[j], upper[j]))
  init.pop <- apply(init.pop, 1, sort)
  if (NCOL(init.pop) == 1) init.pop <- matrix(init.pop) else 
    if (ncol(init.pop) > nrow(init.pop)) init.pop <- t(init.pop)
  ff0 <- apply(init.pop, 1, ff)
  x0 <- init.pop[which.min(ff0), ]
  res[[i]] <- sort(optim(par = x0, fn = ff, method = "BFGS",
                         control = list(reltol = 1e-8, maxit = 100))$par)
}

# Table 2b from Oliver & Ruffhead (1975)
tab <- matrix(nrow = length(res), ncol = nn)
for (i in 1:nn) tab[i, 1:i] <- res[[i]]

## -----------------------------------------------------------------------------
bmean <- round(colMeans(tab, na.rm = TRUE), 2)
print(bmean)
ff(bmean)  # Close to 10
ff(bmean[1])  # Close to 2

## -----------------------------------------------------------------------------
par(mar = c(4, 4, 0, 0) + .1)
plot(NULL, NULL, xlim = c(0, max(res[[i]]) + .5), ylim = c(0, nn+1),
     xlab = "Evaluation points", ylab = expression(n/2-1), bty = "n")
for (i in 1:9) points(c(1, res[[i]]), rep(i, length(res[[i]])+1), pch = 16, type = "b")
for (i in 1:9) points(sapply(res, "[", i), rep(1:nn), pch = 16, type = "b")
abline(v = 0, lty = 2)
abline(v = seq(1, max(bmean)+1), lty = 3, col = "#00000044")

## -----------------------------------------------------------------------------
b <- c(1, bmean)  # Nodes
x <- 1:10
plot(x, b, bty = "n", xlab = expression(n/2-1), ylim = c(0, max(bmean)+1),
     ylab = expression("Approximate"~c[i]))
ft <- nls(b ~ d + c*x^a, start = c(d = 0, c = 2.5, a = 0.75), weights = 1:10)
lines(x, -4.37 + 5.03 * x^0.55)
lines(x, predict(ft), col = 2)
par(oldpar)

