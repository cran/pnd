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
hseq <- 10^seq(-9, -3, length.out = 101)
hopt <- (1.5 * .Machine$double.eps)^(1/3)
e <- function(h) h^2/6 + 0.5*.Machine$double.eps / h
plot(hseq, e(hseq), log = "xy", xlab = "Step size", ylab = "Total error", asp = 1,
     bty = "n", type = "l", main = "Inaccuracy of CD-based derivatives (logarithmic)")
abline(v = hopt, lty = 2)
hseq2 <- seq(hopt - 5e-6, hopt + 5e-6, length.out = 101)
plot(hseq2, e(hseq2), xlab = "Step size", ylab = "Total error", bty = "n",
     type = "l", main = "Inaccuracy of CD-based derivatives (linear)")
abline(v = hopt, lty = 2)

## -----------------------------------------------------------------------------
fdCoef(deriv.order = 2, acc.order = 2) # Same as fdCoef(2)

## -----------------------------------------------------------------------------
w3 <- fdCoef(3)
w4 <- fdCoef(4)
print(w3)
print(w4)

## -----------------------------------------------------------------------------
denom  <- factorial(0:6)
names(denom) <- paste0("f'", 0:6)
num.0  <- c(1, rep(0, 6)) # f(x) = f(x) + 0*f'(x) + 0*f''(x) + ...
num.h  <- rep(1, 7)
num.2h <- 2^(0:6)
# f(x-ch) = f - ch f' + (ch)^2/2 f'' - (ch)^3/6 f''' ...
num.mh  <- suppressWarnings(num.h * c(1, -1)) # Relying on recycling
num.m2h <- suppressWarnings(num.2h * c(1, -1))
num <- colSums(rbind(num.m2h, num.mh, num.0, num.h, num.2h) * w4$weights)
print(round(num / denom, 5))

## -----------------------------------------------------------------------------
sum(abs(w4$weights[c(1, 2, 4, 5)] *
          c(num.m2h[7], num.mh[7], num.h[7], num.2h[7]))) / denom[7]

## -----------------------------------------------------------------------------
sum(abs(w4$weights))

## -----------------------------------------------------------------------------
w <- fdCoef(acc.order = 4)
h2 <- 2^(0:5)
h  <- rep(1, 6)
hm <- h * c(1, -1)
h2m <- h2 * c(1, -1)
coef.tab <- rbind(h2m, hm, h, h2) # Here, using rbind is more convenient
rownames(coef.tab) <- names(w$weights)
colnames(coef.tab) <- paste0("f'", seq_len(ncol(coef.tab)) - 1)
print(coef.tab * w$weights)
print(colSums(coef.tab * w$weights))

## -----------------------------------------------------------------------------
0.5*sum(abs(w$weights))

## -----------------------------------------------------------------------------
m <- 7  # Example order
w <- fdCoef(m)
k <- sum(w$stencil > 0) # ±h, ±2h, ..., ±kh
num.pos <- sapply(1:k, function(i) i^(0:(m+2)))
num.neg <- apply(num.pos, 2, function(x) x * c(1, -1))
num.neg <- num.neg[, rev(seq_len(ncol(num.neg)))]
nz <- abs(w$stencil) > 1e-12 # Non-zero function weights
coef.tab <- sweep(cbind(num.neg, num.pos), 2, w$weights[nz], "*")
rownames(coef.tab) <- paste0("f'", seq_len(nrow(coef.tab))-1)
colnames(coef.tab) <- names(w$weights[nz])
print(coef.tab)

## -----------------------------------------------------------------------------
rs <- rowSums(coef.tab)
print(round(rs, 4))

## -----------------------------------------------------------------------------
print(c1 <- sum(abs(coef.tab[nrow(coef.tab), ])) / factorial(m+2))

## -----------------------------------------------------------------------------
print(c2 <- 0.5*sum(abs(w$weights)))

## -----------------------------------------------------------------------------
getCoefs <- function(x, ord) do.call(expand.grid, replicate(ord, x, simplify = FALSE))
rowProd <- function(x) apply(x, 1, prod)
getMults <- function(ord) {
  steps  <- list(c(1, 1), c(-1, 1), c(1, -1), c(-1, -1))
  coefs <- lapply(steps, getCoefs, ord = ord)
  signs <- lapply(coefs, rowProd)
  mults <- signs[[1]] - signs[[2]] - signs[[3]] + signs[[4]]
  ntab <- expand.grid(replicate(ord, c("i", "j"), simplify = FALSE))
  names(mults) <- apply(ntab, 1, paste0, collapse = "")
  return(mults)
}

print(getMults(2))
print(getMults(3))

## -----------------------------------------------------------------------------
mults4 <- getMults(4)
print(mults4[mults4 != 0])

## -----------------------------------------------------------------------------
x   <- -0.2456605107847454  # 16 sig. digs
t   <- 1/59
print(res1 <- (x-t)^4 * (t^-4 / 4), 17)
print(res2 <- (x-t)^4 / (t^4 * 4), 17)
print(res1 - res2)

## -----------------------------------------------------------------------------
f1 <- function(x) expm1(x)^2 + (1/sqrt(1+x^2) - 1)^2
df1 <- function(x) 2 * exp(x) * expm1(x) - 2*x * (1/sqrt(1+x^2) - 1) / (1 + x^2) / sqrt(1 + x^2)
ddf1 <- function(x) (6 * (1/sqrt(1 + x^2) - 1) * x^2)/(1 + x^2)^(5/2) + (2 * x^2)/(1 + x^2)^3 - (2 * (1/sqrt(1 + x^2) - 1))/(1 + x^2)^(3/2) + 2 * exp(2*x) + 2 * exp(x) * expm1(x)
dddf1 <- function(x) (18 * (1/sqrt(1 + x^2) - 1) * x)/(1 + x^2)^(5/2) + (6 * x)/(1 + x^2)^3 - (30 * (1/sqrt(1 + x^2) - 1) * x^3)/(1 + x^2)^(7/2) - (18 * x^3)/(1 + x^2)^4 + 6 * exp(2*x) + 2 * exp(x) * expm1(1)

squish <- function(x, pow = 1/2, shift = 1) ((abs(x) + shift)^pow - shift^pow) * sign(x)
unsquish <- function(y, pow = 1/2, shift = 1) ((abs(y) + shift^pow)^(1/pow) - shift) * sign(y)

xgrid <- seq(-0.2, 1.5, 0.01)
par(mar = c(2, 2, 0, 0) + .1)
plot(xgrid, squish(f1(xgrid)), type = "l", ylim = squish(c(-1, 20)), bty = "n", lwd = 2, xlab = "", ylab = "", yaxt = "n")
lines(xgrid, squish(df1(xgrid)), col = 2)
lines(xgrid, squish(ddf1(xgrid)), col = 3)
lines(xgrid, squish(dddf1(xgrid)), col = 4)
axis(2, squish(ats <- c(-1, 0, 1, 2, 5, 10, 20)), labels = ats, las = 1)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  f1 <- function(x) x^4
#  try1.good  <- step.SW(x = 1, f1) # Starts at h0 = 1e-5
#  try1.small <- step.SW(x = 1, f1, h0 = 1e-9)
#  try1.large <- step.SW(x = 1, f1, h0 = 10)
#  try1.huge <- step.SW(x = 1, f1, h0 = 1000)
#  try1 <- list(try1.good, try1.small, try1.large, try1.huge)
#  
#  hvals1 <- sapply(try1, function(x) x$iterations$h[1])
#  for (i in 1:4) {
#    cat("\n\nDiagnostics for the SW79 algorithm, f(x) = x^4, x = 1, h0 =", hvals1[i], "\n")
#    cairo_pdf(paste0("power-", i, ".pdf"), 6.2, 4)
#    tryCatch(par(family = "Fira Sans"), error = \(e) NULL, warning = \(w) NULL)
#    printDiag(try1[[i]], true.val = 4, main = paste0("f(x) = x^4 with initial h = ", hvals1[i]))
#    dev.off()
#  }
#  
#  f2 <- sin
#  try2.good  <- step.SW(x = pi/4, f2) # Starts at h0 = 1e-5
#  try2.small <- step.SW(x = pi/4, f2, h0 = 1e-9)
#  try2.large <- step.SW(x = pi/4, f2, h0 = 10)
#  try2.huge <- step.SW(x = pi/4, f2, h0 = 1000) # Observe the custom warning
#  try2 <- list(try2.good, try2.small, try2.large, try2.huge)
#  
#  hvals2 <- sapply(try1, function(x) x$iterations$h[1])
#  for (i in 1:4) {
#    cat("\n\nDiagnostics for the SW79 algorithm, f(x) = sin(x), x = pi/4, h0 =", hvals1[i], "\n")
#    cairo_pdf(paste0("sine-", i, ".pdf"), 6.2, 4)
#    tryCatch(par(family = "Fira Sans"), error = \(e) NULL, warning = \(w) NULL)
#    printDiag(try2[[i]], true.val = sqrt(2)/2, main = paste0("f(x) = sin(x) with initial h = ", hvals2[i]))
#    dev.off()
#  }

