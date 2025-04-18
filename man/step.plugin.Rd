% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step-select.R
\name{step.plugin}
\alias{step.plugin}
\title{Plug-in step selection}
\usage{
step.plugin(
  FUN,
  x,
  h0 = 1e-05 * max(abs(x), sqrt(.Machine$double.eps)),
  range = h0/c(10000, 1e-04),
  cores = 1,
  preschedule = getOption("pnd.preschedule", TRUE),
  cl = NULL,
  diagnostics = FALSE,
  ...
)
}
\arguments{
\item{FUN}{Function for which the optimal numerical derivative step size is needed.}

\item{x}{Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.}

\item{h0}{Numeric scalar: initial step size, defaulting to a relative step of
slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}). This step
size for first derivarives is internallt translated into the initial step size for third
derivatives by multiplying it by the machine epsilon raised to the power -2/15.}

\item{range}{Numeric vector of length 2 defining the valid search range for the step size.}

\item{cores}{Integer specifying the number of CPU cores used for parallel computation.
Recommended to be set to the number of physical cores on the machine minus one.}

\item{preschedule}{Logical: if \code{TRUE}, disables pre-scheduling for \code{mclapply()}
or enables load balancing with \code{parLapplyLB()}. Recommended for functions that
take less than 0.1 s per evaluation.}

\item{cl}{An optional user-supplied \code{cluster} object  (created by \code{makeCluster}
or similar functions). If not \code{NULL}, the code uses \code{parLapply()} (if \code{preschedule}
is \code{TRUE}) or \code{parLapplyLB()} on that cluster on Windows, and \code{mclapply}
(fork cluster) on everything else.}

\item{diagnostics}{Logical: if \code{TRUE}, returns the full iteration history
including all function evaluations.
Note: the history tracks the third derivative, not the first.}

\item{...}{Passed to FUN.}
}
\value{
A list similar to the one returned by \code{optim()}: \code{par} -- the optimal
step size found, \code{value} -- the estimated numerical first derivative (central
differences), \code{counts} -- the number of iterations (here, it is 2),
\code{abs.error} -- an estimate of the total approximation error (sum of truncation and
rounding errors),
\code{exitcode} -- an integer code indicating the termination status:
\code{0} indicates termination with checks passed tolerance,
\code{1} means that the third derivative is exactly zero (large step size preferred),
\code{2} signals that the third derivative is too close to zero (large step size preferred),
\code{3} indicates a solution at the boundary of the allowed value range.
\code{message} is a summary message of the exit status.
If \code{diagnostics} is \code{TRUE}, \code{iterations} is a list
including the two-step size search path, argument grids, function values on those grids,
and estimated 3rd derivative values.
}
\description{
Plug-in step selection
}
\details{
This function computes the optimal step size for central differences using the
plug-in approach.
The optimal step size is determined as the minimiser of the total error, which for central
finite differences is (assuming minimal bounds for relative rounding errors)
\deqn{\sqrt[3]{1.5 \frac{f'(x)}{f'''(x) \epsilon_{\mathrm{mach}}}}}{[(1.5 mach.eps * f' / f''')^(1/3)]}
If the estimated third derivative is too small, the function assumes a third
derivative of 1 to prevent division-by-zero errors.
}
\examples{
f <- function(x) x^4
step.plugin(x = 2, f)
step.plugin(x = 0, f, diagnostics = TRUE)  # f''' = 0, setting a large one
}
\references{
\insertAllCited{}
}
