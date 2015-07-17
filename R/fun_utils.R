# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# COMPUTE PARTIAL CORRELATION -----------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Compute the Partial Correlation
#'
#' Computes the partial correlation of two variables X and Y controlling for a
#' third variable Z.
#'
#' @param x A matrix comprising 3 columns containing the variables X, Y, and Z
#' @return A single value, i.e., the partial correlation
#' @seealso Alternatives are, for example \code{\link[psych]{partial.r}}.
partial_cor <- function(x) {
  # x is a matrix of 3 columns containing the variables X, Y, and Z
  x <- cor(x)
  (x[1,2] - x[1,3] * x[2,3]) / (sqrt(1 - x[1,3]^2) * sqrt(1 - x[2,3]^2))
}


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# COMPUTE EFFECT SIZE F^2 ---------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Compute the effect size $f^2$ for nested (hierarchical) regression models
#'
#' detailed description
#'
#' @param model a model of class \code{lm}
#' @return A vector containing the $f^2$ values for every coefficient in the model
#' @seealso \code{\link[lm.support]{lm.sumSquares}}.
#' @export
f.squared <- function(model) {
  r2 <- numeric(length=model$rank - 1)
  for (i in 2:model$rank) {
    r2[i-1] <- summary(do.call("update",
                               args=list(model,
                                         paste(" ~ . -",
                                               names(model$coefficients)[i])))
    )$r.squared
  }
  f2 <- (summary(model)$r.squared - r2) / (1 - summary(model)$r.squared)
  names(f2) <- names(model$coefficients)[-1]
  return(f2)
}


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# PLOT ITEM CURVES ----------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Plot Category Response Curves for a polytomous item
#'
#' detailed description
#'
#' @param loc The overall item difficulty of the item
#' @return thres The thresholds of the different categories
#' @seealso
#' @export
plot_CRC <- function(loc, thres, alp = 1, ...) {

  K <- length(thres)
  categ <- K + 1

  x <- seq(-3 + 0, 3 + 0, length=100)
  p <- matrix(nrow=100, ncol=categ)

  B <- B <- matrix(0:K, ncol=1)
  A <- diag(1, K)
  A[lower.tri(A)] <- 1
  A <- rbind(0, A)
  A <- cbind(0:K, A)
  delta <- c(loc, thres)
  for(i in 1:100) {
    for(j in 1:categ) {
      p[i, j] <- (exp(alp*(B[j] * x[i] - A[j, ] %*% delta)))/
        (sum(exp(alp*(B %*% x[i] - A %*% delta))))
    }
  }

  ellipsis <- list(...)
  ellipsis <- ellipsis[!names(ellipsis) %in% "add"]

  if (all(any(names(list(...)) %in% "add"), list(...)$add == TRUE)) {
    do.call(lines, args =
              c(list(x = x, y = p[, 1], type = "l"), ellipsis))
  } else {
    do.call(plot, args =
              c(list(x = x, y = p[, 1], type = "l", xlim = c(min(x), max(x)),
                     ylim = c(0, 1), ylab = "Probability", xlab = "Theta",
                     main = "Category Response Curves"), ellipsis))
  }

  do.call(points, args =
            c(list(x = (loc + thres), y = rep(0, K),
                   pch = as.character(1:categ)), ellipsis))
  do.call(points, args =
            c(list(x = loc, y = 0, pch = 17), ellipsis))
  for (i in 2:categ) {
    do.call(lines, arg =
              c(list(x = x, y = p[, i]), ellipsis))
  }
}


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# SCATTERPLOT INCL. CORRELATION AS LEGEND -----------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Plot a default scatterplot and add the correlation coefficient as a legend
#'
#' detailed description
#'
#' @param ... See \link{code{plot.default}}.
#' @seealso \link{code{legend}}, \link{code{cor}}
#' @export
plot_cor <- function(...) {
#     dots <- function(...) {
#         eval(substitute(alist(...)))
#     }
    plot.default(...)
    c.1 <- round(do.call(cor, list(x = eval(substitute(alist(...)))[[1]],
                                   y = eval(substitute(alist(...)))[[2]])), 2)
    legend(x = ifelse(c.1 > 0, "topleft", "topright"),
           legend = paste("r =", c.1),
           box.col = "white", inset = .01)
}


