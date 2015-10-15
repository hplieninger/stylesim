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
#' Compute the effect size \eqn{f^2} for nested (hierarchical) regression models
#'
#' detailed description
#'
#' @param model a model of class \code{lm}
#' @return A vector containing the $f^2$ values for every coefficient in the model
#' @seealso \code{\link[lmSupport]{modelEffectSizes}}.
#' @export
f_squared <- function(model, digits = 3) {
    r2 <- matrix(NA, nrow = model$rank - 1, ncol = 2,
                 dimnames = list(names(model$coefficients)[-1], c("R2_1", "R2_2")))
    x1 <- attr(model$terms, "factors")
    x2 <- attr(model$terms, "order")
    x3 <- attr(model$terms, "term.labels")
    if (any(grepl("^", x3, fixed = T))) {
        x4 <- sort(x3[grepl("^", x3, fixed = T)])
        x5 <- substring(x4, first = unlist(gregexpr("(", x4, fixed = T)) + 1,
                        last = unlist(gregexpr("^", x4, fixed = T)) - 1)
        x6 <- as.integer(substring(x4, first = unlist(gregexpr("^", x4, fixed = T)) + 1,
                                   last = unlist(gregexpr("^", x4, fixed = T)) + 1))
        # the code for x5 and x6 fails for higher*higher-interactions (e.g.,
        # 'I(a^2)*I(b^2)'). Therefore, I experimented with x5a and x6a below, but
        # that's not enough, situation gets realy complicated.
#         x5a <-  vector("list", length(x4)) 
#         for (ii in seq_along(x4)) {
#             x5a[[ii]] <- substring(x4[ii], first = unlist(gregexpr("(", x4, fixed = T)[ii]) + 1,
#                       last = unlist(gregexpr("^", x4, fixed = T)[ii]) - 1)
#         }
#         x6a <-  vector("list", length(x4)) 
#         for (ii in seq_along(x4)) {
#             x6a[[ii]] <- as.integer(substring(x4[ii], first = unlist(gregexpr("^", x4, fixed = T)[ii]) + 1,
#                                               last = unlist(gregexpr("^", x4, fixed = T)[ii]) + 1))
#         }

        a1 <- matrix(NA, nrow = (max(x6)-1), ncol = length(x6),
                     dimnames = list(NULL, x4))
        a1[1, ] <- x5
        
        if (nrow(a1) > 1) {
            for (kk in 2:(max(x6)-1)) {
                for (ii in 1:length(x6)) {
                    if (x6[ii] > kk) a1[kk, ii] <- x4[ii -( kk - 2 + 1)]
                }
            }
        }

        for (ii in 1:ncol(a1)) {
            x1[na.omit(a1[, ii]), colnames(a1)[ii]] <- 1
        }

    }
    drop.x <- vector("list", length = length(x3))
    for (ii in 1:length(drop.x)) {
        rows.x <- x1[, x3[ii]] == 1
        cols.x <- (apply(x1[rows.x, , drop = F], 2, function(x) all(x > 0))) &
                      (x2 >= x2[ii])
        drop.x[[ii]] <- colnames(x1[rows.x, cols.x, drop = F])
    }
    # drop.x <- c(1, drop.x)
    for (ii in 1:length(drop.x)) {

        r2[ii, 1] <- summary(do.call("update",
                                     args=list(model,
                                               paste("~ .", paste("-", drop.x[[ii]], collapse = " "))))
        )$r.squared
        if (length(drop.x[ii]) > 1) {
            r2[ii, 2] <- summary(do.call("update",
                                         args=list(model,
                                                   paste("~ .", paste("-", drop.x[-1], collapse = " "))))
            )$r.squared
        } else {
            r2[ii, 2] <- summary(model)$r.squared
        }
    }
    f2 <- (r2[, 2, drop = F] - r2[, 1, drop = F]) / (1 - r2[, 2, drop = F])
    colnames(f2) <- "f2"
    print(round(f2[f2[, 1] > 10^(-digits), , drop = F], digits = digits))
    invisible(f2)
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
# @seealso
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
    do.call(lines, args =
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
#' @param ... See \code{\link{plot.default}}.
#' @seealso \code{\link{legend}}, \code{\link{cor}}
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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# INTERACTION PLOT ----------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Plot a default scatterplot and add simple slope regression lines to visualize interactions
#'
#' detailed description
#'
#' @param data A data frame.
#' @param model An object of class 'lm'.
#' @param y The dependent variable.
#' @param x The independent variables.
#' @param lvl The levels of x[2] for the simple slopes.
#' @param ... See \code{\link{plot.default}}.
#' @seealso \code{\link{legend}}, \code{\link{cor}}
#' @export
plot_x <- function(data, model, y, x, xq = NULL,
                   lvl = NULL, N = 1000,
                   ylab = y, xlab = x[1], lwd = 2, lty = 1,
                   col = "black", pch = 20, alpha = .2, col.s = NULL,
                   xlim = range(data[, x[1]]), ylim = range(data[, y]),
                   seg.len = .5, cex.legend = .8, plot.legend = TRUE,
                   at.x = NULL, at.y = NULL, labels.x = TRUE, labels.y = TRUE,
                   grid = TRUE, grid.h = NULL) {
    xlim; ylim
    opar <- par(no.readonly = T)
    xes <- length(x)
    if (is.null(lvl) & xes > 1) {
        lvl <- as.list(as.data.frame(apply(data[, x[-1], drop = F], 2, function(x) {
            c(mean(x, na.rm = T) - sd(x, na.rm = T), mean(x, na.rm = T), mean(x, na.rm = T) + sd(x, na.rm = T))
        })))
    }
    if (xes == 1) n.lvl <- list(1)
    else n.lvl <- lapply(lvl, length)
    if (is.null(col.s)) col.s <- rainbow(n.lvl[[1]])
    if (length(lty) == 1) lty <- rep(lty, n.lvl[[1]])
    
    x1 <- attr(model$terms, "factors")[x, , drop = F]
    x2 <- attr(model$terms, "order")
    
    x3 <- c("(Intercept)")
    for (kk in 1:xes) {
        for (ii in 1:sum(rep(1:xes, times = choose(xes, 1:xes)) == kk)) {
            x3 <- c(x3,
                    names(which(
                        apply((x1[combn(xes, kk, drop = F)[, ii], x2 == kk, drop = F] == 1), 2, all)
                    )))
        }
    }
    
    # int <- names(which(apply((x1 == 1), 2, function(x) sum(x) > 1)))
    xes <- length(x)
    a1 <- expand.grid(rep(list(1:xes), xes))
    try(a1 <- a1[!apply(t(diff(t(as.matrix(a1)))) < 0, 1, any), ], TRUE)
    a1 <- apply(a1, 1, function(x) unique(x))
    a1 <- unique(unlist(lapply(a1, paste, collapse = "")))
    a1 <- a1[order(a1)]
    a1 <- a1[order(nchar(a1))]
    
    b <- rep(NA, length = 1 + length(a1))
    names(b) <- c("b0", paste0("b", a1))
    
#     b.x <- names(b)
#     b.x <- gsub("b", "", b.x)
#     for (ii in 1:length(x)) {
#         b.x <- gsub(ii, rainbow(length(x))[ii], b.x)
#     }
#     for (ii in 1:length(x)) {
#         b.x <- gsub(rainbow(length(x))[ii], paste0(x[ii], ":"), b.x)
#     }
#     reverse_chars <- function(string) {
#         string_split = strsplit(string, split = "")
#         rev_order = nchar(string):1
#         reversed_chars = string_split[[1]][rev_order]
#         paste(reversed_chars, collapse="")
#     }
#     b.x <- sapply(b.x[-1], reverse_chars)
#     b.x <- sub(":", "", b.x)
#     b.x <- c("(Intercept)", sapply(b.x, reverse_chars))
    
    b[1:length(b)] <- model$coefficients[x3][1:length(b)]
    if (any(is.na(b))) stop("Something went wrong while collecting the regression weights. Possibly, a desired interaction was not included in the 'model' you provided.")
    
    if (length(x) == 1) {
        data <- data[sample(nrow(data), N), ]
        plot.new()
            plot.window(xlim = xlim, ylim = ylim)
            axis(1, at = at.x, labels = labels.x)
            axis(2, at = at.y, labels = labels.y)
            title(ylab = ylab, xlab = xlab)
            box()
            if (grid == TRUE) {
                if (is.null(grid.h)) grid.h <- axTicks(2)
                abline(h = grid.h, lty = 3, col = "grey75")
            }
        # points(jitter(data[, x[1]]), jitter(data[, y]), pch = pch)
        points(data[, x[1]], data[, y], pch = pch,
               col = rgb(t(col2rgb(col))/255, alpha = alpha))
        
        curve(b["b0"]
              + b["b1"]*x
              + ifelse(is.null(xq), 0, xq)*x^2
              , add = T, col = col.s[ii], lwd = lwd, lty = lty,
              from = xlim[1]*1.2, to = xlim[2]*1.2)
#         legend(x = "topleft", lty = lty, bty = "n", title = x[2],
#                legend = round(lvl[[1]], 2),
#                col = col.s, seg.len = seg.len, lwd = lwd, cex = cex.legend)
    }
        
    if (length(x) == 2) {
        data <- data[sample(nrow(data), N), ]
        plot.new()
            plot.window(xlim = xlim, ylim = ylim)
            axis(1, at = at.x, labels = labels.x)
            axis(2, at = at.y, labels = labels.y)
            title(ylab = ylab, xlab = xlab)
            box()
            if (grid == TRUE) {
                if (is.null(grid.h)) grid.h <- axTicks(2)
                abline(h = grid.h, lty = 3, col = "grey75")
            }
        # points(jitter(data[, x[1]]), jitter(data[, y]), pch = pch)
        points(data[, x[1]], data[, y], pch = pch,
               col = rgb(t(col2rgb(col))/255, alpha = alpha))
        
        for (ii in 1:n.lvl[[1]]) {
            curve(b["b0"]
                  + b["b1"]*x
                  + b["b2"]*lvl[[1]][ii]
                  + b["b12"]*lvl[[1]][ii]*x
                  + ifelse(is.null(xq), 0, xq)*x^2
                  ,add = T, col = col.s[ii], lwd = lwd, lty = lty[ii],
                  from = xlim[1]*1.2, to = xlim[2]*1.2)
        }
        legend(x = "topleft", lty = lty, bty = "n", title = x[2],
               legend = round(lvl[[1]], 2),
               col = col.s, seg.len = seg.len, lwd = lwd, cex = cex.legend,
               plot = plot.legend)
    }
    if (length(x) == 3) {
        data <- data[sample(nrow(data), ifelse(N*n.lvl[[2]] > N, N, N*n.lvl[[2]])), ]
        par(mfrow = c(1, n.lvl[[2]]))
        for (kk in 1:n.lvl[[2]]) {
            id.cut <- c(-Inf, (lvl[[2]][2:3] + lvl[[2]][1:2]) / 2, Inf)
            id.idx <- (data[, x[3]] >= id.cut[kk]) & (data[, x[3]] < id.cut[kk+1])
            plot.new()
            plot.window(xlim = xlim, ylim = ylim)
            axis(1, at = at.x, labels = labels.x)
            axis(2, at = at.y, labels = labels.y)
            title(ylab = ylab, xlab = xlab,
                  main = paste0(x[3], " = ", round(lvl[[2]][kk], 2)))
            box()
            if (grid == TRUE) {
                if (is.null(grid.h)) grid.h <- axTicks(2)
                abline(h = grid.h, lty = 3, col = "grey75")
            }
            points(data[id.idx, x[1]], data[id.idx, y], pch = pch,
                   col = rgb(t(col2rgb(col))/255, alpha = alpha))
#             plot(jitter(data[id.idx, x[1]]), jitter(data[id.idx, y]),
#                  ylim = ylim,
#                  pch = pch, ylab = ylab, xlab = xlab,
#                  main = paste0(x[3], " = ", round(mean(c(lvl[[2]][kk], lvl[[2]][kk+1])), 2)))
            for (ii in 1:length(lvl[[1]])) {
                curve(b["b0"]
                      + b["b1"]*x
                      + b["b2"]*lvl[[1]][ii]
                      + b["b3"]*lvl[[2]][kk]
                      + b["b12"]*x*lvl[[1]][ii]
                      + b["b13"]*x*lvl[[2]][kk]
                      + b["b23"]*lvl[[1]][ii]*lvl[[2]][kk]
                      + b["b123"]*x*lvl[[1]][ii]*lvl[[2]][kk]
                      + ifelse(is.null(xq), 0, xq)*x^2
                      , add = T, col = col.s[ii], lwd = lwd, lty = lty,
                      from = xlim[1]*1.2, to = xlim[2]*1.2)
            }
            if (kk == 2) {
                legend(x = "topleft", lty = lty, bty = "n", title = x[2],
                       legend = round(lvl[[1]], 2),
                       col = col.s, seg.len = seg.len, lwd = lwd, cex = cex.legend,
                       plot = plot.legend)
            }
        }
    }
    # on.exit(par(opar))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#  COVARIATE-FREE AND COVARIATE-BASED CRONBACH'S ALPHA ----------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Calculate covariate-free and covariate-dependent Cronbach's Alpha
#'
#' detailed description
#'
#' @param data A data frame or matrix containing both the variables as well as the covariate(s).
#' @param xx A vector of indices giving the variables' columns in \code{data}
#' @param zz A vector of indices giving the covariates' columns in \code{data}. May have length of one.
#' @return Returns a named vector with covariate-free, covariate-dependent, and total Cronbach's Alpha.
#' @export
alpha_cov <- function(data, xx, zz = NA, use = "pairwise") {
    a1 <- rep(NA, 3)
    p <- ncol(data[, xx])
    if (nrow(data) == ncol(data)) {
        EE <- data
    } else {
        EE <- cov(data, use = use)
    }
    den <- sum(EE[xx, xx])
    if (!any(is.na(zz))) {
        EEf <- (EE[xx, xx, drop = F] - EE[xx, zz, drop = F] %*% solve(EE[zz, zz, drop = F]) %*% EE[zz, xx, drop = F])
        EEd <- EE[xx, zz, drop = F] %*% solve(EE[zz, zz, drop = F]) %*% EE[zz, xx, drop = F]
        # EE[xx, xx] == EEf + EEd
        diag(EEf) <- NA
        diag(EEd) <- NA
        a1[2] <- p^2 * (mean(EEd, na.rm = T)) / den
    } else {
        EEf <- EE[xx, xx, drop = F]
        diag(EEf) <- NA
    }
    
    a1[1] <- p^2 * (mean(EEf, na.rm = T)) / den
    a1[3] <- sum(a1[1:2], na.rm = T)
    
    names(a1) <- c("covariate-free", "covariate-dependent", "total")
    return(a1)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#  BARPLOT WITH AUTOMATIC MAIN TITLE ----------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Make a barplot and include the call to \code{\link{barplot}} as the main title
#'
#' detailed description
#'
#' @param ... Parameters passed to \code{\link{barplot}}
#' @export
barplot_main <- function(...) {
    barplot(main = paste(match.call()[2]), ...)
}
