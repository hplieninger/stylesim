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
#' @seealso Alternatives are, for example, \code{\link[psych]{partial.r}}.
partial_cor <- function(x) {
  # x is a matrix of 3 columns containing the variables X, Y, and Z
  x <- cor(x)
  (x[1,2] - x[1,3] * x[2,3]) / (sqrt(1 - x[1,3]^2) * sqrt(1 - x[2,3]^2))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# INTERACTION PLOT ----------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Plot a default scatterplot and add simple slope regression lines in order to
#' visualize interaction effects
#' 
#' This function may help to understand and interpret interaction effects
#' observed in a linear model (\code{\link[stats]{lm}}) by means of simple
#' slopes.
#' 
#' @param data A data frame, the same one that was used to fit the model.
#' @param model An object of class (\code{\link[stats]{lm}}). Make sure to use 
#'   the \code{data} argument of \code{\link[stats]{lm}} (i.e., do not use 
#'   something like \code{lm(d$y ~ d$x1*d$x2)}).
#' @param y Character. The dependent variable in the \code{model}.
#' @param x Character vector. The independent variable(s) in the \code{model}
#'   for which the interaction effects should be plotted. If this is of length >
#'   1, the first element will be used as the x-variable of the plot and the
#'   second element will be used for different simple slopes (and the third
#'   variable will be used for different plot panels).
#' @param xq Character. Optional term in your model, namely, the quadratic 
#'   effect of your first independent variable.
#' @param lvl List of numeric vectors. The list must be of length 1 with two 
#'   predictors and of length 2 with three predictors. Each vector contains the 
#'   values of the corresponding predictor at which the plot is to be drawn.
#'   Defaults to three values, namely \eqn{M - 1SD}, \eqn{M}, and \eqn{M + 1SD}.
#' @param N Numeric. A number of data points to plot; possibly smaller than 
#'   \code{nrow(data)}.
#' @param lwd Numeric. The line width(s) for the line(s) (see
#'   \code{\link[graphics]{par}}).
#' @param lty Numeric. The line type(s) for the line(s) (see
#'   \code{\link[graphics]{par}}).
#' @param col Character. The color of the points.
#' @param alpha Numeric. An alpha transparency value (as an opacity, so 0 means 
#'   fully transparent and 1 means opaque; see \code{\link[grDevices]{rgb}}).
#' @param col.s Character vector. The color(s) of the slope(s). Must be of 
#'   length 1 for one predictor and equal to \code{length(lvl[[1]])} for 2 
#'   predictors.
#' @param cex.legend Numeric. A value giving the amount by which the legend 
#'   should be magnified (see \code{\link[graphics]{par}}).
#' @param plot.legend Logical. Should a legend be plotted.
#' @param at.x,at.y Numeric. The points at which tick-marks for the x- and the 
#'   y-axis, respectively, are to be drawn (see \code{\link[graphics]{axis}}).
#' @param labels.x,labels.y This can either be a logical value specifying
#'   whether (numerical) annotations are to be made at the tickmarks of the x-
#'   and y-axsi, respetively, or a character or expression vector of labels to
#'   be placed at the tickpoints. If this is not logical, \code{at.x} and/or 
#'   \code{at.y} should also be supplied and of the same length (see 
#'   \code{\link[graphics]{axis}}).
#' @param cex.x.axis,cex.y.axis The magnification to be used for axis annotation
#'   relative to the current setting of cex (see \code{\link[graphics]{par}}).
#' @param grid Logical. Should a grid be plotted.
#' @param grid.h Numeric. Values of the y-axis at which the horizontal grid 
#'   lines should be plotted.
#' @param padj.x,padj.y Numeric. Adjustment for each tick label perpendicular to
#'   the reading direction. For labels parallel to the axes, padj = 0 means 
#'   right or top alignment, and padj = 1 means left or bottom alignment (see 
#'   \code{\link[graphics]{axis}}).
#' @param xlim,ylim Numeric vectors of length 2, giving the x and y coordinates ranges.
#' @inheritParams graphics::par
#' @inheritParams graphics::legend
#' @inheritParams graphics::title
# @inheritParams graphics::plot.window
#' @seealso \code{\link[graphics]{par}}, \code{\link[graphics]{plot}}, and
#'   \code{\link{legend}}
plot_x <- function(data, model, y, x, xq = NULL,
                   lvl = NULL, N = NULL,
                   ylab = y, xlab = x[1], lwd = 2, lty = 1,
                   col = "black", pch = 20, alpha = .05, col.s = NULL,
                   xlim = range(data[, x[1]]), ylim = range(data[, y]),
                   seg.len = .5, cex.legend = .8, plot.legend = TRUE,
                   at.x = NULL, at.y = NULL, labels.x = TRUE, labels.y = TRUE,
                   cex.x.axis = 1, cex.y.axis = 1, grid = TRUE, grid.h = NULL,
                   padj.x = NA, padj.y = NA) {
    xlim; ylim
    opar <- par(no.readonly = T)
    if (is.null(N)) N <- min(c(nrow(data), 2000))
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
    
    b[1:length(b)] <- model$coefficients[x3][1:length(b)]
    if (any(is.na(b))) stop("Something went wrong while collecting the regression weights. Possibly, a desired interaction was not included in the 'model' you provided.")
    
    if (length(x) == 1) {
        data <- data[sample(nrow(data), N), ]
        plot.new()
            plot.window(xlim = xlim, ylim = ylim)
            axis(1, at = at.x, labels = labels.x, cex.axis = cex.x.axis,
                 padj = padj.x)
            axis(2, at = at.y, labels = labels.y, cex.axis = cex.y.axis,
                 padj = padj.y)
            title(ylab = ylab, xlab = xlab)
            box()
            if (grid == TRUE) {
                if (is.null(grid.h)) grid.h <- axTicks(2)
                abline(h = grid.h, lty = 3, col = "grey75")
            }
        points(data[, x[1]], data[, y], pch = pch,
               col = rgb(t(col2rgb(col))/255, alpha = alpha))
        
        curve(b["b0"]
              + b["b1"]*x
              + ifelse(is.null(xq), 0, xq)*x^2
              , add = T, col = col.s[ii], lwd = lwd, lty = lty,
              from = xlim[1]*1.2, to = xlim[2]*1.2)
    }
        
    if (length(x) == 2) {
        data <- data[sample(nrow(data), N), ]
        plot.new()
            plot.window(xlim = xlim, ylim = ylim)
            axis(1, at = at.x, labels = labels.x, cex.axis = cex.x.axis,
                 padj = padj.x)
            axis(2, at = at.y, labels = labels.y, cex.axis = cex.y.axis,
                 padj = padj.y)
            title(ylab = ylab, xlab = xlab)
            box()
            if (grid == TRUE) {
                if (is.null(grid.h)) grid.h <- axTicks(2)
                abline(h = grid.h, lty = 3, col = "grey75")
            }
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
            axis(1, at = at.x, labels = labels.x, cex.axis = cex.x.axis,
                 padj = padj.x)
            axis(2, at = at.y, labels = labels.y, cex.axis = cex.y.axis,
                 padj = padj.y)
            title(ylab = ylab, xlab = xlab,
                  main = paste0(x[3], " = ", round(lvl[[2]][kk], 2)))
            box()
            if (grid == TRUE) {
                if (is.null(grid.h)) grid.h <- axTicks(2)
                abline(h = grid.h, lty = 3, col = "grey75")
            }
            points(data[id.idx, x[1]], data[id.idx, y], pch = pch,
                   col = rgb(t(col2rgb(col))/255, alpha = alpha))
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
#  COVARIATE-FREE AND COVARIATE-DEPENDENT CRONBACH'S ALPHA ------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Calculate covariate-free and covariate-dependent Cronbach's alpha
#'
#' If a covariate (or possibly a set of covariates) is partialled out from a set
#' of target variables, their corresponding coefficient alpha can be decomposed 
#' into a part that is dependent on that covariate and a part that is 
#' independent ("free") from that covariate (and the sum equals ordinary
#' coefficient alpha).
#' 
#' @references Bentler, P. (2015). \emph{Covariate-free and covariate-dependent
#'   reliability}. Manuscript submitted for publication. Retrieved from
#'   \url{http://escholarship.org/uc/item/6v61v1pk}
#'
#' @param data A data frame or matrix containing both the variables as well as
#'   the covariate(s). This can also be a covariance matrix.
#' @param xx Numeric. A vector of indices to specify the columns in \code{data}
#'   that pertain to the target variables.
#' @param zz Numeric. A vector of indices to specify the columns in \code{data}
#'   that pertain to the covariate(s).
#' @inheritParams stats::cov
#' @return Returns a named vector with covariate-free, covariate-dependent, and
#'   total Cronbach's alpha (unstandardized).
#' @seealso Cronbach's alpha in the \strong{psych} package: \code{\link[psych]{alpha}}
#' @export
alpha_cov <- function(data, xx, zz = NULL, use = "pairwise") {
    a1 <- rep(NA, 3)
    p <- ncol(data[, xx])
    if (nrow(data) == ncol(data)) {
        EE <- data
    } else {
        EE <- cov(data, use = use)
    }
    den <- sum(EE[xx, xx])
    if (!is.null(zz)) {
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
#' Make a barplot and include the call to \link{barplot} in the main title
#'
#' Internal helper function
#'
#' @param ... Parameters passed to \code{\link{barplot}}
barplot_main <- function(...) {
    barplot(main = paste(match.call()[2]), ...)
}
