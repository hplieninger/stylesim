# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# REPLICATE (CRONBACH'S ALPHA) ----------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Replicate a Style Simulation and Investigate Effect on Cronbach's Alpha
#'
#' This function replicates \code{\link{sim_style_data}} and returns the observed
#' coefficient alpha as well as the response style-free alpha for every
#' replication sample.
#'
#' @param reps Numeric, the desired number of replications.
#' @param n Numeric, the number of persons. If of length one, it's fixed to the 
#'   provided value. If of length two, it's sampled from a uniform distribution 
#'   using the two values as lower and upper limits, respectively.
#' @param items Numeric, the number of items. If of length one, it's fixed to
#'   the provided value. If of length two, it's sampled from a uniform
#'   distribution using the two values as lower and upper limits, respectively.
#' @param categ Numeric, the number of categories per item. If of length one,
#'   it's fixed to the provided value. If of length two, it's sampled from a
#'   uniform distribution using the two values as lower and upper limits,
#'   respectively..
#' @param reversed Numeric, the number of reverse-coded items. If of length one,
#'   it's fixed to the provided value. If of length two, it's sampled from a
#'   uniform distribution using the two values as lower and upper limits,
#'   respectively.
#' @param var.s Numeric, the response style variance. If of length one, it's
#'   fixed to the provided value. If of length two, it's sampled from a uniform
#'   distribution using the two values as lower and upper limits, respectively.
#' @param mu.s Numeric, the response style mean. If of length one, it's fixed to
#'   the provided value. If of length two, it's sampled from a uniform
#'   distribution using the two values as lower and upper limits, respectively.
#' @param df Numeric. The df-parameter of the Wishart distribution from which
#'   the covariance is drawn.
#' @param ... Other parameters passed to \code{\link{sim_style_data}}.
#' @inheritParams sim_style_data
#' @return Returns a matrix of length \code{reps} with the following columns:
#' \item{bias}{\code{alpha} minus \code{true}}
#' \item{true}{Response style-free alpha}
#' \item{alpha}{Observed coefficient alpha}
#' \item{dep}{Response style-dependent alpha, equal to bias}
#' \item{ }{Further columns contain the input parameters such as the number of categories}
#' @seealso The replicated function \code{\link{sim_style_data}}, covariate-free alpha \code{\link{alpha_cov}}
#' @export
rep_style_rel <- function(reps = 1000, n = c(100, 1000), items = c(5, 10),
                          categ = c(3, 7), ndimc = 1, style = NULL,
                          reversed = c(0, .5), mu.s = c(-1, 1), var.s = c(0, 1),
                          df = 10, sig = NULL, emp = TRUE, ...) {
    time.1 <- Sys.time()
    cat(paste0(time.1, " - I'm gonna run that shit for you, dude.\n"))
    
    
    if (ndimc != 1) stop("Only defined for a single content-related dimension, please specify 'ndimc=1'.")
    if (!is.character(style)) {
        if (length(style) %% categ != 0) stop("Arguments 'style' is in disagreement with argument 'categ'")
        ndims <- length(style) / categ
    } else {
        ndims <- length(style)
    }
    if (ndims != 1) stop("Only defined for a single style-related dimension, please modify 'style'.")
    
    ndim <- ndimc + ndims
    
    res <- matrix(nrow = reps, ncol = 11 + ndimc*max(items) + max(categ) - 1)
    p.it <- vector("character")
    for (ii in 1:ndimc) {
        p.it <- c(p.it, paste0("b.item_d", ii, "_", 1:max(items)))
    }
    p.it <- c(p.it, paste0("b.categ", 1:(max(categ) - 1)))
    colnames(res)  <- c("bias", "true", "alpha", "n", "items", "categ", "rev",
                        "mu.s", "var.s", "cor_12", "dep", p.it)
    
    pb <- txtProgressBar(min = 0, max = reps, style = 3, char="-")
    
    for (ii in 1:reps) {
        if (length(n) == 1) {
            n.i <- n
        } else {
            n.i <- round(runif(1, n[1] - .49, n[2] + .49))
        }
        
        if (length(items) == 1) {
            items.i <- items
        } else {
            items.i <- sample(items[1]:items[2], 1)
        }
        
        if (length(categ) == 1) {
            categ.i <- categ
        } else {
            categ.i <- sample(categ[1]:categ[2], 1)
        }
        
        if (length(reversed) == 1) {
            rev.i <- ifelse(reversed >= 1, reversed, trunc(items.i * reversed))
        } else {
            rev.i <- sample(
                ifelse(reversed >= 1, reversed, trunc(items.i * reversed))[1]:
                    ifelse(reversed >= 1, reversed, trunc(items.i * reversed))[2], 1)
        }
        
        if (is.null(sig)) {
            if (length(var.s) == 1) {
                var.s.i <- rep(var.s, ndims)
            } else {
                var.s.i <- sample(seq(var.s[1], var.s[2], .1), ndims)
            }
            
            sig.i <- cov2cor(rWishart(1, df, diag(ndim))[, , 1])
            if (any(var.s.i > 0)) {
                sig.i <- MBESS::cor2cov(sig.i, c(rep(1, ndimc), sqrt(var.s.i)))
            } else {
                sig.i[(ndimc +1):ndim, ] <- 0
                sig.i[, (ndimc +1):ndim] <- 0
            }
        } else {
            sig.i <- sig
            var.s.i <- diag(sig.i)[(ndimc + 1):ndim]
        }
        
        if (length(mu.s) == 1) {
            mu.s.i <- rep(mu.s, ndims)
        } else {
            mu.s.i <- sample(seq(mu.s[1], mu.s[2], .1), ndims)
        }
        mu.vec.i <- c(rep(0, ndimc), mu.s.i)
        
        my.theta.i <- MASS::mvrnorm(n = n.i,
                                    mu = mu.vec.i,
                                    Sigma = sig.i,
                                    empirical = emp)

        d <- sim_style_data(n = n.i, items = items.i, categ = categ.i,
                            ndimc = ndimc, style = style, reversed = rev.i,
                            my.theta = my.theta.i, ...)

        res[ii, c("true", "dep", "alpha")] <-
            alpha_cov(cbind(d$dat[, , 1], d$theta[, (ndimc + 1)]),
                      xx = 1:items.i,
                      zz = ifelse(var.s.i > 0, (items.i + 1):(items.i + ndims), rep(NA, ndims)))
        res[ii, "n"] <- n.i
        res[ii, "items"] <- items.i
        res[ii, "categ"] <- categ.i
        res[ii, "rev"] <- rev.i
        res[ii, "mu.s"] <- mu.s.i
        res[ii, "var.s"] <- sig.i[ndim, ndim]
        suppressWarnings({
            cors.i <- cov2cor(sig.i)[upper.tri(diag(ndim))]
        })
        cors.i[is.nan(cors.i)] <- 0
        res[ii, "cor_12"] <- cors.i
        for (jj in 1:ndimc) {
            res[ii, paste0("b.item_d", jj, "_", 1:items.i)] <- 
                d$item.parameters[(1+items.i*(jj-1)):(items.i*jj)]
        }
        res[ii, c(paste0("b.categ", 1:(categ.i - 1)))] <- d$item.parameters[grep("categ", names(d$item.parameters))]

        setTxtProgressBar(pb, ii)
    }
    close(pb)
    
    res[, "bias"] <- res[, "alpha"] - res[, "true"]
    time.2 <- difftime(Sys.time(), time.1, units = "hours")
    styleweights <- sim_style_data(n = n.i, items = items.i, categ = max(categ),
                                  ndimc = ndimc, style = style, reversed = rev.i,
                                  my.theta = my.theta.i, ...)$response.style
    res.list <- list(data = res, style = styleweights[1:2], time.elapsed = time.2,
                     dv = "rel", setup = list(df_of_Wishart = df, dims.content = ndimc))
    cat("Total time elapsed:", round(time.2[[1]], 3), "hours\n")
    return(res.list)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# REPLICATE (CONTENT-CONTENT CORRELATION) -----------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Replicate a Style Simulation and Investigate Effect on Correlation of Two Scale Scores
#'
#' This function replicates \code{\link{sim_style_data}} and calculates a scale score
#' (i.e., the mean across all items) for every person for every dimension. It
#' mainly returns two values for every replication sample, namely, the observed
#' correlation of the two scales as well as their partial correlation
#' controlling for response style.
#'
#' @inheritParams sim_style_data
#' @inheritParams rep_style_rel
#' @param ... Other parameters passed to \code{\link{sim_style_data}}.
#' @return Returns a matrix of length \code{reps} with the following columns:
#' \item{bias}{\code{obs} minus \code{partial}}
#' \item{partial}{The partial correlation \eqn{r_{\bar{X}_1, \bar{X}_2 \cdot \theta_{\textrm{RS}}}}}
#' \item{obs}{The observed correlation}
#' \item{ }{Further columns contain the input parameters such as the number of categories}
#' @seealso \code{\link{sim_style_data}}
#' @export
rep_style_cor <- function(reps = 1000, n = c(100, 1000), items = c(5, 10),
                          categ = c(3, 7), ndimc = 2, style = NULL,
                          reversed = c(0, .5), mu.s = c(-1, 1), var.s = c(0, 1), 
                          df = 10, sig = NULL, emp = TRUE, ...) {
    time.1 <- Sys.time()
    cat(paste0(time.1, " - I'm gonna run that shit for you, dude.\n"))
    
    # ndimc > 2 would require the partial correlation to be accurately defined.
    if (ndimc != 2) {
        stop("Currently only defined for TWO content-related dimensions, please specify 'ndimc=2'.")
    }
    if (is.null(style)) {
        stop("'style = NULL' is not implemented, but 'var.s = 0' will do the job.")
    }
    if (!is.character(style)) {
        if (length(style) %% categ != 0) stop("Arguments 'style' is in disagreement with argument 'categ'")
        ndims <- length(style) / categ
    } else {
        ndims <- length(style)
    }
    if (ndims != 1) stop("Only defined for a single style-related dimension, please modify 'style'.")
    
    ndims <- length(style)
    ndim <- ndimc + ndims
    sty <- ifelse(is.null(style), FALSE, TRUE)
    
    res <- matrix(nrow = reps, ncol = 11 + ndim * (ndim - 1) / 2 + ndimc*max(items) + max(categ) - 1)
    cors <- which(upper.tri(diag(ndim)) == TRUE, arr.ind = T)
    cors <- paste0("cor_", cors[, 1], cors[, 2])
    p.it <- vector("character")
    for (ii in 1:ndimc) {
        p.it <- c(p.it, paste0("b.item_d", ii, "_", 1:max(items)))
    }
    p.it <- c(p.it, paste0("b.categ", 1:(max(categ) - 1)))
    colnames(res)  <- c("bias", "partial", "obs", "r.TT", "bias2", "n",
                        "items", "categ", "rev", "var.s", "mu.s", cors, p.it)
    
    pb <- txtProgressBar(min = 0, max = reps, style = 3, char="-")
    
    for (ii in 1:reps) {
        if (length(n) == 1) {
            n.i <- n
        } else {
            n.i <- round(runif(1, n[1] - .49, n[2] + .49))
        }
        
        if (length(items) == 1) {
            items.i <- items
        } else {
            items.i <- sample(items[1]:items[2], 1)
        }
        
        if (length(categ) == 1) {
            categ.i <- categ
        } else {
            categ.i <- sample(categ[1]:categ[2], 1)
        }
        
        if (length(reversed) == 1) {
            rev.i <- ifelse(reversed >= 1, reversed, trunc(items.i * reversed))
        } else {
            rev.i <- sample(
                ifelse(reversed >= 1, reversed, trunc(items.i * reversed))[1]:
                    ifelse(reversed >= 1, reversed, trunc(items.i * reversed))[2], 1)
        }
        
        if (is.null(sig)) {
            if (length(var.s) == 1) {
                var.s.i <- rep(var.s, ndims)
            } else {
                var.s.i <- sample(seq(var.s[1], var.s[2], .1), ndims)
            }
            sig.i <- cov2cor(rWishart(1, df, diag(ndim))[, , 1])
            if (any(var.s.i > 0)) {
                sig.i <- MBESS::cor2cov(sig.i, c(rep(1, ndimc), sqrt(var.s.i)))
            } else {
                sig.i[(ndimc +1):ndim, ] <- 0
                sig.i[, (ndimc +1):ndim] <- 0
            }
        } else {
            sig.i <- sig
            var.s.i <- diag(sig.i)[(ndimc + 1):ndim]
        }
        if (length(mu.s) == 1) {
            mu.s.i <- rep(mu.s, ndims)
        } else {
            mu.s.i <- sample(seq(mu.s[1], mu.s[2], .1), ndims)
        }
        mu.vec.i <- c(rep(0, ndimc), mu.s.i)
        
        my.theta.i <- MASS::mvrnorm(n = n.i,
                                    mu = mu.vec.i,
                                    Sigma = sig.i,
                                    empirical = emp)

        d <- sim_style_data(n = n.i, items = items.i, categ = categ.i,
                            ndimc = ndimc, style = style, reversed = rev.i,
                            my.theta = my.theta.i, ...)
        
        if (all(sty, var.s.i > 0)) {
            res[ii, "partial"] <- partial_cor(cbind(rowMeans(d$dat[, , 1, drop=F]),
                                            rowMeans(d$dat[, , 2, drop=F]), d$theta[, 3]))
        } else {
            res[ii, "partial"] <- cor(rowMeans(d$dat[, , 1, drop=F]), rowMeans(d$dat[, , 2, drop=F]))
        }
        res[ii, "obs"] <- cor(rowMeans(d$dat[, , 1, drop=F]), rowMeans(d$dat[, , 2, drop=F]))
        res[ii, "r.TT"] <- cor(d$theta[, 1], d$theta[, 2])
        res[ii, "n"] <- n.i
        res[ii, "items"] <- items.i
        res[ii, "categ"] <- categ.i
        res[ii, "rev"] <- rev.i
        res[ii, "var.s"] <- ifelse(sty, var.s.i, NA)
        res[ii, "mu.s"] <- mu.s.i
        suppressWarnings({
            cors.i <- cov2cor(sig.i)[upper.tri(diag(ndim))]
        })
        cors.i[is.nan(cors.i)] <- 0
        res[ii, colnames(res) %in% cors] <- cors.i
        for (jj in 1:ndimc) {
            res[ii, paste0("b.item_d", jj, "_", 1:items.i)] <- 
                d$item.parameters[(1+items.i*(jj-1)):(items.i*jj)]
        }
        res[ii, c(paste0("b.categ", 1:(categ.i - 1)))] <- d$item.parameters[grep("categ", names(d$item.parameters))]
        
        setTxtProgressBar(pb, ii)
    }
    close(pb)
    
    res[, "bias"] <- res[, "obs"] - res[, "partial"]
    res[, "bias2"] <- res[, "obs"] - res[, "r.TT"]
    time.2 <- difftime(Sys.time(), time.1, units = "hours")
    styleweights <- sim_style_data(n = n.i, items = items.i, categ = max(categ),
                                   ndimc = ndimc, style = style, reversed = rev.i,
                                   my.theta = my.theta.i, ...)$response.style
    res.list <- list(data = res, style = styleweights[1:2], time.elapsed = time.2,
                     dv = "cor", setup = list(df_of_Wishart = df, dims.content = ndimc))
    cat("Total time elapsed:", round(time.2[[1]], 3), "hours\n")
    return(res.list)
}