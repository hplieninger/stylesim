# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# REPLICATE (CRONBACH'S ALPHA) ----------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Replicate a Style Simulation and Investigate Effect on Cronbach's Alpha
#'
#' This function replicates \code{sim_style_data()} and returns the true
#' reliability as well as observed Cronbach's Alpha for every replication sample.
#'
#' @param reps Desired number of replications
#' @param n Optional. Either fixed to the provided value or randomly sampled for
#'   each replication.
#' @param items Optional. Either fixed to the provided value or randomly sampled
#'   for each replication.
#' @param categ Optional. Either fixed to the provided value or randomly sampled
#'   for each replication.
#' @param reversed Optional. Either fixed to the provided value or randomly
#'   sampled for each replication.
#' @param var.s Optional. Either fixed to the provided value or randomly sampled
#'   for each replication.
#' @inheritParams sim_style_data
#' @return Returns a matrix of length \code{reps} containing the true
#'   reliability in the first column and Cronbach's Alpha in the second column.
#'   The  true reliability is calculated via \eqn{r^2_{\theta_1, \theta_2}}.
#'   Further columns contain the input parameters such as the number of
#'   categories.
#' @seealso The replicated function \code{\link{sim_style_data}}. Cronbach's
#'   alpha: \code{\link[ltm]{cronbach.alpha}}
#' @export
#' @importFrom psych alpha
rep_style_rel <- function(reps = 1000, n = c(100, 1000), items = c(5, 10),
                          categ = c(3, 7), ndimc = 1, style = NULL,
                          reversed = c(0, .5), var.s = c(0, 1), mu.s = c(-1, 1),
                          ...) {
    
    if (ndimc != 1) stop("Only defined for a single content-related dimension,
                       please specify 'ndimc=1'.")
    if (length(style) != 1) stop("Only defined for a single style-related dimension, please modify 'style'.")
    # ndims <- ifelse(is.character(style), length(style), length(style)/categ)
    
    res <- matrix(nrow=reps, ncol=9)
    colnames(res)  <- c("bias", "true", "alpha", "n", "items", "categ", "rev",
                        "var.s", "mu.s")
    
    pb <- txtProgressBar(min = 0, max = reps, style = 3, char="-")
    
    for (i in 1:reps) {
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
        
        if (length(var.s) == 1) {
            var.s.i <- var.s
        } else {
            var.s.i <- sample(seq(var.s[1], var.s[2], .1), 1)
        }
        
        sig.x <- diag(1 + ndims)
        sig.x[2, 2] <- var.s.i
        sig.i <- cov2cor(rWishart(1, 10, sig.x)[, , 1])
        diag(sig.i) <- diag(sig.x)
        
        
        mu.vec <- c(rep(0, ndimc), mu.s)
        my.theta.i <- t(MASS::mvrnorm(n = n.i, mu = mu.vec, Sigma = sig, empirical = emp))
        
        
        if (length(mu.s) == 1) {
            mu.s.i <- mu.s
        } else {
            mu.s.i <- sample(seq(mu.s[1], mu.s[2], .1), 1)
        }
        
        #     if (args.prov[2]) {
        #       items <- sample(6:12, 1)
        #     }
        #     if (args.prov[3]) {
        #       categ <- sample(3:8, 1)
        #     }
        #     if (args.prov[4]) {
        #       reversed <- sample(0:trunc(items / 2), 1)
        #     }
        #     if (args.prov[5]) {
        #       var.s <- sample(seq(0, 1, .1), 1)
        #     }
        #     if (args.prov[6]) {
        #       mu.s <- sample(seq(0, 1, .1), 1)
        #     }
        
        d <- sim_style_data(n=n.i, items=items.i, categ=categ.i, ndimc=ndimc, style=style,
                            reversed=rev.i, var.s=var.s.i, mu.s=mu.s.i, ...)
        
        res[i, 2] <- cor(rowMeans(d$dat[, , 1, drop=F]), d$theta[, 1])^2
        # res[i, 3] <- cronbach.alpha(d$dat[, , 1])$alpha
        res[i, 3] <- alpha(as.data.frame(d$dat[, , 1]))$total[, 1]
        res[i, c(4:9)] <- c(n.i, items.i, categ.i, rev.i, var.s.i, mu.s.i)
        
        setTxtProgressBar(pb, i)
    }
    
    res[, 1] <- res[, 3] - res[, 2]
    return(res)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# REPLICATE (CONTENT-CONTENT CORRELATION) -----------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Replicate a Style Simulation and Investigate Effect on Correlation of Two Scales
#'
#' This function replicates \code{sim_style_data()} and calculates a scale score
#' (i.e., the mean across all items) for every person for every dimension. It
#' mainly returns two values for every replication sample, namely, the observed
#' correlation of the two scales as well as their partial correlation
#' controlling for response style.
#'
#' @param reps Desired number of replications
#' @param cor.cc Optional. Either fixed to the provided value or randomly
#'   sampled for each replication.
#' @param n Optional. Either fixed to the provided value or randomly sampled for
#'   each replication.
#' @param items Optional. Either fixed to the provided value or randomly sampled
#'   for each replication.
#' @param categ Optional. Either fixed to the provided value or randomly sampled
#'   for each replication.
#' @param reversed Optional. Either fixed to the provided value or randomly
#'   sampled for each replication.
#' @param var.s Optional. Either fixed to the provided value or randomly sampled
#'   for each replication.
#' @inheritParams sim_style_data
#' @return Returns a matrix of length \code{reps} containing the partial
#'   correlation \eqn{r_{X_1, X_2 \cdot \theta_{RS}}} in the first column and
#'   the observed correlation in the second column. Further columns contain the
#'   input parameters such as the number of categories.
#' @seealso \code{\link{sim_style_data}}
#' @export
rep_style_cor <- function(reps=1000, n=c(100, 1000), items=c(5, 10),
                          categ=c(3, 7), ndimc=2, style=NULL, reversed=c(0, .5),
                          var.s=c(0, 1), mu.s=c(-1, 1), cor.cc=c(-.7, .7), ...) {

  if (ndimc != 2) {
    stop("### Currently only defined for TWO content-related dimensions,
  ### please specify 'ndimc=2'.")
  }
  if (length(style) > 1) {
    stop("### Currently only defined for a single response style dimension.")
  }

  sty <- ifelse(is.null(style), FALSE, TRUE)

#   args.prov <- c(missing(n), missing(items), missing(categ), missing(reversed),
#                  missing(var.s), missing(mu.s), missing(cor.cc))

  res <- matrix(nrow=reps, ncol=12)
  colnames(res)  <- c("bias", "partial", "obs", "r.tt", "bias2", "n", "items", "categ", "rev",
                      "var.s", "mu.s", "cor.cc")


  pb <- txtProgressBar(min = 0, max = reps, style = 3, char="-")

  for (i in 1:reps) {
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

    if (sty) {
      if (length(var.s) == 1) {
        var.s.i <- var.s
      } else {
        var.s.i <- sample(seq(var.s[1], var.s[2], .1), 1)
      }

      if (length(mu.s) == 1) {
        mu.s.i <- mu.s
      } else {
        mu.s.i <- sample(seq(mu.s[1], mu.s[2], .1), 1)
      }
    }

    if (length(cor.cc) == 1) {
      cor.cc.i <- cor.cc
    } else {
      cor.cc.i <- round(sample(seq(cor.cc[1], cor.cc[2], .01), 1), 2)
    }

    if (sty) {
      d <- sim_style_data(n=n.i, items=items.i, categ=categ.i, ndimc=ndimc, style=style,
                          reversed=rev.i, var.s=var.s.i, mu.s=mu.s.i,
                          cor.cc=cor.cc.i, ...)
    } else {
      d <- sim_style_data(n=n.i, items=items.i, categ=categ.i, ndimc=ndimc, style=style,
                          reversed=rev.i, cor.cc=cor.cc.i, ...)
    }

    if (all(sty, d$response.style$var.style > 0)) {
      res[i, 2] <- partial_cor(cbind(rowMeans(d$dat[, , 1, drop=F]),
                                     rowMeans(d$dat[, , 2, drop=F]), d$theta[, 3]))
    } else {
      res[i, 2] <- cor(rowMeans(d$dat[, , 1, drop=F]), rowMeans(d$dat[, , 2, drop=F]))
    }
    res[i, 3] <- cor(rowMeans(d$dat[, , 1, drop=F]), rowMeans(d$dat[, , 2, drop=F]))
    res[i, 4] <- cor(d$theta[, 1], d$theta[, 2])
    res[i, c(6:12)] <- c(n.i, items.i, categ.i, rev.i, ifelse(sty, var.s.i, NA),
                         ifelse(sty, mu.s.i, NA), cor.cc.i)

    setTxtProgressBar(pb, i)
  }

  res[, 1] <- res[, 3] - res[, 2]
  res[, 5] <- res[, 3] - res[, 4]
  return(res)
}