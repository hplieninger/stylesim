# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# SIMULATE DATA FOR ROC ANALYSES --------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Simulate Data Distorted by Response Style and Calculate, for Example, the
#' True Positive Rate for Subsequent ROC Analyses
#' 
#' This function repeatedly simulates data that are distorted by response style. It makes
#' use of the scale scores \eqn{\bar{X}}---which are influenced by both the attribute and
#' response style---as well as the true person parameters \eqn{\theta_1}, that is, the
#' persons' standing on the attribute, which is unaffected by response style.
#' These two scores are used to classify each individual at certain cutoff
#' values as either above (i.e., positive) or below (i.e., negative) the cutoff.
#' Subsequently, the true positive rate (TPR), the false positive rate (FPR),
#' the precision, and the accuracy are calculated at every cutoff value.
#' 
#' @param reps Desired number of replications
#' @param cutoff Optional numeric vector. Desired cutoff- (alpha-) values at which
#'   to calculate TPR, FPR, etc; probably a sequence, for example from .95 to
#'   .05 (must be within the open interval of (0, 1))
#' @param ... Other parameters passed to \code{\link{sim_style_data}}.
#' @inheritParams sim_style_data
#' @inheritParams rep_style_rel
#' @return
#'   \item{alpha}{The cutoff-/alpha-values}
#'   \item{TPR}{True positive rate}
#'   \item{FPR}{False positive rate}
#'   \item{PREC}{Precision}
#'   \item{ACC}{Accuracy}
#'   \item{}{Further elements containing the input parameters}
#' @seealso The replicated function \code{\link{sim_style_data}}.
#' @export
predict_style <- function(reps = 100, n = 200, items = 10, categ = 5, ndimc = 1,
                          style = NULL, reversed = 1/3, mu.s = 0, var.s,
                          sig = NULL, emp = TRUE, cutoff, ...) {
    time.1 <- Sys.time()
    
    pb <- txtProgressBar(min = 0, max = reps, style = 3, char="-")  
    
    if (ndimc != 1) stop("Only defined for a single content-related dimension,
    please specify ndimc=1.")
    if (missing(style)) stop("Please specify the argument 'style'")
    if ("df" %in% names(list(...))) stop("Argument 'df' currently not implemented, i.e., it has no effect.")
    
    d <- array(dim = c(n, reps, 2))  
    
    for(ii in 1:reps){
        
        dx <- sim_style_data(n = n, items = items, categ = categ, ndimc = ndimc,
                             style = style, reversed = reversed, var.s = var.s,
                             sig = sig, emp = emp, ...)
        
        d[, ii, 1] <- dx$theta[, 1]
        d[, ii, 2] <- rowMeans(dx$dat[, , 1])
        
        setTxtProgressBar(pb, ii)
    }
    close(pb)
    
    if (missing(cutoff)) {
        cutoff <- seq(1, 0, length = items * (categ - 1) + 1)
        cutoff <- cutoff[-c(1, length(cutoff))]
    }
    cutoff.l <- length(cutoff)
    cutoff <- array(rep(cutoff, each = (n*reps)), dim = c(n, reps, cutoff.l))
    
    d.the <- apply(d[, , 1], 2, function(x) (ecdf(x)(x)))
    d.obs <- apply(d[, , 2], 2, function(x) (ecdf(x)(x)))
    
    d.the <- array(d.the, dim = c(n, reps, cutoff.l))  # ORIGINAL DATA
    d.obs <- array(d.obs, dim = c(n, reps, cutoff.l))  # OBSERVED DATA
    
    c.the <- (d.the >= cutoff) + 0  # ORIGINAL CLASSES
    c.obs <- (d.obs >= cutoff) + 0  # OBSERVED CLASSES
    
    tp <- colSums(c.the==1 & c.obs==1)
    fp <- colSums(c.the==0 & c.obs==1)
    fn <- colSums(c.the==1 & c.obs==0)
    tn <- colSums(c.the==0 & c.obs==0)
    
    pos <- tp + fn
    neg <- fp + tn
    
    tpr <- tp/pos
    fpr <- fp/neg
    pre <- tp/(tp + fp)
    acc <- (tp + tn)/(pos + neg)

    styleweights <- sim_style_data(n = n, items = items, categ = categ, ndimc = ndimc,
                                   style = style, reversed = reversed, var.s = var.s,
                                   ...)$response.style
    time.2 <- difftime(Sys.time(), time.1, units = "hours")
    res.list <- list(alpha.values = cutoff[1, , ],
                     TPR = tpr, FPR = fpr, PREC = pre, ACC = acc,
                     style = styleweights[1:2],
                     setup = list(df_of_Wishart = df, dims.content = ndimc),
                     time.elapsed = time.2)
        
    cat("Total time elapsed:", round(time.2[[1]], 3), "hours\n")
    return(res.list)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# REARRANGE SIMULATION RESULTS FOR ROC PLOTS --------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Rearrange Simulation Results in Order to Plot Them
#' 
#' This function takes as input the output from \code{\link{predict_style}} and
#' rearranges this output so that it can subsequently be used by the plotting
#' function of the \strong{ROCR} package (\code{\link[ROCR]{plot.performance}}). Conversely to traditional
#' ROC-analyses, which employ two measures (e.g., TPR and FPR), this function
#' returns output to plot only one of these measures (e.g., TPR).
#' 
#' @param x Character. One of \code{"TPR"}, \code{"FPR"}, \code{"ACC"}, or \code{"PREC"}.
#' @param dat List. The output from \code{\link{predict_style}}.
#' @return An S4 object of the ROCR-class performance.
#' @examples
#' \dontrun{
#' res_1 <- predict_style(reps = 5, style = "ARS")
#' res_2 <- performance_style("ACC", res_1)
#' ROCR::plot(res_2)
#' ROCR::plot(res_2, avg="vertical", spread.estimate="stddev",
#'            show.spread.at=seq(.1, .9, length=9))
#' }
#' @seealso The package \strong{ROCR} and its \code{\link[ROCR]{performance-class}}.
#' @export
# @importClassesFrom ROCR performance
performance_style <- function(x, dat){
  
  res <- new("performance")
  res@alpha.name <- "Cutoff"
  res@alpha.values <- lapply(seq_len(nrow(dat$alpha.values)),
                             function(i) dat$alpha.values[i,])  
  res@x.name <- "Cutoff"
  res@x.values <- lapply(seq_len(nrow(dat$alpha.values)),
                         function(i) dat$alpha.values[i,])
  
  if(x == "FPR"){
    res@y.name <- "False positive rate"
    res@y.values <- lapply(seq_len(nrow(dat$FPR)), function(i) dat$FPR[i,])
  }
  if(x == "ACC"){
    res@y.name <- "Accuracy"
    res@y.values <- lapply(seq_len(nrow(dat$ACC)), function(i) dat$ACC[i,]) 
  }
  if(x == "TPR"){
    res@y.name <- "True positive rate"
    res@y.values <- lapply(seq_len(nrow(dat$TPR)), function(i) dat$TPR[i,]) 
  }
  if(x == "PREC"){
    res@y.name <- "Precision"
    res@y.values <- lapply(seq_len(nrow(dat$PREC)), function(i) dat$PREC[i,]) 
  }
  
  return(res)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# DISCARDED: SIMULATE DATA FOR ROC ANALYSES ---------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# In an earlier stage, the function 'predict_style()' internally called this
# function, but now 'predict_style()' simulates the data itself/directly.

# Simulate Data for Subsequent ROC Analyses
# 
# This function simulates data that are distorted by response style. It outputs
# the scale scores X---which are influenced by both the attribute and response
# style---as well as the true person parameters theta_1, that is, the persons'
# standing on the attribute, which is unaffected by response style.
# 
# @param reps Desired number of replications
# @inheritParams sim_style_data
# @return Returns an array of dimension \code{n} x \code{reps} x 2. The first 
#   matrix contains the true person parameters theta_1 for every person for
#   every replication, and the second contains the corresponding scale scores
#   (i.e., the mean across all items).
# @seealso The replicated function \code{\link{sim_style_data}}.
# @export
# rep_style_roc <- function(reps = 100, n = 200, items = 10, categ = 5, ndimc = 1,
#                           style = NULL, reversed = 1/3, mu.s = 0, var.s,
#                           df = 10, sig = NULL, emp = TRUE, ...) {
#     
#     pb <- txtProgressBar(min = 0, max = reps, style = 3, char="-")  
#     
#     if (ndimc != 1) stop("Only defined for a single content-related dimension,
#     please specify ndimc=1.")
#     if (missing(style)) stop("Please specify the argument 'style'")
#     
#     res <- array(dim=c(n, reps, 2))  
#     
#     for(i in 1:reps){
#         
#         d <- sim_style_data(n = n, items = items, categ = categ, ndimc = ndimc,
#                             style = style, reversed = reversed, var.s = var.s,
#                             ...)
#         
#         res[, i, 1] <- d$theta[, 1]
#         res[, i, 2] <- rowMeans(d$dat[, , 1])
#         
#         setTxtProgressBar(pb, i)
#     }
#     return(res)
# }