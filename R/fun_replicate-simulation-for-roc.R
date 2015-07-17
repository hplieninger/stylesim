# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# SIMULATE DATA FOR ROC ANALYSES --------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Simulate Data for Subsequent ROC Analyses
#' 
#' This function simulates data under two parallel conditions, one without
#' response styles and one with response style(s). The aim is to compare
#' individual scale score shifts due to response styles.
#' 
#' @param reps Desired number of replications
#' @inheritParams sim_style_data
#' @return Returns an array of dimension \code{n} × \code{reps} × 2. The first
#'   matrix contains a scale score (i.e., the mean across all items) for every
#'   person for every replication if response styles were ABSENT, and the second
#'   matrix contains the data when response styles are PRESENT.
#' @seealso The replicated function \code{\link{sim_style_data}}. The package
#'   \code{\link[ROCR]{}}.
#' @export
rep_style_roc <- function(reps=100, n=200, items=10, categ=5, ndimc=1,
                          style=NULL, irtmodel="RSM", reversed=1/3, var.s,
                          cor.cc, ...) {
  
  pb <- txtProgressBar(min = 0, max = reps, style = 3, char="-")  
  
  if (ndimc != 1) stop("Only defined for a single content-related dimension,
    please specify ndimc=1.")
  if (missing(style)) stop("Please specify the argument 'style'")
  
  res <- array(dim=c(n, reps, 2))  
  
  for(i in 1:reps){

    d <- sim_style_data(n=n, items=items, categ=categ, ndimc=ndimc, style=style,
                        irtmodel=irtmodel, reversed=reversed, var.s=var.s)
    
    res[, i, 1] <- d$theta[, 1]
    res[, i, 2] <- rowMeans(d$dat[, , 1])
    
    setTxtProgressBar(pb, i)
  }
  return(res)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# SIMULATE DATA FOR ROC ANALYSES --------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#' Compute TPR, FPR, etc. Comparing Data With and Without Response Style
#' 
#' This function simulates data under two parallel conditions, one without 
#' response styles and one with response style(s). The aim is to compare 
#' individual scale score shifts due to response styles.
#' 
#' @param reps Desired number of replications
#' @param length.alpha Desired number of cutoff- (alpha-) values at which to
#'   calculate TPR, FPR, etc.
#' @inheritParams sim_style_data
#' @return Returns a list.
#' @seealso The replicated function \code{\link{sim_style_data}}. The package 
#'   \code{\link[ROCR]{}}.
#' @export
predict_style <- function(reps=100, n=200, items=10, categ=5, ndimc=1,
                          style=NULL, irtmodel="RSM", reversed=1/3, var.s,
                          ...){
  
  if (missing(style)) stop("Please specify the argument 'style'")
  
  d <- rep_style_roc(reps=reps, n=n, items=items, categ=categ, ndimc=ndimc,
                     style=style, irtmodel=irtmodel, reversed=reversed,
                     var.s=var.s)
  
  rate <- seq(1, 0, length = items * (categ - 1) + 1)
  rate <- rate[-c(1, length(rate))]
  rate.l <- length(rate)
  rate <- array(rep(rate, each=(n*reps)), dim=c(n, reps, length(rate)))
  
#   q1 <- apply(d[, , 2], 2, quantile, rate) # QUANTILES FROM OBSERVED DATA
#   q1 <- array(rep(t(q1), each=n), dim=c(n, reps, length(rate)))  # AS ARRAY
  
  d.the <- apply(d[, , 1], 2, function(x) (ecdf(x)(x)))
  d.obs <- apply(d[, , 2], 2, function(x) (ecdf(x)(x)))

  d.the <- array(d.the, dim=c(n, reps, rate.l))  # ORIGINAL DATA
  d.obs <- array(d.obs, dim=c(n, reps, rate.l))  # OBSERVED DATA  
  
  c.the <- ifelse(d.the >= rate, 1, 0)  # ORIGINAL CLASSES
  c.obs <- ifelse(d.obs >= rate, 1, 0)  # OBSERVED CLASSES
  
  tp <- colSums(c.the==1 & c.obs==1)
  fp <- colSums(c.the==0 & c.obs==1)
  fn <- colSums(c.the==1 & c.obs==0)
  tn <- colSums(c.the==0 & c.obs==0)
  
  p <- tp + fn
  n <- fp + tn
  
  tpr <- tp/p
  fpr <- fp/n
  pre <- tp/(tp + fp)
  acc <- (tp + tn)/(p + n)  
  
  return(list(alpha.values = rate[1, , ],
              TPR=tpr, FPR=fpr, PREC=pre, ACC=acc))
}

PerformanceForROC <- function(x, dat){
  library("ROCR")
  
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