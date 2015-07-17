\name{predict_style}
\alias{predict_style}
\title{Compute TPR, FPR, etc. Comparing Data With and Without Response Style}
\usage{
predict_style(reps = 100, n = 250, items = 12, categ = 7, ndimc = 1,
  style = NULL, irtmodel = "RSM", reversed = 1/3, var.s, cor.cc,
  length.alpha = 99)
}
\arguments{
  \item{reps}{Desired number of replications}

  \item{length.alpha}{Desired number of cutoff- (alpha-)
  values at which to calculate TPR, FPR, etc.}

  \item{n}{Desired number of observations / persons}

  \item{items}{Desired number of items / stimuli. If
  \code{ndimc} > 1, this specifies the number of items per
  dimension.}

  \item{categ}{Desired number of categories per item}

  \item{ndimc}{Desired number of content-related latent
  variables (irrespective of number of style-related latent
  variables)}

  \item{style}{Parameter to specify which response style(s)
  influence the data, can be either numeric or character.
  Users may choose one or more among \code{"ERS1"} (e.g., 1
  / 0 / 0 / 0 / 1), \code{"ERS2"} (e.g., 2 / 1 / 0 / 1 /
  2), \code{"ARS"} (e.g., 0 / 0 / 0 / 1 / 1), \code{"ADRS"}
  (e.g, -1 / -1 / 0 / 1 / 1), and \code{"MRS"} (e.g., 0 / 0
  / 1 / 0 / 0). Alternatively, a user-specified vector of
  weights can be employed. Can also be \code{NULL}
  indicating complete abscence of response styles.}

  \item{irtmodel}{A character string indicating which model
  to use. Currently available is \code{"RSM"} (i.e., rating
  scale model).}

  \item{reversed}{Indicates the number of reverse-coded
  items. Can be either a ratio between 0 and 1 indicating
  the ratio of reverse-coded items or an integer between 0
  and \code{items}) indicating the number of reverse-coded
  items.}

  \item{var.s}{A vector indicating the variance of the
  style-related latent variables}

  \item{cor.cc}{An optional vector indicating the
  correlation between the content-related variables (if
  \code{ndimc} > 1). If \code{ndimc} > 2, \code{cor.cc} is
  recycled if it is not of length \code{ndimc*(ndimc-1) * 1
  / 2}.}
}
\value{
Returns a list.
}
\description{
This function simulates data under two parallel conditions,
one without response styles and one with response style(s).
The aim is to compare individual scale score shifts due to
response styles.
}
\seealso{
The replicated function \code{\link{sim_style_data}}. The
package \code{\link[ROCR]{}}.
}

