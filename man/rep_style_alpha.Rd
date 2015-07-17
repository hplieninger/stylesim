\name{rep_style_alpha}
\alias{rep_style_alpha}
\title{Replicate a Style Simulation and Investigate Effect on Cronbach's Alpha}
\usage{
rep_style_alpha(reps = 100, n = 250, items = 12, categ = 7, ndimc = 1,
  style = NULL, irtmodel = "RSM", reversed = 1/3, var.s, cor.cc)
}
\arguments{
  \item{reps}{Desired number of replications}

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
Returns a matrix of length \code{reps} containing the true
reliability in the first column and Cronbach's Alpha in the
second column.  The true reliability is calculated via
\eqn{r^2_{\theta_1, \theta_2}}.
}
\description{
This function replicates \code{sim_style_data()} and
returns the true reliability as well as observed Cronbach's
Alpha for every replication sample.
}
\seealso{
The replicated function \code{\link{sim_style_data}}.
Cronbach's alpha: \code{\link[ltm]{cronbach.alpha}}
}

