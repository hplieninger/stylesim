\name{sim_style_data}
\alias{sim_style_data}
\title{Simulate Data From a Multidimensional Rasch Model}
\usage{
sim_style_data(n = 250, items = 12, categ = 7, ndimc = 1,
  style = NULL, irtmodel = "RSM", reversed = 1/3, var.s, cor.cc, seed,
  pop.thres = FALSE, normed = TRUE)
}
\arguments{
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

  \item{seed}{Optional argument to specify a seed for all
  random number generators (e.g., for the latent
  variables)}

  \item{pop.thres}{Logical. If \code{TRUE}, the thresholds
  of the RSM are set to their expected value (i.e.,
  \code{seq(-1.5, 1.5, length=categ - 1)})}

  \item{normed}{Logical. If \code{TRUE}, the weights for
  the person parameters (e.g., 0 / 1 / 2 / 3 / 4 or 2 / 1 /
  0 / 1 / 2) are centered and normed (i.e., \eqn{||x|| =
  1}).}
}
\value{
Returns a list containing the data and basically a summary
of the input specifications with elements \describe{

\item{dat}{an array containing the data / responses of
format \code{n} x \code{items} x \code{ndimc}}

\item{theta}{a matrix containing the true person parameters
of format \code{n} x \code{(ndimc+ ength(style))}. For
example, if \code{ndimc=1}, and \code{style=c("ERS1",
"ARS")}, then the second column corresponds to "ERS1" and
the third column to "ARS".}

\item{item.parameters}{a vector containing the item
parameters. The first \code{items} elements are the
location parameters and the last \code{(categ-1)} elements
are the threshold parameters.}

\item{n}{number of observations}

\item{items.per.dimension}{number of items per dimension}

\item{reverse.coded.items}{number of reverse-coded items
per dimension}

\item{categories}{number of categories}

\item{irtmodel}{type of IRT model}

\item{dims.content}{content-related dimensions}

\item{var.style}{optional. Variance of the style-related
latent variable(s).}

\item{c.c.cor}{optional. Correlation of the content-related
latent variables.}

\item{response.style}{\describe{ \item{the specified
response style(s)}{} \item{the implied coding scheme}{}}}

}
}
\description{
This function allows to simulate data from a
multidimensional Rasch model, in which one (or more) latent
variables represent response styles such as extreme
response style or acquiescence.
}

