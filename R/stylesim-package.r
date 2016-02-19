#' Simulate (and Analyse) Data Distorted by Response Styles
#' 
#' The most important function is \code{\link{sim_style_data}}, which allows to
#' simulate data distorted by response styles.
#' 
#' @references
#' Plieninger, H. (in press). Mountain or molehill: A simulation study on the impact of reponse styles. \emph{Educational and Psychological Measurement}.
#' 
#' @examples
#' \dontrun{
#' res_1 <- sim_style_data(n = 123, items = 12, categ = 5, style = "ARS", 
#'                         reversed = 6)
#' }
#'
#' @name stylesim
#' @docType package
#' @author Hansjoerg Plieninger
#' @importClassesFrom ROCR performance
NULL
