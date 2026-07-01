
###############################################################################
## fun_pava_comb.R
## Bivariate isotonic regression for J x K combination matrices
## Uses Iso::biviso for bivariate PAVA
###############################################################################

pava_comb <- function(y_mat, w_mat) {
  require(Iso)
  w_mat[w_mat <= 0] <- 0.01
  result <- Iso::biviso(y_mat, w_mat, warn = TRUE)[, ]
  return(result)
}
