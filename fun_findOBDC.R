###############################################################################
## fun_findOBDC.R
## Identify the true OBDC in a J x K dose matrix using RDS utility
###############################################################################

findOBDC_RDS <- function(toxM, effM, pT, qE, u11, u00) {
  J <- nrow(toxM); K <- ncol(toxM)
  RDS <- matrix(NA, J, K)
  admissible <- (toxM <= pT)
  for (j in 1:J) for (k in 1:K) {
    if (admissible[j, k]) {
      RDS[j, k] <- 100 * effM[j, k] * (1 - toxM[j, k]) +
        u00 * (1 - toxM[j, k]) * (1 - effM[j, k]) +
        u11 * toxM[j, k] * effM[j, k]
    }
  }
  if (all(is.na(RDS))) return(matrix(c(-1, -1), nrow = 1))
  if (max(effM[admissible], na.rm = TRUE) < qE) return(matrix(c(-1, -1), nrow = 1))
  best <- which(RDS == max(RDS, na.rm = TRUE), arr.ind = TRUE)
  storage.mode(best) <- "integer"
  return(best)
}