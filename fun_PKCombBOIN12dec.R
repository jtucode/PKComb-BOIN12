###############################################################################
## fun_PKCombBOIN12dec.R
## Decision table for dose elimination boundaries
## Identical structure to PKBOIN-12 but used per-combination
###############################################################################

fun_PKCombBOIN12dec <- function(pT, qE, lambda_e, lambda_d, csize, cN,
                                cutoff_tox = 0.95, cutoff_eff = 0.90) {
  # Returns data.frame with columns: n, DU_T, DU_E, D, E
  # DU_T: min #DLTs to eliminate dose for toxicity
  # DU_E: max #responses below which dose eliminated for futility
  
  result <- data.frame(
    n = (1:cN) * csize,
    DU_T = rep(NA, cN),
    DU_E = rep(NA, cN),
    D = rep(NA, cN),
    E = rep(NA, cN)
  )
  result$D <- ceiling(result$n * lambda_d)
  result$E <- floor(result$n * lambda_e)
  
  for (i in 1:cN) {
    #cutoff_tox<-ifelse(i<=1,0.85,cutoff_tox) #only for the first cohort
    # Toxicity elimination boundary
    if (1 - pbeta(pT, result$n[i] + 1, 1) >= cutoff_tox) {
      result$DU_T[i] <- min(which(
        sapply(0:result$n[i], function(x)
          1 - pbeta(pT, x + 1, result$n[i] + 1 - x)) >= cutoff_tox
      )) - 1
    }
    # Efficacy elimination boundary
    if (pbeta(qE, 1, result$n[i] + 1) >= cutoff_eff) {
      result$DU_E[i] <- max(which(
        sapply(0:result$n[i], function(x)
          pbeta(qE, x + 1, result$n[i] + 1 - x)) >= cutoff_eff
      )) - 1
    }
  }
  return(result)
}
