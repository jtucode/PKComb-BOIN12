
###############################################################################
## fun_PKCombBOIN12_OBDC.R
## Final OBDC selection after trial completion
## Uses isotonic regression on toxicity and PK, logistic regression for efficacy
###############################################################################

fun_PKCombBOIN12_OBDC <- function(doseDT, J, K, pT, u11, u00, zeta1,true_admissible) {
  n_mat <- x_mat <- y_mat <- r_mat <- keep_mat <- matrix(0, J, K)
  for (idx in 1:nrow(doseDT)) {
    jj <- doseDT$j[idx]; kk <- doseDT$k[idx]
    n_mat[jj, kk] <- doseDT$n[idx]
    x_mat[jj, kk] <- doseDT$x[idx]
    y_mat[jj, kk] <- doseDT$y[idx]
    r_mat[jj, kk] <- doseDT$r_d[idx]
    keep_mat[jj, kk] <- doseDT$keep[idx]
  }
  
  ## Step 1: MTD via isotonic regression on toxicity
  phat <- (x_mat + 0.05) / (n_mat + 0.1)
  phat[n_mat == 0] <- pT
  ptilde <- pava_comb(phat, n_mat + 0.1)
  # tie-breaking: small increment for higher doses
  ptilde <- ptilde + 0.001 * (row(ptilde) + col(ptilde))
  ptilde_search <- ptilde; ptilde_search[n_mat == 0] <- 10
  MTD_idx <- which(abs(ptilde_search - pT) == min(abs(ptilde_search - pT)), arr.ind = TRUE)
  if (nrow(MTD_idx) > 1) MTD_idx <- MTD_idx[1, , drop = FALSE]
  MTD_j <- MTD_idx[1]; MTD_k <- MTD_idx[2]
  
  ## Step 2: PK isotonic regression
  r_iso <- r_mat; r_iso[n_mat == 0] <- 0
  w_pk <- n_mat; w_pk[w_pk == 0] <- 0.01
  rtilde <- pava_comb(r_iso, w_pk)
  
  ## Step 3: Utility for OBDC selection
  ubar <- matrix(-100, J, K)
  for (jj in 1:J) for (kk in 1:K) {
    if (n_mat[jj, kk] > 0 && keep_mat[jj, kk] == 1) {
      ubar[jj, kk] <- (u11 * y_mat[jj, kk] + u00 * (n_mat[jj, kk] - x_mat[jj, kk])) / 100
      ubar[jj, kk] <- (ubar[jj, kk] + 1) / (n_mat[jj, kk] + 2)
    }
  }
  
  ## Step 4: Restrict admissible set
  Addm<-matrix(rep(1,J*K), nrow = J, ncol = K, byrow = TRUE) #Check featured admissiblility not vanishing true admissible set
  for (jj in 1:J) for (kk in 1:K) {
    # Exclude above MTD (doses where BOTH indices exceed MTD)
    above_mtd <- (jj > MTD_j & kk >= MTD_k) | (jj >= MTD_j & kk > MTD_k)
    if (above_mtd && !(jj == MTD_j && kk == MTD_k)) {
      if (n_mat[jj, kk] > 0 && (ptilde[jj, kk] - ptilde[MTD_j, MTD_k]) >= 0.05) {
        ubar[jj, kk] <- -100
        Addm[jj, kk] <- 0
      }
    }
    # Exclude PK-insufficient doses
    if (n_mat[jj, kk] >= 3 && rtilde[jj, kk] > 0 && rtilde[jj, kk] < zeta1) {
      ubar[jj, kk] <- -100
      Addm[jj, kk] <- 0
    }
    if (keep_mat[jj, kk] == 0) 
    {
      ubar[jj, kk] <- -100
      Addm[jj, kk] <- 0
    }
    if (n_mat[jj, kk] == 0) 
    {
      ubar[jj, kk] <- -100
      Addm[jj, kk] <- 0
    }
  }
  


  add_score<-sum(true_admissible & (Addm!=0))
  
  if (max(ubar) <= -100) return(list(OBDC=c(99, 99),Addmissible_score=add_score,  Utility_vals=ubar))
  best <- which(ubar == max(ubar), arr.ind = TRUE)
  if (nrow(best) > 1) best <- best[1, , drop = FALSE]
  return(list(OBDC=as.integer(best),Addmissible_score=add_score, Utility_vals=ubar))
}

