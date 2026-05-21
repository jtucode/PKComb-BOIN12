###############################################################################
## fun_PKCombBOIN12_one_v2.R
## Single-step dosing decision for PKComb-BOIN12
## Updated: Adjacent-only admissible sets with simplified PK branching
## Fixed: Untried doses compete via prior utility; removed random exploration
##
## v2 Fix: PK elimination now uses zeta1 (= 0.8 * r_P) instead of r_P as
##         the threshold. Previously, pnorm(r_P, ...) tested against the full
##         target PK (psi0PK), incorrectly eliminating doses whose PK was
##         above zeta1 but below psi0PK. This caused catastrophic failures
##         in scenarios where the OBDC had PK between zeta1 and psi0PK
##         (e.g., Scenarios 13 and 17).
###############################################################################
PKCombBOIN12_one <- function(doseDT, current_jk, pT, qE, J, K,
                             lambda_e, lambda_d, zeta1, N_star,
                             csize, decisionM, ub, u11, u00) {

  j <- current_jk[1]; k <- current_jk[2]
  idx_cur <- which(doseDT$j == j & doseDT$k == k)

  n <- doseDT$n[idx_cur]
  phat <- doseDT$pi_t_hat[idx_cur]
  rhat <- doseDT$r_d[idx_cur]
  cn <- max(1, round(n / csize))
  cn <- min(cn, nrow(decisionM))
  dN <- nrow(doseDT)

  find_idx <- function(jj, kk) {
    if (jj < 1 || jj > J || kk < 1 || kk > K) return(NA)
    idx <- which(doseDT$j == jj & doseDT$k == kk)
    if (length(idx) == 0) return(NA)
    if (doseDT$keep[idx] == 0) return(NA)
    return(idx)
  }

  ## Posterior utility: untried doses (n=0) get prior Beta(1,1) value
  post_util <- function(idx) {
    if (is.na(idx)) return(-Inf)
    xd <- doseDT$x_d[idx]; nn <- doseDT$n[idx]
    ## When n=0, x_d=0: returns 1 - pbeta(ub, 1, 1) = 1 - ub (prior)
    return(1 - pbeta(ub, 1 + xd, nn + 1 - xd))
  }

  ## ---- Safety elimination ---- ##
  if (!is.na(decisionM$DU_T[cn])) {
    if (phat * n >= decisionM$DU_T[cn]) {
      for (ii in 1:dN) {
        if (doseDT$j[ii] >= j & doseDT$k[ii] >= k) {
          if (is.na(doseDT$elim_reason[ii])) doseDT$elim_reason[ii] <- "SAFETY"
          doseDT$keep[ii] <- 0
        }
      }
      cand_idx <- c(find_idx(j - 1, k), find_idx(j, k - 1))
      cand_idx <- cand_idx[!is.na(cand_idx)]
      if (length(cand_idx) == 0) return(list(doseDT = doseDT, newdose = NA))
      posts <- sapply(cand_idx, post_util)
      best <- cand_idx[which.max(posts + seq_along(posts) * 1e-6)]
      return(list(doseDT = doseDT, newdose = c(doseDT$j[best], doseDT$k[best])))
    }
  }

  ## ---- Efficacy elimination ---- ##
  if (!is.na(decisionM$DU_E[cn])) {
    qhat <- doseDT$pi_e_hat[idx_cur]
    if (qhat * n <= decisionM$DU_E[cn]) {
      if (is.na(doseDT$elim_reason[idx_cur])) doseDT$elim_reason[idx_cur] <- "EFFICACY"
      doseDT$keep[idx_cur] <- 0
    }
  }

  ## ---- PK elimination ---- ##
  ## FIX (v2): Use zeta1 as the elimination threshold, not r_P = psi0PK.
  ## zeta1 = 0.8 * psi0PK is the PK sufficiency boundary used throughout
  ## the design (flowchart dose-escalation and OBDC admissible set).
  ## The elimination should be consistent: eliminate only when we are
  ## confident the dose's PK is below zeta1.
  if (n >= 6 && doseDT$r_sd[idx_cur] > 0) {
    pk_elim <- pnorm(zeta1, mean = rhat, sd = doseDT$r_sd[idx_cur] / sqrt(n)) > 0.95

    if (pk_elim) {
      if (j == J && k == K) {
        ## d = D: eliminate all dose levels
        for (ii in 1:dN) {
          if (is.na(doseDT$elim_reason[ii])) doseDT$elim_reason[ii] <- "PK"
        }
        doseDT$keep <- 0
        return(list(doseDT = doseDT, newdose = NA))
      } else {
        ## 2 <= d < D: eliminate the lowest available dose level
        for (ii in 1:dN) {
          if (doseDT$keep[ii] == 1) {
            if (is.na(doseDT$elim_reason[ii])) doseDT$elim_reason[ii] <- "PK"
            doseDT$keep[ii] <- 0; break
          }
        }
      }
    }
  }

  if (sum(doseDT$keep) == 0) return(list(doseDT = doseDT, newdose = NA))

  ## ---- Build candidate set (Updated Flowchart) ---- ##
  cands <- list()

  if (phat >= lambda_d) {
    ## CASE 1: Toxic — de-escalate (no PK branch)
    cands <- list(c(j - 1, k), c(j, k - 1))

  } else if (phat > lambda_e && n >= N_star) {
    ## CASE 2: Intermediate + saturated — exploit (no PK branch)
    cands <- list(c(j - 1, k), c(j, k), c(j, k - 1))

  } else if (phat > lambda_e && n < N_star) {
    ## CASE 3: Intermediate + unsaturated — full 5 neighbors (no PK branch)
    cands <- list(c(j-1,k), c(j,k), c(j,k-1), c(j+1,k), c(j,k+1))

  } else {
    ## CASE 4: Safe (p_hat <= lambda_e) — PK branch
    if (rhat > zeta1) {
      cands <- list(c(j-1,k), c(j,k), c(j,k-1), c(j+1,k), c(j,k+1))
    } else {
      cands <- list(c(j,k), c(j+1,k), c(j,k+1))
    }
  }

  ## ---- Filter valid candidates ---- ##
  valid <- list(); seen <- c()
  for (cc in cands) {
    idx <- find_idx(cc[1], cc[2])
    key <- paste(cc[1], cc[2])
    if (!is.na(idx) && !(key %in% seen)) {
      valid[[length(valid) + 1]] <- cc
      seen <- c(seen, key)
    }
  }

  if (length(valid) == 0) {
    if (doseDT$keep[idx_cur] == 1) return(list(doseDT = doseDT, newdose = current_jk))
    rem <- which(doseDT$keep == 1)
    if (length(rem) == 0) return(list(doseDT = doseDT, newdose = NA))
    return(list(doseDT = doseDT, newdose = c(doseDT$j[rem[1]], doseDT$k[rem[1]])))
  }

  ## ---- Select dose using RDS (all candidates including untried) ---- ##
  posts <- rep(-Inf, length(valid))
  for (i in seq_along(valid)) {
    idx <- find_idx(valid[[i]][1], valid[[i]][2])
    if (!is.na(idx)) posts[i] <- post_util(idx) + runif(1, 0, 1e-6)
  }

  if (max(posts) == -Inf) {
    if (doseDT$keep[idx_cur] == 1) return(list(doseDT = doseDT, newdose = current_jk))
    return(list(doseDT = doseDT, newdose = NA))
  }
  best_i <- which.max(posts)
  return(list(doseDT = doseDT, newdose = valid[[best_i]]))
}
