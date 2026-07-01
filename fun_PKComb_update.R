###############################################################################
## fun_PKComb_update.R
## Update cohort: generate PK-linked outcomes, compute TITE statistics
## Extends fun_TITE_PK_update to J x K combination dose matrix
##
## v3 FIX: (1) Removed double-counted accrual gap in time_next.
##         (2) Aligned enrollment generation with CB12-TITE: first patient
##             enrolls at time_current, subsequent patients spaced randomly.
##         (3) Non-TITE gap moved inside max() for consistency with TITE.
###############################################################################

fun_PKComb_update <- function(cid, CV, g_P, patDT, doseDT, current_jk,
                              time_current, tox_win, eff_win, csize,
                              tox_dist, eff_dist, accrual, susp = 0.5,
                              u11 = 60, u10 = 0, u01 = 100, u00 = 40,
                              use_susp = TRUE, accrual_random = FALSE,
                              pqcorr = 0) {
  # current_jk: c(j,k) current dose combination
  # doseDT: data.frame with columns: j, k, PK, tox, eff, n, x, y, keep,
  #         ESS_t, ESS_e, pi_t_hat, pi_e_hat, x_d, r_d, r_sd
  # patDT: patient-level data.frame
  
  #require(truncnorm)
  dN <- nrow(doseDT)
  j_cur <- current_jk[1]; k_cur <- current_jk[2]
  # Find row index in doseDT for current combination
  idx_cur <- which(doseDT$j == j_cur & doseDT$k == k_cur)
  
  ### Step 1: Generate outcomes for current cohort ###
  
  # 1.1 Enrollment time
  # First patient enrolls at time_current (matching CB12-TITE convention).
  # Subsequent patients spaced by inter-arrival times.
  rows_cid <- which(patDT$cid == cid)
  if (accrual_random) {
    inter_arrival <- c(0, cumsum(runif(csize - 1, 0, 2 * accrual)))
    patDT$enroll[rows_cid] <- time_current + inter_arrival
  } else {
    patDT$enroll[rows_cid] <- time_current + (0:(csize - 1)) * accrual
  }
  
  # 1.2 Assessment end times
  patDT$toxend[rows_cid] <- patDT$enroll[rows_cid] + tox_win
  patDT$effend[rows_cid] <- patDT$enroll[rows_cid] + eff_win
  
  # 1.3 Generate individual PK values
  pk_mean_cur <- doseDT$PK[idx_cur]
  patDT$PK[rows_cid] <- truncnorm::rtruncnorm(csize, a = 0, b = Inf,
                                              mean = pk_mean_cur,
                                              sd = pk_mean_cur * CV)
  
  # 1.4 Individual-level toxicity/efficacy probabilities (PK-linked)
  pk_vals <- patDT$PK[rows_cid]
  patDT$tox_prob[rows_cid] <- pmin(pmax(
    doseDT$tox[idx_cur] * (1 + g_P * (pk_vals - pk_mean_cur) / pk_mean_cur), 0), 1)
  patDT$eff_prob[rows_cid] <- pmin(pmax(
    doseDT$eff[idx_cur] * (1 + g_P * (pk_vals - pk_mean_cur) / pk_mean_cur), 0), 1)
  
  # 1.5 Generate binary outcomes
  if (pqcorr == 0) {
    patDT$toxobs[rows_cid] <- rbinom(csize, 1, patDT$tox_prob[rows_cid])
    patDT$effobs[rows_cid] <- rbinom(csize, 1, patDT$eff_prob[rows_cid])
  } else {
    jp <- pqcorr * sqrt(patDT$tox_prob[rows_cid] * (1 - patDT$tox_prob[rows_cid]) *
                          patDT$eff_prob[rows_cid] * (1 - patDT$eff_prob[rows_cid])) +
      patDT$tox_prob[rows_cid] * patDT$eff_prob[rows_cid]
    rand <- runif(csize)
    patDT$toxobs[rows_cid] <- as.numeric(rand <= patDT$tox_prob[rows_cid])
    patDT$effobs[rows_cid] <- as.numeric(
      (rand <= jp) | (rand > 1 - patDT$eff_prob[rows_cid] + jp))
  }
  
  # 1.6 Dose assignment
  patDT$dj[rows_cid] <- j_cur
  patDT$dk[rows_cid] <- k_cur
  
  # 1.7 Generate event times
  new_tox <- sum(patDT$toxobs[rows_cid] == 1)
  new_eff <- sum(patDT$effobs[rows_cid] == 1)
  
  tox_idx <- rows_cid[patDT$toxobs[rows_cid] == 1]
  eff_idx <- rows_cid[patDT$effobs[rows_cid] == 1]
  
  if (tox_dist == "Uniform" && new_tox > 0) {
    patDT$toxdat[tox_idx] <- patDT$enroll[tox_idx] + runif(new_tox, max = tox_win)
  }
  if (eff_dist == "Uniform" && new_eff > 0) {
    patDT$effdat[eff_idx] <- patDT$enroll[eff_idx] + runif(new_eff, max = eff_win)
  }
  
  # 1.8 Confirmation times
  patDT$tox_confirm[rows_cid] <- pmin(patDT$toxend[rows_cid], patDT$toxdat[rows_cid])
  patDT$eff_confirm[rows_cid] <- pmin(patDT$effend[rows_cid], patDT$effdat[rows_cid])
  
  ### Step 2: Determine time_next ###
  #
  # time_next is the earliest time at which the next dose decision can be made.
  # It must be >= the last enrollment time (can't decide before all current
  # patients are enrolled). Under TITE (use_susp=TRUE), it must also be >=
  # the time at which the suspension fraction of patients at the current dose
  # have completed their assessment. Under non-TITE (use_susp=FALSE), it must
  # be >= the time at which ALL patients at the current dose have completed.
  #
  # NOTE: No additional accrual gap is added here. The next cohort's enrollment
  # times are generated in the NEXT call to this function starting at time_next,
  # with their own inter-arrival spacing. Adding a gap here would double-count.
  
  enrolled <- patDT$enroll[!is.na(patDT$enroll)]
  last_enroll <- max(enrolled)
  
  tmp <- patDT[!is.na(patDT$dj) & patDT$dj == j_cur & patDT$dk == k_cur, ]
  n_d <- nrow(tmp)
  
  if (use_susp) {
    # TITE: wait until susp fraction of patients at current dose have completed
    min_n <- min(floor(n_d * susp + 1), n_d)
    tox_next <- sort(tmp$tox_confirm)[min_n]
    eff_next <- sort(tmp$eff_confirm)[min_n]
    time_next <- max(tox_next, eff_next, last_enroll)
  } else {
    # Non-TITE: wait until ALL patients at current dose have completed
    tox_next <- max(tmp$tox_confirm)
    eff_next <- max(tmp$eff_confirm)
    time_next <- max(tox_next, eff_next, last_enroll)
  }
  
  ### Step 3: Update doseDT statistics ###
  
  patDT_current <- patDT[!is.na(patDT$dj), ]
  patDT_current$delta_t <- -1
  patDT_current$delta_e <- -1
  patDT_current$t_i <- time_next - patDT_current$enroll
  
  # Toxicity observation status
  patDT_current$delta_t[patDT_current$toxdat <= time_next] <- 1
  patDT_current$delta_t[patDT_current$toxend <= time_next &
                          patDT_current$toxdat > time_next] <- 0
  
  # Efficacy observation status
  patDT_current$delta_e[patDT_current$effdat <= time_next] <- 1
  patDT_current$delta_e[patDT_current$effend <= time_next &
                          patDT_current$effdat > time_next] <- 0
  
  # TITE weights
  patDT_current$w_t <- ifelse(
    patDT_current$t_i < tox_win & patDT_current$delta_t < 0,
    patDT_current$t_i / tox_win, 0)
  patDT_current$w_e <- ifelse(
    patDT_current$t_i < eff_win & patDT_current$delta_e < 0,
    patDT_current$t_i / eff_win, 0)
  
  # Update per-dose statistics
  for (idx in 1:dN) {
    jj <- doseDT$j[idx]; kk <- doseDT$k[idx]
    mask <- which(patDT_current$dj == jj & patDT_current$dk == kk)
    if (length(mask) > 0) {
      doseDT$n[idx] <- length(mask)
      delta_t_i <- patDT_current$delta_t[mask]
      delta_e_i <- patDT_current$delta_e[mask]
      
      doseDT$x[idx] <- sum(delta_t_i == 1)
      doseDT$y[idx] <- sum(delta_e_i == 1)
      
      w_t_i <- patDT_current$w_t[mask]
      w_e_i <- patDT_current$w_e[mask]
      
      doseDT$ESS_t[idx] <- sum(delta_t_i == 1) + sum(delta_t_i == 0) + sum(w_t_i)
      doseDT$pi_t_hat[idx] <- doseDT$x[idx] / doseDT$ESS_t[idx]
      doseDT$ESS_e[idx] <- sum(delta_e_i == 1) + sum(delta_e_i == 0) + sum(w_e_i)
      doseDT$pi_e_hat[idx] <- doseDT$y[idx] / doseDT$ESS_e[idx]
      
      # TITE imputation for quasi-events
      tox_exp <- delta_t_i
      eff_exp <- delta_e_i
      
      pending_t <- which(delta_t_i == -1)
      if (length(pending_t) > 0) {
        wt_p <- w_t_i[pending_t]
        pi_t <- doseDT$pi_t_hat[idx]
        tox_exp[pending_t] <- pi_t * (1 - wt_p) / (1 - pi_t * wt_p)
      }
      
      pending_e <- which(delta_e_i == -1)
      if (length(pending_e) > 0) {
        we_p <- w_e_i[pending_e]
        pi_e <- doseDT$pi_e_hat[idx]
        eff_exp[pending_e] <- pi_e * (1 - we_p) / (1 - pi_e * we_p)
      }
      
      # Quasi-events
      doseDT$x_d[idx] <- (u11 * sum(tox_exp * eff_exp) +
                            u10 * sum(tox_exp * (1 - eff_exp)) +
                            u01 * sum((1 - tox_exp) * eff_exp) +
                            u00 * sum((1 - tox_exp) * (1 - eff_exp))) / 100
      
      # PK statistics
      pk_i <- patDT_current$PK[mask]
      doseDT$r_d[idx] <- mean(pk_i)
      doseDT$r_sd[idx] <- ifelse(length(pk_i) > 1, sd(pk_i), pk_i[1] * CV)
    }
  }
  
  return(list(patDT = patDT, doseDT = doseDT, time_next = time_next))
}
