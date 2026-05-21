
###############################################################################
## fun_PKCombBOIN12_v2.R
## Simulate one complete PKComb-BOIN12 trial
## v2: Adds elim_reason tracking (SAFETY / EFFICACY / PK / NA)
###############################################################################

fun_PKCombBOIN12 <- function(index, pkM, toxM, effM, J, K, trueOBDC,
                             pT, qE, pqcorr, lambda_e, lambda_d, zeta1,
                             CV, g_P, csize, cN, N_star, decisionM,
                             ub, u11 = 60, u00 = 40,
                             current = c(1, 1), doselimit = Inf,
                             accrual, susp, tox_win, eff_win,
                             tox_dist, eff_dist,
                             use_susp = TRUE, accrual_random = FALSE,
                             considerPK = TRUE) {
  set.seed(index)
  dN <- J * K
  time_current <- 0

  ## doseDT: row-major ordering (j varies first within each k)
  doseDT <- expand.grid(j = 1:J, k = 1:K)
  doseDT$id <- 1:nrow(doseDT)
  doseDT$PK  <- pkM[cbind(doseDT$j, doseDT$k)]
  doseDT$tox <- toxM[cbind(doseDT$j, doseDT$k)]
  doseDT$eff <- effM[cbind(doseDT$j, doseDT$k)]
  doseDT$n <- 0; doseDT$x <- 0; doseDT$y <- 0; doseDT$keep <- 1
  doseDT$ESS_t <- -1; doseDT$ESS_e <- -1
  doseDT$pi_t_hat <- pT/2; doseDT$pi_e_hat <- qE
  doseDT$x_d <- 0; doseDT$r_d <- 0; doseDT$r_sd <- 0
  doseDT$elim_reason <- NA_character_

  total_pts <- csize * cN
  patDT <- data.frame(
    cid = rep(1:cN, each = csize), pid = rep(1:csize, cN), id = 1:total_pts,
    PK = NA, tox_prob = NA, eff_prob = NA,
    enroll = NA, toxend = NA, effend = NA, toxobs = -1, effobs = -1,
    toxdat = 9999, effdat = 9999, tox_confirm = NA, eff_confirm = NA,
    dj = NA, dk = NA
  )

  earlystop <- FALSE
  record_j <- record_k <- rep(-1, cN)
  current_jk <- current

  for (i in 1:cN) {
    upd <- fun_PKComb_update(
      cid=i, CV=CV, g_P=g_P, patDT=patDT, doseDT=doseDT,
      current_jk=current_jk, time_current=time_current,
      tox_win=tox_win, eff_win=eff_win, csize=csize,
      tox_dist=tox_dist, eff_dist=eff_dist,
      accrual=accrual, susp=susp,
      u11=u11, u10=0, u01=100, u00=u00,
      use_susp=use_susp, accrual_random=accrual_random, pqcorr=pqcorr
    )
    patDT <- upd$patDT; doseDT <- upd$doseDT; time_current <- upd$time_next

    decision <- PKCombBOIN12_one(
      doseDT=doseDT, current_jk=current_jk,
      pT=pT, qE=qE, J=J, K=K,
      lambda_e=lambda_e, lambda_d=lambda_d, zeta1=zeta1,
      N_star=N_star, csize=csize, decisionM=decisionM,
      ub=ub, u11=u11, u00=u00
    )
    doseDT <- decision$doseDT
    record_j[i] <- current_jk[1]; record_k[i] <- current_jk[2]
    current_jk <- decision$newdose

    if (is.null(current_jk) || any(is.na(current_jk))) { earlystop <- TRUE; break }
    idx_nxt <- which(doseDT$j == current_jk[1] & doseDT$k == current_jk[2])
    if (length(idx_nxt) > 0 && doseDT$n[idx_nxt] + csize > doselimit) break
  }

  ## Final data with fully observed outcomes
  doseDT_final <- doseDT
  patDT_obs <- patDT[!is.na(patDT$dj), ]
  for (idx in 1:nrow(doseDT_final)) {
    jj <- doseDT_final$j[idx]; kk <- doseDT_final$k[idx]
    mask <- which(patDT_obs$dj == jj & patDT_obs$dk == kk)
    if (length(mask) > 0) {
      doseDT_final$x[idx] <- sum(patDT_obs$toxobs[mask] == 1)
      doseDT_final$y[idx] <- sum(patDT_obs$effobs[mask] == 1)
      doseDT_final$r_d[idx] <- mean(patDT_obs$PK[mask])
    }
  }

  if (!earlystop) {
    OBDC <- fun_PKCombBOIN12_OBDC(doseDT_final, J, K, pT, u11, u00, zeta1)
    rN <- sum(doseDT_final$keep == 1 & doseDT_final$n > 0)
  } else {
    OBDC <- c(99, 99); rN <- 0
  }

  ## Metrics
  ## trueOBDC is a matrix (n_obdc x 2); handle multiple tied OBDCs
  if (is.null(dim(trueOBDC))) trueOBDC <- matrix(trueOBDC, nrow = 1)
  no_obdc <- all(trueOBDC[1, ] < 0)

  select_OBDC <- if (no_obdc) {
    as.integer(all(OBDC == c(99, 99)))
  } else {
    as.integer(any(apply(trueOBDC, 1, function(row) all(OBDC == row))))
  }

  num_at_OBDC <- if (!no_obdc) {
    total <- 0
    for (r in 1:nrow(trueOBDC)) {
      tid <- which(doseDT_final$j == trueOBDC[r, 1] & doseDT_final$k == trueOBDC[r, 2])
      if (length(tid) > 0) total <- total + doseDT_final$n[tid]
    }
    total
  } else NA

  overdose_mask <- doseDT_final$tox > pT
  num_overdose <- ifelse(any(overdose_mask), sum(doseDT_final$n[overdose_mask]), 0)
  select_overdose <- as.integer(
    !all(OBDC == c(99, 99)) &&
    any(overdose_mask & doseDT_final$j == OBDC[1] & doseDT_final$k == OBDC[2])
  )

  ## Duration
  patDT_final <- patDT[!is.na(patDT$tox_confirm), ]
  duration <- if (nrow(patDT_final) > 0) {
    max(c(patDT_final$toxend, patDT_final$effend), na.rm = TRUE)
  } else NA

  ## Per-dose n
  n_vec <- numeric(J * K)
  for (idx in 1:nrow(doseDT_final)) {
    pos <- (doseDT_final$k[idx] - 1) * J + doseDT_final$j[idx]
    n_vec[pos] <- doseDT_final$n[idx]
  }
  names(n_vec) <- paste0("n_", rep(1:J, K), "_", rep(1:K, each = J))

  ## Per-dose elimination reason
  elim_vec <- character(J * K)
  for (idx in 1:nrow(doseDT_final)) {
    pos <- (doseDT_final$k[idx] - 1) * J + doseDT_final$j[idx]
    elim_vec[pos] <- ifelse(is.na(doseDT_final$elim_reason[idx]), "NONE",
                            doseDT_final$elim_reason[idx])
  }
  names(elim_vec) <- paste0("elim_", rep(1:J, K), "_", rep(1:K, each = J))

  return(data.frame(
    earlystop = earlystop,
    OBDC_j = OBDC[1], OBDC_k = OBDC[2],
    rN = rN,
    trueOBDC_j = trueOBDC[1, 1], trueOBDC_k = trueOBDC[1, 2],
    select_OBDC = select_OBDC,
    select_overdose = select_overdose,
    num_at_OBDC = num_at_OBDC,
    num_overdose = num_overdose,
    duration = duration,
    t(n_vec),
    t(elim_vec)
  ))
}
