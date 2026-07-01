###############################################################################
## fun_PKComb_fixsimu_v2.R
## Main simulation driver for PKComb-BOIN12
## v2: Returns three J×K elimination-reason matrices (safety/efficacy/PK %)
###############################################################################

fun_PKComb_fixsimu <- function(J, K, pkV, toxV, effV, pT, qE,
                               pqcorr = 0, psi0PK=6000, zeta1, CV = 0.25, g_P = 1,
                               csize = 3, cN = 17, N_star = 6,
                               utility = TRUE, u11 = 60, u00 = 40,
                               cutoff_tox = 0.95, cutoff_eff = 0.9,cutoff_pk=0.95,
                               repsize = 5000, n_cores = 1,
                               accrual = 10, susp = 0.5,
                               tox_win = 30, eff_win = 60,
                               tox_dist = "Uniform", eff_dist = "Uniform",
                               use_susp = TRUE, accrual_random = FALSE,
                               considerPK = TRUE) {
  # require(parallel)
  # require(dplyr)

  toxM <- matrix(toxV, nrow = J, ncol = K, byrow = TRUE)
  effM <- matrix(effV, nrow = J, ncol = K, byrow = TRUE)
  pkM  <- matrix(pkV,  nrow = J, ncol = K, byrow = TRUE)

  trueOBDC <- findOBDC_RDS(toxM, effM, pkM, pT, qE, zeta1, u11, u00)
  true_admissible<-fun_true_admissible(J,K,toxV,effV,pkV,pT,qE,zeta1)
  True_Utility <- fun_true_utility(J,K,toxV,effV,pkV, u00, u11)

  phi1 <- 0.6 * pT; phi2 <- 1.4 * pT
  lambda_e <- log((1 - phi1)/(1 - pT)) / log(pT*(1 - phi1)/phi1/(1 - pT))
  lambda_d <- log((1 - pT)/(1 - phi2)) / log(phi2*(1 - pT)/pT/(1 - phi2))

  psi1PK <- 0.6 * psi0PK
  #zeta1 <- (psi0PK + psi1PK) / 2

  u <- 100 * qE * (1 - pT) + u00 * (1 - pT) * (1 - qE) + u11 * pT * qE
  ub <- (u + (100 - u) / 2) / 100

  decisionM <- fun_PKCombBOIN12dec(pT, qE, lambda_e, lambda_d, csize, cN, cutoff_tox, cutoff_eff)

  cat(paste0("True OBDC: ",
             paste(apply(trueOBDC, 1, function(r) paste0("(", r[1], ",", r[2], ")")),
                   collapse = " "), "\n"))
  cat(paste0("Lambda_e = ", round(lambda_e, 4), ", Lambda_d = ", round(lambda_d, 4),
             ", Zeta1 = ", round(zeta1, 1), ", ub = ", round(ub, 4), "\n"))

  ResultDF <- mclapply(1:repsize, FUN = function(x)
    fun_PKComb_core_para(
      index=x, pkM=pkM, toxM=toxM, effM=effM,
      J=J, K=K, trueOBDC=trueOBDC, true_admissible=true_admissible, decisionM=decisionM,
      pT=pT, qE=qE, pqcorr=pqcorr, csize=csize, cN=cN, N_star=N_star,
      lambda_d=lambda_d, lambda_e=lambda_e, zeta1=zeta1, cutoff_pk=cutoff_pk,
      CV=CV, g_P=g_P,
      ub=ub, u11=u11, u00=u00, current=c(1,1), doselimit=Inf,
      accrual=accrual, susp=susp, tox_win=tox_win, eff_win=eff_win,
      tox_dist=tox_dist, eff_dist=eff_dist,
      use_susp=use_susp, accrual_random=accrual_random, considerPK=considerPK
    ), mc.cores = n_cores
  ) %>% do.call(rbind, .)

  ## Aggregate
  No_OBDC <- mean(ResultDF$no_obdc, na.rm=T)*100
  Admissiblity <- mean(ResultDF$admissible_score,na.rm=T)*100
  coherence <- mean(ResultDF$coherence,na.rm=T)*100
  sel_pct <- mean(ResultDF$select_OBDC) * 100
  sel_overdose_pct <- mean(ResultDF$select_overdose) * 100
  early_pct <- mean(ResultDF$earlystop) * 100
  avg_dur <- mean(ResultDF$duration, na.rm = TRUE)
  avg_OBDC <- mean(ResultDF$num_at_OBDC, na.rm = TRUE)
  avg_od <- mean(ResultDF$num_overdose, na.rm = TRUE)

  U_mat <- matrix(0, J, K)# mean(ResultDF$Utility_vals,na.rm=T)
  sel_mat <- matrix(0, J, K)
  pts_mat <- matrix(0, J, K)
  pct_mat <- matrix(0, J, K)
  elim_safety_mat   <- matrix(0, J, K)
  elim_efficacy_mat <- matrix(0, J, K)
  elim_pk_mat       <- matrix(0, J, K)

  total_planned_n <- csize * cN

  for (jj in 1:J) for (kk in 1:K) {
    sel_mat[jj, kk] <- mean(ResultDF$OBDC_j == jj & ResultDF$OBDC_k == kk) * 100
    cn <- paste0("n_", jj, "_", kk)
    if (cn %in% names(ResultDF)) {
      pts_mat[jj, kk] <- mean(ResultDF[[cn]])
    }
    en <- paste0("elim_", jj, "_", kk)
    if (en %in% names(ResultDF)) {
      elim_safety_mat[jj, kk]   <- mean(ResultDF[[en]] == "SAFETY")   * 100
      elim_efficacy_mat[jj, kk] <- mean(ResultDF[[en]] == "EFFICACY") * 100
      elim_pk_mat[jj, kk]       <- mean(ResultDF[[en]] == "PK")       * 100
    }
    un <- paste0("u_", jj, "_", kk)
    if (un %in% names(ResultDF)) {
      U_mat[jj, kk] <- mean(ResultDF[[un]],na.rm=T)+100
    }
  }

  ## Compute percentage using actual treated patients as denominator
  avg_total_treated <- sum(pts_mat)
  if (avg_total_treated > 0) {
    pct_mat <- (pts_mat / avg_total_treated) * 100
  }

  ## Percentage of patients at OBDC
  no_obdc <- all(trueOBDC[1, ] < 0)
  avg_total_treated <- sum(pts_mat)
  pct_OBDC <- if (!no_obdc && avg_total_treated > 0) {
    (avg_OBDC / avg_total_treated) * 100
  } else {
    NA
  }

  cat("\n=== PKComb-BOIN12 Results ===\n")
  cat("No OBDC %", round(No_OBDC,1),"\n")
  cat("Admissibility %:", round(Admissiblity, 1), "\n")
  cat("Coherence %:", round(coherence, 1), "\n")
  cat("Selection % of OBDC:", round(sel_pct, 1), "\n")
  cat("Selection % of overdose:", round(sel_overdose_pct, 1), "\n")
  cat("Early termination %:", round(early_pct, 1), "\n")
  cat("Avg patients at OBDC:", round(avg_OBDC, 1), "\n")
  cat("% patients at OBDC:", round(pct_OBDC, 1), "\n")
  cat("Avg patients at overdose:", round(avg_od, 1), "\n")
  cat("Avg duration:", round(avg_dur, 1), "\n")

  cat("\nSelection Probability (%):\n")
  print(round(sel_mat, 1))

  cat("\nAvg Number of Patients:\n")
  print(round(pts_mat, 1))

  cat("\nPercentage of Patients (%):\n")
  print(round(pct_mat, 1))

  cat("\nElimination for SAFETY (%):\n")
  print(round(elim_safety_mat, 1))

  cat("\nElimination for EFFICACY (%):\n")
  print(round(elim_efficacy_mat, 1))

  cat("\nElimination for PK (%):\n")
  print(round(elim_pk_mat, 1))

  return(list(
    NOOBDC = No_OBDC, #No OBDC %
    trueOBDC=trueOBDC, #TRUE OBDC
    trueAdmissible =true_admissible, #TRUE ADMISSIBILITY SET OF DOSE COMBINATIONS
    trueUtility=True_Utility, #The RDS score matrix, true values
    utility_score_mat =U_mat, # Avg. utility score for OBDC selection
    sel_mat=sel_mat, # Selection probability of dose combination
    pts_mat=pts_mat, # Avg. no. of patients allocated to the dose combinations
    pct_mat=pct_mat, #Avg. % age of patients allocated to the dose combinations
    sel_OBDC=sel_pct, # Selection % of correct OBDC
    coherence=coherence, #Coherence metric average of ( escalate when All DLT in cohort, de-escalate when no DLT in cohort)/ (total cohorts with no DLT or all DLT)
    admissibility=Admissiblity, # Admissibility : %age of times, atleast one dose combination was admissible in the true admissible set, when OBDC is present.
    sel_overdose=sel_overdose_pct, #%age of overdose combinations were selected
    early_stop=early_pct, #early stop %age
    avg_duration=avg_dur, #Average trial duration
    avg_num_OBDC=avg_OBDC, # Avg. no of patients at OBDC
    pct_OBDC=pct_OBDC, # percentage of patients at OBDC
    avg_overdose=avg_od, # avg. % of overdose patients
    elim_safety_mat=elim_safety_mat, #Elimination for SAFETY (%)
    elim_efficacy_mat=elim_efficacy_mat, #Elimination for EFFICACY (%)
    elim_pk_mat=elim_pk_mat, #Elimination for PK (%)
    raw=ResultDF
  ))
}
