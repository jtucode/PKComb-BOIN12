###############################################################################
## fun_combboin12_tite.R
## TITE version of Comb-BOIN12 with time-to-event imputation
## Uses conditional probability TITE imputation (same as PKComb-BOIN12):
##   Pending tox: E[tox_i] = pi_t * (1 - w_t) / (1 - pi_t * w_t)
##   Pending eff: E[eff_i] = pi_e * (1 - w_e) / (1 - pi_e * w_e)
## Suspension rule: wait until susp fraction of current-dose patients complete
###############################################################################

combboin12_tite <- function(rseed, p.true.tox, p.true.eff, dose1, dose2, ntrial,
                            wt, wu, target.tox, target.eff,
                            ncohort, cohortsize, n.earlystop, N1,
                            cutoff.eli.tox, cutoff.eli.eff, extrasafe, offset,
                            accrual.rate, A, maxpen) {
  set.seed(rseed)
  library(Iso)
  library(rjags)
  library(coda)

  startdose <- c(1, 1)
  p.saf <- 0.6 * target.tox
  p.tox <- 1.4 * target.tox
  ndose <- length(p.true.tox)
  npts <- ncohort * cohortsize
  J <- dim(p.true.tox)[1]; K <- dim(p.true.tox)[2]

  lambda1 <- log((1 - p.saf)/(1 - target.tox)) / log(target.tox * (1 - p.saf)/(p.saf * (1 - target.tox)))
  lambda2 <- log((1 - target.tox)/(1 - p.tox)) / log(p.tox * (1 - target.tox)/(target.tox * (1 - p.tox)))
  u.true <- p.true.eff - wt * p.true.tox
  uu <- target.eff - wt * target.tox
  uu <- wu * uu + (1 - wu) * 1

  ## --- Sub-functions (same as original) ---
  select.mtd.comb <- function(target.tox, npts, ntox, cutoff.eli.tox, extrasafe, offset) {
    y <- ntox; n <- npts
    if (nrow(n) > ncol(n) | nrow(y) > ncol(y)) {
      cat("Error: npts and ntox should be rotated.\n"); return()
    }
    elimi <- matrix(0, dim(n)[1], dim(n)[2])
    if (extrasafe) {
      if (n[1,1] >= 3) {
        if (1 - pbeta(target.tox, y[1,1]+1, n[1,1]-y[1,1]+1) > cutoff.eli.tox - offset) {
          elimi[,] <- 1
        }
      }
    }
    for (i in 1:dim(n)[1]) {
      for (j in 1:dim(n)[2]) {
        if (n[i,j] >= 3) {
          if (1 - pbeta(target.tox, y[i,j]+1, n[i,j]-y[i,j]+1) > cutoff.eli.tox) {
            for (id in min(i, dim(n)[1]):dim(n)[1]) {
              for (jd in min(j, dim(n)[2]):dim(n)[2]) { elimi[id, jd] <- 1 }
            }
          }
        }
      }
    }
    if (elimi[1] == 1) {
      selectdose <- t(c(99, 99)); selectdoses <- t(c(99, 99))
    } else {
      phat <- (y + 0.05) / (n + 0.1)
      phat <- Iso::biviso(phat, n + 0.1, warn = TRUE)[,]
      phat.out <- phat; phat.out[n == 0] <- NA
      phat[elimi == 1] <- 1.1
      phat <- phat * (n != 0) + (1E-5) * (matrix(rep(1:dim(n)[1], each = dim(n)[2], len = length(n)), dim(n)[1], byrow = T) +
                                            matrix(rep(1:dim(n)[2], each = dim(n)[1], len = length(n)), dim(n)[1]))
      phat[n == 0] <- 10
      selectdose <- which(abs(phat - target.tox) == min(abs(phat - target.tox)), arr.ind = TRUE)
      if (length(selectdose) > 2) selectdose <- selectdose[1,]
      selectdoses <- matrix(99, nrow = 1, ncol = 2)
      selectdoses[1,] <- selectdose
      selectdoses <- matrix(selectdoses[selectdoses[,2] != 99,], ncol = 2)
      colnames(selectdoses) <- c('DoseA', 'DoseB')
    }
    if (selectdoses[1,1] == 99 & selectdoses[1,2] == 99) {
      return(list(target.tox = target.tox, MTD = 99, p_est = matrix(NA, nrow = dim(npts)[1], ncol = dim(npts)[2])))
    } else {
      return(list(target.tox = target.tox, MTD = selectdoses, p_est = round(phat.out, 2)))
    }
  }

  estimation <- function(p.true.tox, y, x, n) {
    J <- dim(p.true.tox)[1]; K <- dim(p.true.tox)[2]
    tried <- which(n != 0, arr.ind = T)
    yT <- y[tried]; yE <- x[tried]; n1 <- n[tried]; N1 <- length(n1)
    dosem1 <- dosem2 <- NULL
    dose11 <- log(dose1) - mean(log(dose1))
    dose22 <- log(dose2) - mean(log(dose2))
    for (i in 1:length(dose1)) { dosem1 <- rbind(dosem1, rep(dose11[i], dim(p.true.tox)[2])) }
    for (i in 1:length(dose2)) { dosem2 <- cbind(dosem2, rep(dose22[i], dim(p.true.tox)[1])) }
    jags.data <- list("yE" = yE, "n" = n1, "N" = N1, "dose" = tried,
                      "dosem1" = dosem1, "dosem2" = dosem2, "J" = J, "K" = K)
    model_text <- '
      model {
        for(i in 1:N){
          yE[i] ~ dbinom(q[i],n[i])
          logit(q[i]) <- gamma0+gamma1*dosem1[dose[i,1],dose[i,2]]+gamma2*dosem2[dose[i,1],dose[i,2]]+gamma3*dosem1[dose[i,1],dose[i,2]]*dosem1[dose[i,1],dose[i,2]]+gamma4*dosem2[dose[i,1],dose[i,2]]*dosem2[dose[i,1],dose[i,2]]
        }
        gamma0 ~ dt(0, 1/sqrt(10), 1)
        gamma1 ~ dt(0, 1/sqrt(2.5), 1)
        gamma2 ~ dt(0, 1/sqrt(2.5), 1)
        gamma3 ~ dt(0, 1/sqrt(2.5), 1)
        gamma4 ~ dt(0, 1/sqrt(2.5), 1)
        for(i in 1:J){
          for(j in 1:K){
            logit(qp[i,j]) <- gamma0+gamma1*dosem1[i,j]+gamma2*dosem2[i,j]+gamma3*dosem1[i,j]*dosem1[i,j]+gamma4*dosem2[i,j]*dosem2[i,j]
          }
        }
      }'
    modstr <- textConnection(model_text)
    jags.init <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = rseed)
    jags.fit <- jags.model(file = modstr, data = jags.data, inits = jags.init,
                           n.adapt = 5000, n.chains = 1, quiet = TRUE)
    output <- coda.samples(jags.fit, variable.names = c("qp"), n.iter = 5000,
                           progress.bar = "none", thin = 10)
    output <- as.matrix(output)
    qp <- matrix(colMeans(output), nrow = dim(p.true.tox)[1])
    return(list(qp = qp))
  }

  ## --- Storage ---
  Y <- array(0, dim = c(J, K, ntrial))
  X <- array(0, dim = c(J, K, ntrial))
  N <- array(0, dim = c(J, K, ntrial))
  dselect <- matrix(0, ntrial, 2)
  durationV <- rep(0, ntrial)

  ################## simulate trials ###################
  for (trial in 1:ntrial) {
    set.seed(trial + rseed)

    y <- matrix(0, J, K)      # confirmed DLTs
    x <- matrix(0, J, K)      # confirmed efficacy responses
    n <- matrix(0, J, K)      # total enrolled patients
    pu <- matrix(1 - (uu + wt)/(1 + wt), J, K)
    earlystop <- 0
    d <- startdose
    elimi <- matrix(0, J, K)
    elimi.tox <- matrix(0, J, K)

    # Patient-level tracking
    pat_dose_j <- pat_dose_k <- integer(0)
    pat_enroll <- pat_tox_time <- pat_eff_time <- numeric(0)
    pat_tox_out <- pat_eff_out <- integer(0)

    time_current <- 0

    for (pp in 1:ncohort) {
      # Generate outcomes for this cohort
      T.out <- as.integer(runif(cohortsize) < p.true.tox[d[1], d[2]])
      E.out <- as.integer(runif(cohortsize) < p.true.eff[d[1], d[2]])

      # Enroll patients with inter-arrival times
      for (ipt in 1:cohortsize) {
        if (ipt == 1) {
          enroll_t <- time_current
        } else {
          enroll_t <- pat_enroll[length(pat_enroll)] + runif(1, 0, 2/accrual.rate)
        }
        # Event times
        tox_t <- ifelse(T.out[ipt] == 1, enroll_t + runif(1, 0, A[1]), enroll_t + A[1] + 0.001)
        eff_t <- ifelse(E.out[ipt] == 1, enroll_t + runif(1, 0, A[2]), enroll_t + A[2] + 0.001)

        pat_dose_j <- c(pat_dose_j, d[1])
        pat_dose_k <- c(pat_dose_k, d[2])
        pat_enroll <- c(pat_enroll, enroll_t)
        pat_tox_time <- c(pat_tox_time, tox_t)
        pat_eff_time <- c(pat_eff_time, eff_t)
        pat_tox_out <- c(pat_tox_out, T.out[ipt])
        pat_eff_out <- c(pat_eff_out, E.out[ipt])
      }

      n[d[1], d[2]] <- n[d[1], d[2]] + cohortsize

      # Determine decision time using suspension rule:
      # Wait until susp fraction of current-dose patients have completed
      cur_mask <- which(pat_dose_j == d[1] & pat_dose_k == d[2])
      n_cur <- length(cur_mask)
      min_complete <- min(floor(n_cur * maxpen + 1), n_cur)
      # Completion time = max(tox_end, eff_end) for each patient
      completion_times <- pmax(
        pmin(pat_tox_time[cur_mask], pat_enroll[cur_mask] + A[1]),
        pmin(pat_eff_time[cur_mask], pat_enroll[cur_mask] + A[2])
      )
      time_current <- sort(completion_times)[min_complete]
      # Also ensure we're past the last enrollment
      time_current <- max(time_current, pat_enroll[length(pat_enroll)] + 0.001)

      if (n[d[1], d[2]] >= n.earlystop) break

      # Compute TITE statistics at time_current for ALL doses
      y_tite <- matrix(0, J, K)
      x_tite <- matrix(0, J, K)
      ess_t <- matrix(0, J, K)
      ess_e <- matrix(0, J, K)

      for (ii in 1:J) for (jj in 1:K) {
        mask <- which(pat_dose_j == ii & pat_dose_k == jj)
        if (length(mask) == 0) next

        for (idx in mask) {
          t_elapsed <- time_current - pat_enroll[idx]

          # Toxicity status
          if (pat_tox_out[idx] == 1 && pat_tox_time[idx] <= time_current) {
            # Confirmed DLT
            y_tite[ii, jj] <- y_tite[ii, jj] + 1
            ess_t[ii, jj] <- ess_t[ii, jj] + 1
          } else if (t_elapsed >= A[1]) {
            # Completed tox window without DLT
            ess_t[ii, jj] <- ess_t[ii, jj] + 1
          } else {
            # Pending - TITE weight
            w_t <- t_elapsed / A[1]
            ess_t[ii, jj] <- ess_t[ii, jj] + w_t
          }

          # Efficacy status
          if (pat_eff_out[idx] == 1 && pat_eff_time[idx] <= time_current) {
            # Confirmed response
            x_tite[ii, jj] <- x_tite[ii, jj] + 1
            ess_e[ii, jj] <- ess_e[ii, jj] + 1
          } else if (t_elapsed >= A[2]) {
            # Completed eff window without response
            ess_e[ii, jj] <- ess_e[ii, jj] + 1
          } else {
            # Pending
            w_e <- t_elapsed / A[2]
            ess_e[ii, jj] <- ess_e[ii, jj] + w_e
          }
        }
      }

      # Also track confirmed-only counts for elimination
      y_confirmed <- matrix(0, J, K)
      n_confirmed_t <- matrix(0, J, K)
      x_confirmed <- matrix(0, J, K)
      n_confirmed_e <- matrix(0, J, K)

      for (ii in 1:J) for (jj in 1:K) {
        mask <- which(pat_dose_j == ii & pat_dose_k == jj)
        if (length(mask) == 0) next
        for (idx in mask) {
          t_elapsed <- time_current - pat_enroll[idx]
          if (pat_tox_out[idx] == 1 && pat_tox_time[idx] <= time_current) {
            y_confirmed[ii, jj] <- y_confirmed[ii, jj] + 1
            n_confirmed_t[ii, jj] <- n_confirmed_t[ii, jj] + 1
          } else if (t_elapsed >= A[1]) {
            n_confirmed_t[ii, jj] <- n_confirmed_t[ii, jj] + 1
          }
          if (pat_eff_out[idx] == 1 && pat_eff_time[idx] <= time_current) {
            x_confirmed[ii, jj] <- x_confirmed[ii, jj] + 1
            n_confirmed_e[ii, jj] <- n_confirmed_e[ii, jj] + 1
          } else if (t_elapsed >= A[2]) {
            n_confirmed_e[ii, jj] <- n_confirmed_e[ii, jj] + 1
          }
        }
      }

      # Toxicity rate using TITE ESS
      phat_cur <- ifelse(ess_t[d[1], d[2]] > 0, y_tite[d[1], d[2]] / ess_t[d[1], d[2]], 0)

      # Safety elimination (use confirmed counts for Bayesian test)
      nc <- n_confirmed_t[d[1], d[2]]
      if (nc >= 3) {
        if ((1 - pbeta(target.tox, 1 + y_confirmed[d[1], d[2]], 1 + nc - y_confirmed[d[1], d[2]])) >= cutoff.eli.tox) {
          for (i in min(d[1], J):J) {
            for (j in min(d[2], K):K) {
              elimi[i, j] <- 1; elimi.tox[i, j] <- 1
            }
          }
          if (d[1] == 1 && d[2] == 1) { d <- c(99, 99); earlystop <- 1; break }
        }

        # Efficacy elimination (use confirmed counts)
        nc_e <- n_confirmed_e[d[1], d[2]]
        if (nc_e >= 3) {
          if (pbeta(target.eff, 1 + x_confirmed[d[1], d[2]], nc_e - x_confirmed[d[1], d[2]] + 1) >= cutoff.eli.eff) {
            elimi[d[1], d[2]] <- 1
          }
        }

        # Extra safe rule
        if (extrasafe) {
          if (d[1] == 1 && d[2] == 1 && n_confirmed_t[1,1] >= 3) {
            if (1 - pbeta(target.tox, y_confirmed[1,1] + 1, n_confirmed_t[1,1] - y_confirmed[1,1] + 1) > cutoff.eli.tox - offset) {
              d <- c(99, 99); earlystop <- 1; break
            }
          }
        }
      }

      if (sum(elimi == 1) == ndose) { d <- c(99, 99); earlystop <- 1; break }

      # Compute utility using TITE-imputed counts
      # Impute expected outcomes for pending patients at current dose
      pi_t_hat <- ifelse(ess_t[d[1], d[2]] > 0, y_tite[d[1], d[2]] / ess_t[d[1], d[2]], target.tox/2)
      pi_e_hat <- ifelse(ess_e[d[1], d[2]] > 0, x_tite[d[1], d[2]] / ess_e[d[1], d[2]], target.eff)

      # Compute imputed utility for current dose
      # Use imputed x, y based on TITE
      mask_cur <- which(pat_dose_j == d[1] & pat_dose_k == d[2])
      tox_imp <- eff_imp <- numeric(length(mask_cur))
      for (idx in seq_along(mask_cur)) {
        pidx <- mask_cur[idx]
        t_elapsed <- time_current - pat_enroll[pidx]
        # Tox imputation
        if (pat_tox_out[pidx] == 1 && pat_tox_time[pidx] <= time_current) {
          tox_imp[idx] <- 1
        } else if (t_elapsed >= A[1]) {
          tox_imp[idx] <- 0
        } else {
          w_t <- t_elapsed / A[1]
          tox_imp[idx] <- pi_t_hat * (1 - w_t) / (1 - pi_t_hat * w_t)
        }
        # Eff imputation
        if (pat_eff_out[pidx] == 1 && pat_eff_time[pidx] <= time_current) {
          eff_imp[idx] <- 1
        } else if (t_elapsed >= A[2]) {
          eff_imp[idx] <- 0
        } else {
          w_e <- t_elapsed / A[2]
          eff_imp[idx] <- pi_e_hat * (1 - w_e) / (1 - pi_e_hat * w_e)
        }
      }

      y_imp_cur <- sum(tox_imp)
      x_imp_cur <- sum(eff_imp)
      n_imp_cur <- length(mask_cur)

      u_curr <- (x_imp_cur + wt * (n_imp_cur - y_imp_cur)) / (1 + wt)
      pu[d[1], d[2]] <- 1 - pbeta((uu + wt)/(1 + wt), 1 + u_curr, n_imp_cur - u_curr + 1)
      pu <- pu * (1 - elimi)

      if (n[d[1], d[2]] >= N1) { safe <- 1 } else { safe <- 0 }

      ## Dose escalation/de-escalation (using TITE phat)
      if (phat_cur <= lambda1) {
        elevel <- matrix(c(1,0,0,1,-1,0,0,-1,0,0), 2)
        pr_H0 <- rep(0, ncol(elevel))
        for (i in 1:ncol(elevel)) {
          di <- d + elevel[,i]
          if (di[1] >= 1 && di[1] <= J && di[2] >= 1 && di[2] <= K) {
            if (elimi[di[1], di[2]] == 0) {
              pr_H0[i] <- pu[di[1], di[2]] + 1e-6 * (elevel[1,i] + elevel[2,i])
            }
          }
        }
        if (max(pr_H0) == 0) { d <- d } else {
          k <- which(pr_H0 == max(pr_H0))[as.integer(runif(1) * length(which(pr_H0 == max(pr_H0))) + 1)]
          d <- d + elevel[, k]
        }

      } else if (phat_cur > lambda2) {
        delevel <- matrix(c(-1,0,0,-1), 2)
        pr_H0 <- rep(0, ncol(delevel))
        for (i in 1:ncol(delevel)) {
          di <- d + delevel[,i]
          if (di[1] >= 1 && di[2] >= 1) {
            if (elimi[di[1], di[2]] == 0) {
              pr_H0[i] <- pu[di[1], di[2]]
            }
          }
        }
        if (max(pr_H0) == 0) { d <- d } else {
          k <- which(pr_H0 == max(pr_H0))[as.integer(runif(1) * length(which(pr_H0 == max(pr_H0))) + 1)]
          d <- d + delevel[, k]
        }

      } else {
        elevel <- matrix(c(-1,0,0,-1,0,0), 2)
        if (safe == 0) { elevel <- matrix(c(1,0,0,1,-1,0,0,-1,0,0), 2) }
        pr_H0 <- rep(0, ncol(elevel))
        for (i in 1:ncol(elevel)) {
          di <- d + elevel[,i]
          if (di[1] >= 1 && di[1] <= J && di[2] >= 1 && di[2] <= K) {
            if (elimi[di[1], di[2]] == 0) {
              pr_H0[i] <- pu[di[1], di[2]] + 1e-6 * (elevel[1,i] + elevel[2,i])
            }
          }
        }
        if (max(pr_H0) == 0) { d <- d } else {
          k <- which(pr_H0 == max(pr_H0))[as.integer(runif(1) * length(which(pr_H0 == max(pr_H0))) + 1)]
          d <- d + elevel[, k]
        }
      }

      # If selected dose is eliminated, find nearest non-eliminated
      if (d[1] >= 1 && d[1] <= J && d[2] >= 1 && d[2] <= K) {
        if (elimi[d[1], d[2]] == 1) {
          candi <- which(elimi == 0, arr.ind = T)
          if (length(candi) == 0) { earlystop <- 1; break }
          if (length(candi) == 2) { d <- as.vector(candi) }
          if (length(candi) > 2) {
            lowest <- which(apply(candi, 1, sum) == min(apply(candi, 1, sum)))
            if (length(lowest) == 1) { d <- as.vector(candi[lowest,]) } else {
              lowe <- sample(lowest, 1); d <- as.vector(candi[lowe,]) }
          }
        }
      }
    } # end cohort loop

    # Final fully-observed counts
    y_final <- matrix(0, J, K)
    x_final <- matrix(0, J, K)
    n_final <- matrix(0, J, K)
    for (idx in seq_along(pat_dose_j)) {
      ii <- pat_dose_j[idx]; jj <- pat_dose_k[idx]
      n_final[ii, jj] <- n_final[ii, jj] + 1
      y_final[ii, jj] <- y_final[ii, jj] + pat_tox_out[idx]
      x_final[ii, jj] <- x_final[ii, jj] + pat_eff_out[idx]
    }

    Y[,,trial] <- y_final
    X[,,trial] <- x_final
    N[,,trial] <- n_final

    # Duration: time when last patient's outcomes are fully observed
    if (length(pat_enroll) > 0) {
      last_complete <- max(pmax(pat_enroll + A[1], pat_enroll + A[2]))
      durationV[trial] <- last_complete
    }

    # OBDC selection (same as original - uses fully observed data)
    if (earlystop == 1) {
      dselect[trial,] <- c(99, 99)
    } else {
      selcomb <- select.mtd.comb(target.tox, n_final, y_final, cutoff.eli.tox, extrasafe, offset)
      phat.est <- selcomb$p_est
      modelEst <- estimation(p.true.tox, y_final, x_final, n_final)
      modelEstE <- modelEst$qp
      if (is.na(phat.est[1])) {
        dselect[trial,] <- c(99, 99)
      } else {
        selcomb <- selcomb$MTD
        u.mean <- modelEstE - wt * phat.est
        u.mean[n_final == 0] <- -100
        u.mean[elimi == 1] <- -100
        phat.est[n_final == 0] <- 100
        for (i in min(selcomb[1], J):J) {
          for (j in min(selcomb[2], K):K) {
            if (i != selcomb[1] | j != selcomb[2]) {
              if ((phat.est[i,j] - phat.est[selcomb[1], selcomb[2]]) >= 0.05) {
                u.mean[i, j] <- -100
              }
            }
          }
        }
        dopt <- which(u.mean == max(u.mean), arr.ind = TRUE)
        dselect[trial, 1] <- dopt[1]
        dselect[trial, 2] <- dopt[2]
      }
    }
  } # end trial loop

  # Output (same format as original)
  sel <- matrix(0, J, K)
  pts <- apply(N, c(1, 2), mean)
  for (i in 1:J) for (j in 1:K) {
    sel[i, j] <- sum(dselect[,1] == i & dselect[,2] == j) / ntrial * 100
  }
  earlystop_pct <- sum(dselect[,1] == 99) / ntrial * 100

  sel.overdose <- sum((p.true.tox > target.tox) * sel)
  pts.overdose <- sum((p.true.tox > target.tox) * pts)
  ppts.overdose <- sum((p.true.tox > target.tox) * pts) * 100 / sum(pts)

  candi.OBD <- which((p.true.eff >= target.eff) * (p.true.tox <= target.tox) == 1)
  sel.OBD <- pts.OBD <- ppts.OBD <- 0
  if (length(candi.OBD) > 0) {
    true.OBD <- candi.OBD[which(u.true[candi.OBD] == max(u.true[candi.OBD]))]
    sel.OBD <- sum(sel[true.OBD])
    pts.OBD <- sum(pts[true.OBD])
    ppts.OBD <- sum(pts[true.OBD]) * 100 / sum(pts)
  }

  results <- list(
    p.true.tox = p.true.tox, p.true.eff = p.true.eff, u.true = u.true,
    sel = sel, pts = pts, earlystop = earlystop_pct,
    present.vec = c(sel.OBD, pts.OBD, ppts.OBD, 0, 0, 0,
                    sel.overdose, pts.overdose, ppts.overdose, round(mean(durationV), 2))
  )
  return(results)
}
