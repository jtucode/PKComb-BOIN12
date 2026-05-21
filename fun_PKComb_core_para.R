###############################################################################
## fun_PKComb_core_para.R
## Parallel wrapper for a single simulation replicate
###############################################################################

fun_PKComb_core_para <- function(index, pkM, toxM, effM, J, K, trueOBDC,
                                 decisionM, pT, qE, pqcorr, csize, cN, N_star,
                                 lambda_d, lambda_e, zeta1, CV, g_P,
                                 ub, u11, u00, current, doselimit,
                                 accrual, susp, tox_win, eff_win,
                                 tox_dist, eff_dist,
                                 use_susp, accrual_random, considerPK) {
  fun_PKCombBOIN12(
    index=index, pkM=pkM, toxM=toxM, effM=effM,
    J=J, K=K, trueOBDC=trueOBDC,
    pT=pT, qE=qE, pqcorr=pqcorr,
    lambda_e=lambda_e, lambda_d=lambda_d, zeta1=zeta1,
    CV=CV, g_P=g_P, csize=csize, cN=cN, N_star=N_star,
    decisionM=decisionM, ub=ub, u11=u11, u00=u00,
    current=current, doselimit=doselimit,
    accrual=accrual, susp=susp,
    tox_win=tox_win, eff_win=eff_win,
    tox_dist=tox_dist, eff_dist=eff_dist,
    use_susp=use_susp, accrual_random=accrual_random,
    considerPK=considerPK
  )
}