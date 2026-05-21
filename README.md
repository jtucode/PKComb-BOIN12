# PKComb-BOIN12

**PKComb-BOIN12: A Pharmacokinetics-Informed Bayesian Optimal Interval Design for Dose Optimization in Cancer Drug-Combination Trials**

PKComb-BOIN12 extends [PKBOIN-12](https://github.com/EugeneHao/PKBOIN-12) (Sun and Tu, 2024) to the two-agent combination setting via the [Comb-BOIN12](https://doi.org/10.1080/19466315.2024.2370403) framework (Lu et al., 2025). The design integrates pharmacokinetic (PK) information with toxicity and efficacy outcomes to identify the optimal biological dose combination (OBDC) in early-phase oncology trials.

## Features

- **PK-informed dose escalation**: PK data guides the candidate set in the safe zone, directing escalation toward combinations with adequate drug exposure
- **PK-informed OBDC selection**: Dose combinations with insufficient PK exposure are excluded from the final admissible set
- **PK elimination rule**: Dose combinations with strong evidence of sub-therapeutic exposure are removed during the trial
- **TITE extension**: TITE-PKComb-BOIN12 handles late-onset toxicity and efficacy outcomes, reducing trial duration by 33%–47%
- **Model-assisted**: Uses the rank-based desirability score (RDS) from BOIN12 — no real-time model fitting required

## File Description

| File | Description |
|------|-------------|
| `fun_PKCombBOIN12.R` | Trial-level orchestrator for the PKComb-BOIN12 design |
| `fun_PKCombBOIN12_one.R` | Single-step dose decision logic (4-case flowchart + elimination rules) |
| `fun_PKCombBOIN12_OBDC.R` | OBDC selection at trial completion (posterior mean utility + PK restriction) |
| `fun_PKCombBOIN12dec.R` | BOIN escalation/de-escalation boundary computation |
| `fun_PKComb_update.R` | Cohort update: generates PK-linked outcomes, computes TITE statistics |
| `fun_PKComb_fixsimu.R` | Simulation driver for fixed (non-TITE) scenarios |
| `fun_PKComb_core_para.R` | Parallel simulation wrapper |
| `fun_findOBDC.R` | True OBDC identification from scenario parameters |
| `fun_pava_comb.R` | Bivariate isotonic regression via pool-adjacent-violators algorithm (PAVA) |

## Quick Start

```r
# Source all functions
source("fun_PKCombBOIN12.R")
source("fun_PKCombBOIN12_one.R")
source("fun_PKCombBOIN12_OBDC.R")
source("fun_PKCombBOIN12dec.R")
source("fun_PKComb_update.R")
source("fun_PKComb_fixsimu.R")
source("fun_PKComb_core_para.R")
source("fun_findOBDC.R")
source("fun_pava_comb.R")

# Define a scenario (e.g., 4x4 dose matrix)
J <- 4; K <- 4
toxV <- c(0.05, 0.10, 0.15, 0.20,
          0.10, 0.15, 0.20, 0.25,
          0.15, 0.20, 0.25, 0.30,
          0.20, 0.25, 0.30, 0.35)
effV <- c(0.10, 0.20, 0.30, 0.40,
          0.20, 0.30, 0.40, 0.50,
          0.30, 0.40, 0.55, 0.55,
          0.40, 0.50, 0.55, 0.55)
pkV  <- c(2000, 3000, 4000, 5000,
          3000, 4000, 5500, 6500,
          4000, 5500, 7000, 8000,
          5000, 6500, 8000, 9000)

# Run simulation
results <- fun_PKComb_fixsimu(
  J = J, K = K,
  toxV = toxV, effV = effV, pkV = pkV,
  target_tox = 0.35, target_eff = 0.20,
  r_P = 6000, CV = 0.25, g_P = 1,
  ncohort = 17, cohortsize = 3,
  ntrial = 2000
)
```

## Design Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `r_P` | Target PK value | Trial-specific |
| `zeta_1` | PK-sufficiency cutoff ($0.8 \times r_P$) | Derived |
| `CV` | Coefficient of variation for individual PK | 0.25 |
| `g_P` | PK-outcome coupling ratio | 1.0 |
| `C_T` | Safety elimination cutoff | 0.95 |
| `C_E` | Efficacy elimination cutoff | 0.90 |
| `C_P` | PK elimination cutoff | 0.95 |
| `N*` | Per-cell sample size cutoff | 6 |
| `lambda_1`, `lambda_2` | BOIN escalation/de-escalation boundaries | Derived from `p_T` |

## Dependencies

- R (≥ 4.0)
- `truncnorm` — truncated normal distribution
- `arm` — Bayesian generalized linear models (`bayesglm`)
- `parallel` — parallel simulation (optional)

## References

- **PKComb-BOIN12**: Tu J, Lin R, Mukherjee A (2026). PKComb-BOIN12: A Pharmacokinetics-Informed Bayesian Optimal Interval Design for Dose Optimization in Cancer Drug-Combination Trials. *Manuscript in preparation*.
- **PKBOIN-12**: Sun H, Tu J (2024). PKBOIN-12: A Bayesian Optimal Interval Phase I/II Design Incorporating Pharmacokinetics Outcomes to Find the Optimal Biological Dose. *Pharmaceutical Statistics*. [DOI: 10.1002/pst.2444](https://doi.org/10.1002/pst.2444)
- **Comb-BOIN12**: Lu M, Zhang J, Yuan Y, Lin R (2025). Comb-BOIN12: A Utility-Based Bayesian Optimal Interval Design for Dose Optimization in Cancer Drug-Combination Trials. *Statistics in Biopharmaceutical Research*, 17(2), 266–276. [DOI: 10.1080/19466315.2024.2370403](https://doi.org/10.1080/19466315.2024.2370403)
- **BOIN12**: Lin R, Zhou Y, Yan F, Li D, Yuan Y (2020). BOIN12: Bayesian Optimal Interval Phase I/II Trial Design for Utility-Based Dose Finding in Immunotherapy and Targeted Therapies. *JCO Precision Oncology*, 4, 1393–1402.

## License

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License](https://creativecommons.org/licenses/by-nc/4.0/).

For commercial use, please contact the author: Jieqi Tu (jieqi.tu@lilly.com).
