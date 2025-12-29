# rocvb

## Overview

`rocvb` provides point estimates and confidence intervals for ROC-based diagnostic 
accuracy measures when disease verification is subject to **verification bias**.

`rocvb` currently supports inference for the following accuracy metrics for two-class 
continuous tests under missing-at-random (MAR) disease verification:

- Area under the ROC curve (AUC)
- Sensitivity at a fixed level of specificity
- Maximum Youden index

Multiple types of confidence intervals are provided, including bootstrap-based,
MOVER-based, and empirical likelihood–based intervals.

For each estimand, results are returned **simultaneously** using four bias-corrected 
estimators for sensitivity and specificity proposed by Alonzo and Pepe (2005):

- **FI**: Full imputation  
- **MSI**: Mean score imputation  
- **IPW**: Inverse probability weighting  
- **SPE**: Semiparametric efficient estimation  

---

## Installation

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("swang1021/rocvb")
```

---

## Technical details

For each target, the package returns:
- point estimates based on FI, MSI, IPW, and SPE
- confidence intervals based on bootstrap and/or hybrid empirical likelihood methods
- sample size and missingness summaries

The methods implemented in `rocvb` assume that disease verification is
**missing at random (MAR)** conditional on the observed test result and covariate.

Let  
- \(T_i\) denote the continuous test measurement for subject \(i\),  
- \(A_i\) denote a covariate,  
- \(D_i \in \{0,1\}\) denote the true disease status, and  
- \(V_i \in \{0,1\}\) indicate whether the disease status is verified.

### Disease model

The disease process is modeled using a **probit regression** based on verified
subjects only. The conditional probability of disease is given by

\[
\rho_i = P(D_i = 1 \mid T_i, A_i)
       = \Phi(\alpha + \beta T_i + \gamma A_i),
\]

where \(\Phi(\cdot)\) denotes the standard normal cumulative distribution function.
This model is used to impute disease probabilities for unverified subjects and
forms the basis for the FI, MSI, and SPE estimators.

### Verification model

The verification model is estimated using a logit regression model linear in
\(T\) and \(A\) based on all subjects, given by

\[
\operatorname{logit}(\pi_i)
= \log\!\left( \frac{\pi_i}{1 - \pi_i} \right)
= \alpha + \beta T_i + \gamma A_i,
\quad i = 1, \ldots, n.
\]

Here, \(\pi_i = P(V_i = 1 \mid T_i, A_i)\) denotes the probability that the disease
status of subject \(i\) is verified.

This model characterizes the verification process under the MAR assumption and
is used to construct inverse probability weights for the IPW and SPE estimators.

---

## Usage

```r
library(rocvb)

set.seed(123)
n <- 100
T <- abs(rnorm(n))
A <- abs(rnorm(n))

score <- 0.3 * T + 0.3 * A + rnorm(n, sd = 1)
D <- as.logical(score > stats::quantile(score, 0.7))
D[sample(n, 30)] <- NA

auc.ci.mar(T, D, A, n.boot = 50, plot = FALSE)
sen.ci.mar(T, D, A, p = 0.8, n.boot = 50, plot = FALSE)
yi.ci.mar(T, D, A, n.boot = 50, plot = FALSE)

#For more details, see
?auc.ci.mar
?sen.ci.mar
?yi.ci.mar
```

---

## References

Alonzo, T. A. and Pepe, M. S. (2005). Assessing accuracy of a continuous
screening test in the presence of verification bias. *Journal of the Royal
Statistical Society: Series C (Applied Statistics)*.

Wang, S., Shi, S., and Qin, G. (2025). Interval estimation for the Youden index
of a continuous diagnostic test with verification biased data. *Statistical Methods in Medical Research*.

Wang, S., Shi, S., and Qin, G. (2025). Empirical likelihood inference for the area
under the ROC curve with verification-biased data. Manuscript under peer review.

Wang, S., Shi, S., and Qin, G. (2025). Empirical likelihood-based confidence
intervals for sensitivity of a continuous test at a fixed level of specificity
with verification bias. Manuscript under peer review.
