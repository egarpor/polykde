# polykde

[![License:
GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R build
status](https://github.com/egarpor/polykde/workflows/R-CMD-check/badge.svg)](https://github.com/egarpor/polykde/actions)
[![R build
status](https://github.com/egarpor/polykde/workflows/test-coverage/badge.svg)](https://github.com/egarpor/polykde/actions)
[![](https://codecov.io/gh/egarpor/polykde/branch/main/graph/badge.svg)](https://app.codecov.io/gh/egarpor/polykde)
[![](https://www.r-pkg.org/badges/version/polykde?color=green)](https://cran.r-project.org/package=polykde)
[![](http://cranlogs.r-pkg.org/badges/grand-total/polykde)](https://cran.r-project.org/package=polykde)
[![](http://cranlogs.r-pkg.org/badges/last-month/polykde)](https://cran.r-project.org/package=polykde)

## Overview

Companion package for the article *Kernel density estimation with
polyspherical data and its applications* (García-Portugués and
Meilán-Vila, 2025).

## Replicability

The folder `/replication` contains the scripts to replicate the
numerical experiments and real data application of the paper and its
Supplementary Material (SM):

- The script `kde-sims.R` reproduces the asymptotic normality experiment
  (Figures 5–8 in the SM).
- The script `kde-effic.R` computes the kernel efficiency table (Table 1
  in the SM) and the kernel and kernel efficiency graphs (Figure 1 in
  the paper).
- The scripts `jsd-sims-k2-S2.R`, `jsd-sims-hippo.R`, and
  `jsd-sims-k3-S10^2.R` reproduce two simulation experiments for the
  $`k`$-sample test in (Figures 9–12 in the SM).
- The scripts `kde-spoke-dirs.R` and `test-spoke-dirs.R` reproduce the
  real data application on the hippocampus shape analysis (Figure 3 in
  the paper and Figure 13 in the SM, and Figure 4 in the paper,
  respectively).

## References

García-Portugués, E. and Meilán-Vila, A. (2025). Kernel density
estimation with polyspherical data and its applications. *Journal of the
American Statistical Association*, to appear.
[doi:10.1080/01621459.2025.2521898](https://doi.org/10.1080/01621459.2025.2521898).

García-Portugués, E. and Meilán-Vila, A. (2023). Hippocampus shape
analysis via skeletal models and kernel smoothing. In Larriba, Y. (Ed.),
*Statistical Methods at the Forefront of Biomedical Advances*,
pp. 63–82. Springer, Cham.
[doi:10.1007/978-3-031-32729-2_4](https://doi.org/10.1007/978-3-031-32729-2_4).
