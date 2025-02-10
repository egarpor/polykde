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
Meilán-Vila, 2024).

## Replicability

The folder `/replication` contains the scripts to replicate the
simulations and numerical experiments of the Supplementary Material
(SM):

-   The script `kde_sims.R` reproduces the asymptotic normality
    experiment in Section B.1 of the SM.
-   The script `kde_efic.R` computes the kernel efficiency table in
    Section B.2 of the SM.
-   The two scripts `jsd-sims-k2-S2.R` and `jsd-sims-k3-S10^2.R`
    reproduce two simulation experiments for the *k*-sample test in
    Section B.3 of the SM.

## References

García-Portugués, E. and Meilán-Vila, A. (2024). Kernel density
estimation with polyspherical data and its applications.
*arXiv:2411.04166*.
[doi:10.48550/arXiv.2411.04166](https://doi.org/10.48550/arXiv.2411.04166).

García-Portugués, E. and Meilán-Vila, A. (2023). Hippocampus shape
analysis via skeletal models and kernel smoothing. In Larriba, Y. (Ed.),
*Statistical Methods at the Forefront of Biomedical Advances*,
pp. 63–82. Springer, Cham.
[doi:10.1007/978-3-031-32729-2_4](https://doi.org/10.1007/978-3-031-32729-2_4).
