# polykde 1.0.0

* Initial release.

# polykde 1.1.3

* Add cross-validation bandwidths with exact loss for vMF distributions in `bw_cv_polysph()`.
* Add `fast_log_c_vMF()` as a faster version of `rotasym::c_vMF()`.
* Add exact MISE for vMF distributions.
* Fix error in r-devel-linux-x86_64-fedora-clang reported on 2025-04-15.

# polykde 1.1.4

* Add hippocampus dataset to give replicability to García-Portugués and Meilán-Vila (2024) (doi:10.48550/arXiv.2411.04166).

# polykde 1.1.7

* Add `d_mvmf_polysph()` and `r_mvmf_polysph()` to handle mixtures of vMF distributions on polyspheres.
* Add exact ISE for vMF distributions.
* Expand support of `fast_log_c_vMF()` to large orders.
* Add Hammer projection convenience `hammer_to_sph()` and `fib_latt()`.
* Update reference García-Portugués and Meilán-Vila (2025) to the published version.

# polykde 1.2.0

* Add `arcsinh` to `bw_cv_polysph()` to handle high-dimensional data.
* Move `Bessel` from `Suggests` to `Imports`.
* Fix the Monte Carlo LSCV loss in `bw_cv_polysph(exact_vmf = FALSE)`, which mixed densities with respect to the Lebesgue and uniform measures.
* Extend and speed up the spline approximation of the Bessel functions used by `fast_log_c_vMF()`: the tabulated orders now reach `nu = 24.5`, interpolation is done on `log(x)`, small arguments use an exact series, and non-tabulated orders fall back to `besselI()` instead of erroring.
* Robustify the JSD test in `hom_test_polysph()` against non-finite values; use the standard permutation p-value `(1 + #{T* >= T}) / (B + 1)`.
* Fix the normalizing constants of the Epanechnikov and softplus product kernels in `grad_hess_kde_polysph()`.
* Fix `proj_grad_kde_polysph()` with `sparse = TRUE`, which could request too few eigenvalues for small `r`.
* Accept a length-`r` vector `kappa` for one-component mixtures in `d_mvmf_polysph()`, `r_mvmf_polysph()`, and `mise_vmf_polysph()`.
* Restore the user's random seed on exit whenever `seed_mc`, `seed_psi`, or `seed_jsd` are given.
* Validate inputs in `log_cv_kde_polysph()`, `bw_lcv_min_epa()`, and `proj_polysph()`.
* Add a package logo.
* And the reference García-Portugués and Meilán-Vila (2023) (doi:10.1007/978-3-031-32729-2_4) to CITATION.
* Documentation improvements.

# polykde 1.2.1

* Unrandomize the tests in `tests_bwd.R`, `tests_distr.R`, `tests_grad_hess.R`, `tests_kde.R`, `tests_samplers.R`, and `tests_tests.R` to avoid spurious failures on some platforms.
