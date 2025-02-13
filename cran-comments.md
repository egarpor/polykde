## Test environments

* local R installation, R 4.2.2
* win-builder (release, devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

Fixed the following issues signaled by CRAN team:

- Please omit the redundant "Tools for" from the description.
- Missing Rd-tags:
  r_kde_polysph.Rd: \value
  r_kern_polysph.Rd: \value
  r_unif_polysph.Rd: \value
  r_vmf_polysph.Rd: \value
- Please omit one colon. -> Warning: Used ::: in documentation:
      man/AP.Rd:
         polykde:::AP(x = x, v = v, ind_dj = comp_ind_dj(d))
      man/beta0_R.Rd:
         polykde:::beta0_R(d = d, R = R)
      man/bind_lists.Rd:
         polykde:::bind_lists(lists = lists, bind = "rbind")
      man/diamond_crossprod.Rd:
         polykde:::diamond_crossprod(X = X, ind_dj = comp_ind_dj(d))
      man/hd_mc.Rd:
         polykde:::hd_mc(log_f = log_f, log_g = log_f, d = d)
      man/hd_mc.Rd:
         polykde:::hd_mc(log_f = log_f, log_g = log_g, d = d)
      man/hd_mc.Rd:
         polykde:::hd_mc(log_f = log_f, log_g = log_f, d = d, bhatta = TRUE)
      man/hd_mc.Rd:
         polykde:::hd_mc(log_f = log_f, log_g = log_g, d = d, bhatta = TRUE)
      man/J_d_k.Rd:
         polykde:::J_d_k(d = 1:5, k = 10)
      man/polylog_minus_exp_mu.Rd:
         polykde:::polylog_minus_exp_mu(mu = 1:5, s = 1)
      man/polylog_minus_exp_mu.Rd:
         polykde:::polylog_minus_exp_mu(mu = 1, s = 1:5)
      man/polylog_minus_exp_mu.Rd:
         polykde:::polylog_min [... truncated]
