## Test environments

* local R installation, R 4.5.2
* check_win_devel()

## R CMD check results

0 errors | 1 warnings | 1 notes

─  checking whether package ‘polykde’ can be installed ... [11s/12s] WARNING (12.2s)
   Found the following significant warnings:
     /Library/Frameworks/R.framework/Resources/include/R_ext/Boolean.h:62:36: warning: unknown warning group '-Wfixed-enum-extension', ignored [-Wunknown-warning-option]
   See ‘/private/var/folders/q8/d13pdb9s7b31x54qxcywwy1w0000gn/T/RtmpCa1ipF/file63f724b9f72f/polykde.Rcheck/00install.out’ for details.

❯ checking CRAN incoming feasibility ... [3s/22s] NOTE
  Maintainer: ‘Eduardo García-Portugués <edgarcia@est-econ.uc3m.es>’

  Days since last update: 1

## Comments

Submission to fix a failing test in <https://cran.r-project.org/web/checks/check_results_polykde.html>, in Additional Issues.