## Test environments

* local R installation, R 4.5.2
* check_win_devel()

## R CMD check results

0 errors | 1 warnings | 0 notes

W  checking whether package ‘rotasym’ can be installed (4s)
   Found the following significant warnings:
     /Library/Frameworks/R.framework/Resources/include/R_ext/Boolean.h:62:36: warning: unknown warning group '-Wfixed-enum-extension', ignored [-Wunknown-warning-option]
   See ‘/private/var/folders/q8/d13pdb9s7b31x54qxcywwy1w0000gn/T/RtmpOH8LUn/file154312a14c26c/rotasym.Rcheck/00install.out’ for details.

## Comments

The URL https://doi.org/10.1080/01621459.2025.2521898 has been double-checked.

The warning from Boolean.h is out of my control.

Revdep tested with revdepcheck::revdep_check(pkg, num_workers = 1).
