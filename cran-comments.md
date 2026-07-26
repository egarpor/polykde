## Test environments

* local R installation, R 4.5.2
* check_win_devel()

## R CMD check results

0 errors | 1 warnings | 0 notes

─  checking whether package ‘polykde’ can be installed ... [11s/12s] WARNING (12.2s)
   Found the following significant warnings:
     /Library/Frameworks/R.framework/Resources/include/R_ext/Boolean.h:62:36: warning: unknown warning group '-Wfixed-enum-extension', ignored [-Wunknown-warning-option]
   See ‘/private/var/folders/q8/d13pdb9s7b31x54qxcywwy1w0000gn/T/RtmpCa1ipF/file63f724b9f72f/polykde.Rcheck/00install.out’ for details.


## Comments

In this release the flaky tests that were causing problems with the update of rotasym have been fixed.

The warning is benign.