This directory contains scripts used to generate some extra results. These are versions of the standard trio of isoQTL scripts -- isoqtl.py (the user-facing script), cis\_pass.pyx (cython script for running regressions and computing results), and setup.py (used to compile cis\_pass.pyx) -- in addition to a few helper scripts. The files are as follows:

- cis\_pass\_opposites.pyx + isoqtl\_opposites.py + setup\_opposites.py: assessment of how often genes found by isoform-aware methods but not sum/average methods (as determined by compute\_diff.py) had isoforms with opposite-sign effect sizes (as determined by compute\_opposite.py)

- cis\_pass\_printall.pyx + isoqtl\_printall.py + setup\_printall.py: same as cis\_pass.pyx except it prints all isoform-SNP associations (not just the most significant), which was used to create QQ plots (new\_qq\_plot.ipynb)

- cis\_pass\_timing.pyx + isoqtl\_timing.py + setup\_timing.py: used to obtain timing results for non-QTLtools methods

Please see the README in the top level directory and the manuscript for more details.

