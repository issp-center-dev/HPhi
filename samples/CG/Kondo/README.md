# 6-site Kondo lattice model on the triangular lattice

Compute the ground state and plot correlation function.

``` bash
$ HPhi -s stan.in
$ gnuplot lattice.gp
$ fourier namelist.def geometry.dat
$ corplot output/zvo_corr_eigen0.dat
```
