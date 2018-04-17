# 8-site Hubbard model on square lattice

Compute the ground state and plot correlation function.

``` bash
$ HPhi -s stan.in
$ gnuplot lattice.gp
$ echo "4 20
G 0 0 0
X 0.5 0 0
M 0.5 0.5 0
G 0 0 0
16 16 1" >> geometry.dat
$ greenr2k namelist.def geometry.dat
$ gnuplot
```

``` gnuplot
gnuplot> load "kpath.gp"
gnuplot> plot "output/zvo_corr_eigen0.dat" u 1:12 w l
```
