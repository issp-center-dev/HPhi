# 14-site Heisenberg model on Chain lattice

Compute the Sz-Sz dynamical correlation function
and plot the imaginary part.

``` bash
$ HPhi -s stan.in
$ echo "2 6
G 0.0 0.0 0.0
X 0.0 0.5 0.0
" >> geometry.dat
$ dynamicalr2k namelist.def geometry.dat
$ gnuplot
```

``` gnuplot
gnuplot> load "kpath.gp"
gnuplot> splot [][][0:] "output/zvo_dyn.dat" u 1:2:(-$4) w l
```
