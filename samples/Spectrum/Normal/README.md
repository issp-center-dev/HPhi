# 8-site Hubbard model on the square lattice

Compute <SzSz> spectrum at q = (pi, pi) with Shifted BiCG method.

``` bash
$ HPhi -s stan1.in
$ HPhi -s stan2.in
$ gnuplot
```

``` gnuplot
gnuplot> set xlabel "Energy"
gnuplot> set ylabel "G_{SzSz}(E)"
gnuplot> set xzeroaxis
gnuplot> plot "output/zvo_DynamicalGreen.dat" u 1:3 w l tit "Real", \
> "output/zvo_DynamicalGreen.dat" u 1:4 w l tit "Imaginary"
```
