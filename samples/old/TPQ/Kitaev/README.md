# 12-site Kitaev model

Compute the specific heat with the TPQ method.

``` bash
$ HPhi -s stan.in
$ python AveSSrand.py -n 5
$ gnuplot
```

``` gnuplot
gnuplot> set xlabel "Temperature"
gnuplot> set ylabel "Spesific heat"
gnuplot> plot [0:1] "ave_TPQ.dat" u 1:5:2:6 w xyerrorlines
```
