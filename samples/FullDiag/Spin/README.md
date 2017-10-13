# 12-site Heisenberg model on the Kagome lattice

Compute all states with the full diagonalization

``` bash
$ HPhi -s stan.in
$ histgram.sh 0.1 > hist.dat
$ gnuplot
```

``` gnuplot
gnuplot> set xlabel "Energy"
gnuplot> set ylabel "Number of states"
gnuplot> plot "hist.dat" u 1:2 w histeps
```
