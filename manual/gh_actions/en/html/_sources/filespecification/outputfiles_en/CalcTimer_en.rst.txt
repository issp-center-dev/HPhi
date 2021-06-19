.. highlight:: none

CalcTimer.dat
-------------

| The name of the calculation process, the process number, and the
  calculation process time are outputted in order at each line in the
  CalcTime.dat file. An example of the file format for the TPQ method is
  as follows.

::

    All                                                [0000]     12.94052
      sz                                               [1000]      0.01795
      diagonalcalc                                     [2000]      0.00693
      CalcByTPQ                                        [3000]     12.90670
        FirstMultiply                                  [3100]      0.08416
          rand   in FirstMultiply                      [3101]      0.00172
          mltply in FirstMultiply                      [3102]      0.07707
        expec_energy_flct                              [3200]      9.06255
          calc flctuation in expec_energy_flct         [3201]      1.67779
          mltply in expec_energy_flct                  [3202]      7.31207
        expec_onebody                                  [3300]      0.11640
        expec_twobody                                  [3400]      3.28796
        Multiply                                       [3500]      0.14840
        FileIO                                         [3600]      0.20493
    ================================================
    All mltply                                         [0001]      7.38883
      diagonal                                         [0100]      0.04153
      Hubbard                                          [0300]      7.34636
        trans    in Hubbard                            [0310]      7.34595
          double                                       [0311]      0.00000
          single                                       [0312]      0.00000
          inner                                        [0313]      7.34299
        interall in Hubbard                            [0320]      0.00008
          interPE                                      [0321]      0.00000
          inner                                        [0322]      0.00000
        pairhopp in Hubbard                            [0330]      0.00006
          interPE                                      [0331]      0.00000
          inner                                        [0332]      0.00000
        exchange in Hubbard                            [0340]      0.00004
          interPE                                      [0341]      0.00000
          inner                                        [0342]      0.00000
    ================================================

.. raw:: latex

   \newpage