.. highlight:: none

recalcvec.dat
-------------

This file is the outputted file for two vectors to recalculate the
dynamical Green’s function by the Lanczos method. The file format is of
the binary type. An example of the file format is as follows.

.. _file_name_21:

File name
~~~~~~~~~

*  ##_recalcvec_rank_$$.dat

## indicates [string02] in a ModPara file and $$ is the number of rank.

.. _file_format_45:

File format
~~~~~~~~~~~

*  Line 1: [int01]

*  Line 2: [int02]

*  Lines 3 - 3+ [int02]:
   [double01]  [double02]

*  Lines 4+ [int02] -
   4+2 :math:`\times` [int02]:
   [double03]  [double04].

.. _parameters_43:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The step for calculating dynamical Green’s
   functions by the Lanczos method :math:`N_d`.

*  [int02]

   **Type :** Long Int

   **Description :** The total number of targets of the Hilbert spaces.

*  [double01], [double02]

   **Type :** Double

   | **Description :** The value of the vector :math:`\boldsymbol{v}_{k+1}` for
     recalculating dynamical Green’s functions by the Lanczos method.
   | [double01] and [double02]
     are a real part and an imaginary part of :math:`\boldsymbol{v}_{k+1}`,
     respectively. The fist component is not used for calculation.

*  [double03], [double04]

   **Type :** Double

   | **Description :** The value of the vector :math:`\boldsymbol{v}_{k}` for
     recalculating dynamical Green’s functions by the Lanczos method.
   | [double03] and [double04]
     are a real part and an imaginary part of :math:`\boldsymbol{v}_{k}`,
     respectively. The fist component is not used for calculation.

.. raw:: latex

   \newpage