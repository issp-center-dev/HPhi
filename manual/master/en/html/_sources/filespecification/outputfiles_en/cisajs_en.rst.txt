.. highlight:: none

.. _Subsec:cgcisajs:

cisajs.dat
----------

This file is the outputted files for one-body Green’s function
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`. The target
components are set in the input file with the keyword “OneBodyG". An
example of the file format is as follows.

::

        0    0    0    0 0.4452776740 0.0000000000
        0    1    0    1 0.4452776740 0.0000000000
        1    0    1    0 0.5000000000 0.0000000000
        1    1    1    1 0.5000000000 0.0000000000
        2    0    2    0 0.4452776740 0.0000000000
        2    1    2    1 0.4452776740 0.0000000000
        3    0    3    0 0.5000000000 0.0000000000
        3    1    3    1 0.5000000000 0.0000000000
        ...

.. _file_name_16:

File name
~~~~~~~~~

*  Lanczos method: ##_cisajs.dat

*  TPQ method: ##_cisajs_set??step%%.dat

*  Full diagonalization method, LOBCG method: ##_cisajs_eigen&&.dat.

*  Real time evolution method: ## cisajs step%%.dat

##, ??, %%, and && indicate [string02] in a ModPara file, the number of
runs in calculation in the TPQ method, the number of steps in the TPQ
method, and the index of the eigenvalues, respectively.

.. _file_format_40:

File format
~~~~~~~~~~~

*  [int01]  [int02]  [int03]  [int04]  [double01]  [double02]

.. _parameters_40:

Parameters
~~~~~~~~~~

*  [int01], [int03]

   **Type :** Int

   **Description :** The integer of the site number.
   [int01] and [int03] show the
   :math:`i` and :math:`j` site numbers, respectively.

*  [int02], [int04]

   **Type :** Int

   | **Description :** The integer of the spin index:
   | 0: Up-spin
   | 1: Down-spin.
   | [int02] and [int04] show
     :math:`\sigma_1` and :math:`\sigma_2`, respectively.

*  [double01], [double02]

   **Type :** Double

   | **Description :** The value of
     :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`.
   | [double01] and [double02]
     show the real and imaginary part of
     :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`,
     respectively.

.. raw:: latex

   \newpage