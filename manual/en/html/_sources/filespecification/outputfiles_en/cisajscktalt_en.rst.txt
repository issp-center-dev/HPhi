.. highlight:: none

.. _Subsec:cisajscktalt:

cisajscktalt.dat
----------------

This file is the outputted files for the two-body Green’s function
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle`.
The target components are set in the input file with the keyword
“TwoBodyG". An example of the file format is as follows.

::

        0    0    0    0    0    0    0    0 0.4452776740 0.0000000000
        0    0    0    0    0    1    0    1 0.1843355815 0.0000000000
        0    0    0    0    1    0    1    0 0.1812412105 0.0000000000
        0    0    0    0    1    1    1    1 0.2640364635 0.0000000000
        0    0    0    0    2    0    2    0 0.0279690007 0.0000000000
        0    0    0    0    2    1    2    1 0.2009271524 0.0000000000
        0    0    0    0    3    0    3    0 0.2512810778 0.0000000000
        0    0    0    0    3    1    3    1 0.1939965962 0.0000000000
        ...

.. _file_name_17:

File name
~~~~~~~~~

*  Lanczos method: ##_cisajscktalt.dat

*  TPQ method: ##_cisajscktalt_set??step%%.dat

*  Full diagonalization method, LOBCG method:
   ##_cisajscktalt_eigen&&.dat

*  Real time evolution method: ## cisajscktalt step%%.dat

##, ??, %%, and && indicate [string02] in a ModPara file, the number of
runs in calculation in the TPQ method, the number of steps in the TPQ
method, and the index of the eigenvalues, respectively.

.. _file_format_41:

File format
~~~~~~~~~~~

*  [int01]  [int02]  [int03]  [int04]  [int05]  [int06]  [int07]  [int08]  [double01]  [double02].

.. _parameters_41:

Parameters
~~~~~~~~~~

*  [int01],
   [int03],[int05],
   [int07]

   **Type :** Int

   **Description :** The integer of the site number.
   [int01], [int03],
   [int05], and [int07] show the
   :math:`i`, :math:`j`, :math:`k`, and :math:`l` site numbers,
   respectively.

*  [int02],
   [int04],[int06],
   [int08]

   **Type :** Int

   | **Description :** The integer of the spin index:
   | 0: Up-spin
   | 1: Down-spin.
   | [int02], [int04],
     [int06], and [int08] show
     :math:`\sigma_1`, :math:`\sigma_2`, :math:`\sigma_3`, and
     :math:`\sigma_4`, respectively.

*  [double01], [double02]

   **Type :** Double

   | **Description :** The value of
     :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle`.
   | [double01] and [double02]
     show the real and imaginary part of
     :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle`,
     respectively.

.. raw:: latex

   \newpage