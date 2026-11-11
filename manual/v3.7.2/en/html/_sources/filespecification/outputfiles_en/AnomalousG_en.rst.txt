.. highlight:: none

.. _Subsec:anomalousgoutput:

AnomalousG.dat
--------------

This file is the output file for anomalous pair Green's functions specified
by the input file with the keyword ``AnomalousG``. A line with type ``0`` is
the expectation value of :math:`c_{i\sigma_1}c_{j\sigma_2}`, and a line with
type ``1`` is the expectation value of
:math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger}`. An example of the
file format is as follows.

::

    0 0 0 0 1 -0.5000000000 0.0000000000
    1 0 1 0 0 -0.5000000000 0.0000000000

File name
~~~~~~~~~

*  Lanczos method, LOBCG method: ##_AnomalousG.dat

*  TPQ method, cTPQ method: ##_AnomalousG_set??step%%.dat

*  TPQ/cTPQ method with ``OutputGreenFormat=1``:
   ##_AnomalousG_tpq.dat

*  Full diagonalization method: ##_AnomalousG_eigen&&.dat

*  Full diagonalization method, LOBCG method with ``OutputGreenFormat=1``:
   ##_AnomalousG_eigen.dat

*  Real time evolution method: ##_AnomalousG_step%%.dat

*  Real time evolution method with ``OutputGreenFormat=1``:
   ##_AnomalousG_te.dat

##, ??, %%, and && indicate [string02] in a ModPara file, the number of
runs in calculation in the TPQ method, the number of steps in the TPQ or
time evolution method, and the index of the eigenvalues, respectively.
In aggregate files, each data row has leading index columns: ``set step``
for TPQ/cTPQ, ``step`` for real-time evolution, and ``eigen`` for Full
diagonalization/LOBCG.

File format
~~~~~~~~~~~

*  [int01] [int02] [int03] [int04] [int05] [double01] [double02].

For aggregate files generated with ``OutputGreenFormat=1``, the same
columns are preceded by ``set step``, ``step``, or ``eigen`` according to
the calculation method.

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   | **Description :** The anomalous pair type:
   | 0: :math:`\langle c_{i\sigma_1}c_{j\sigma_2}\rangle`
   | 1: :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger}\rangle`.

*  [int02], [int04]

   **Type :** Int

   **Description :** The integer of the site number.
   [int02] and [int04] show the :math:`i` and :math:`j` site numbers,
   respectively.

*  [int03], [int05]

   **Type :** Int

   | **Description :** The integer of the spin index:
   | 0: Up-spin
   | 1: Down-spin.
   | [int03] and [int05] show
     :math:`\sigma_1` and :math:`\sigma_2`, respectively.

*  [double01], [double02]

   **Type :** Double

   | **Description :** The value of the anomalous pair Green's function.
   | [double01] and [double02] show the real and imaginary parts,
     respectively.

.. raw:: latex

   \newpage
