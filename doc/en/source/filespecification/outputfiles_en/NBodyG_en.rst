.. highlight:: none

.. _Subsec:nbodygdat:

NBodyG.dat
----------

This file contains the generic N-body Green's functions specified by the
NBodyG input file. A line with :math:`N` factors stores

.. math::

   \left\langle \prod_{p=1}^{N}
   c_{i_p\sigma'_p}^{\dagger} c_{j_p\sigma_p} \right\rangle .

An example of the file format is as follows.

::

    1    0    0    1    0 0.1250000000 0.0000000000
    2    0    0    0    0    1    1    1    1 0.2500000000 0.0000000000
    3    0    0    0    0    1    1    1    1    2    0    2    0 0.0312500000 0.0000000000

.. _file_name_nbodygdat:

File name
~~~~~~~~~

*  Lanczos method: ##_NBodyG.dat

*  TPQ method: ##_NBodyG_set??step%%.dat

*  TPQ/cTPQ method with ``OutputGreenFormat=1``: ##_NBodyG_tpq.dat

*  Full diagonalization method, LOBCG method: ##_NBodyG_eigen&&.dat

*  Full diagonalization method, LOBCG method with ``OutputGreenFormat=1``:
   ##_NBodyG_eigen.dat

*  Real time evolution method: ##_NBodyG_step%%.dat

*  Real time evolution method with ``OutputGreenFormat=1``:
   ##_NBodyG_te.dat

##, ??, %%, and && indicate [string02] in a ModPara file, the number of
runs in calculation in the TPQ method, the number of steps in the TPQ or
real-time evolution method, and the index of the eigenvalues,
respectively. In aggregate files, each data row has leading index
columns: ``set step`` for TPQ/cTPQ, ``step`` for real-time evolution, and
``eigen`` for Full diagonalization/LOBCG.

.. _file_format_nbodygdat:

File format
~~~~~~~~~~~

*  [int01] ([int02] [int03] [int04] [int05]) ... [double01] [double02].
   The group in parentheses is repeated [int01] times.

For aggregate files generated with ``OutputGreenFormat=1``, the same
columns are preceded by ``set step``, ``step``, or ``eigen`` according to
the calculation method.

.. _parameters_nbodygdat:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The number of factors :math:`N` in this component.

*  [int02], [int04] of each factor

   **Type :** Int

   **Description :** Site indices. [int02] and [int04] show
   :math:`i_p` and :math:`j_p`, respectively.

*  [int03], [int05] of each factor

   **Type :** Int

   | **Description :** Spin or local-state indices.
   | For Hubbard, tJ, Kondo, and their GC/NConserved variants:
   | 0: Up-spin
   | 1: Down-spin.
   | For SpinlessFermion/SpinlessFermionGC, only 0 is output.

*  [double01], [double02]

   **Type :** Double

   | **Description :** The value of the specified generic N-body
     Green's function.
   | [double01] and [double02] show the real and imaginary parts,
     respectively.

.. raw:: latex

   \newpage
