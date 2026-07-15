.. highlight:: none

.. _Subsec:dynamicalG:

DynamicalGreen.dat
------------------

This file is the outputted file for calculating the dynamical Green’s
function. An example of the file format is as follows.

.. _file_name_20:

File name
~~~~~~~~~

*  ##_DynamicalGreen.dat

## indicates [string02] in a ModPara file.

When the finite-temperature eigenstate loop and the multiple-operator /
multiple-bra modes are enabled (``SpectrumLoopExct``, ``SpectrumNumOp``,
``SpectrumNumBra`` in the :ref:`ModPara <Subsec:modpara>` file), one file is
written per ``(eigenstate, ket operator, bra)`` combination and the
zero-based indices are appended to the file name:

*  ##_DynamicalGreen\_\ *idx*\ .dat — eigenstate *idx*
   (``SpectrumLoopExct``).

*  ##_DynamicalGreen\_\ *idx*\ \_\ *op*\ .dat — eigenstate *idx*, ket operator
   set *op* (``SpectrumNumOp``).

*  ##_DynamicalGreen\_\ *idx*\ \_\ *op*\ \_\ *bra*\ .dat — eigenstate *idx*,
   ket operator set *op*, bra operator set *bra* (``SpectrumNumBra``).

The file format described below is identical for all of these variants.

.. _file_format_44:

File format
~~~~~~~~~~~

*  Lines 1-:
   [double01]  [double02]  [double03]  [double04]

.. _parameters_42:

Parameters
~~~~~~~~~~

*  [double01], [double02]

   **Type :** Double

   | **Description :** The value of the frequency :math:`z`.
   | [double01] and [double02]
     are a real and an imaginary part of :math:`z`, respectively.

*  [double03], [double04]

   **Type :** Double

   | **Description:** The value of dynamical Green’s functions
     :math:`G(z)`.
   | [double03] and [double04]
     are a real and an imaginary part of :math:`G(z)`, respectively.

.. raw:: latex

   \newpage