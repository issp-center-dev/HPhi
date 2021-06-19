.. highlight:: none

.. _Subsec:spectrumvec:

SpectrumVec File 
----------------

The header of the input file for the initial vector to calculate the
dynamical Green’s functions is defined. The file name and the file
format (binary type) are as follows.

File name
~~~~~~~~~

*  ##_rank_$$.dat

## is the name of the head indicated by the key word ``SpectrumVec`` in
the ``CalcMod`` file, and $$ is the number of the rank.

.. _file_format_17:

File format
~~~~~~~~~~~

*  Line 1: [int01]

*  Line 2: [int02]

*  Lines 3 -:
   [double01]  [double02].

.. _parameters_17:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The total number of the targets of the Hilbert
   spaces.

*  [int02]

   **Type :** Int

   **Description :** The number of Lanczos or TPQ steps.

*  [double01], [double02]

   **Type :** double

   | **Description :** The coefficient value of the input vector.
   | [double01] is a real part and
     [double02] is an imaginary part.

.. raw:: latex

   \newpage