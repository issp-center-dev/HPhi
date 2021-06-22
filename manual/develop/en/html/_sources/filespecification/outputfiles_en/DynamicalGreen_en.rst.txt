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