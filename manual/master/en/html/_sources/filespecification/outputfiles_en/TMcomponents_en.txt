.. highlight:: none

TMcomponents.dat
----------------

This file is the outputted files for the components of the tridiagonal
matrix and the norm of the excited state to recalculate the dynamical
Green’s function by the Lanczos method. The file format is of the binary
type. An example of the file format is as follows.

.. _file_name_22:

File name
~~~~~~~~~

*  ##_TMcomponents.dat

## indicates [string02] in a ModPara file and $$ is the number of rank.

.. _file_format_46:

File format
~~~~~~~~~~~

*  Line 1: [int01]

*  Line 2: [double01]

*  Lines 3-:
   [double02]  [double03].

.. _parameters_44:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The step for calculating dynamical Green’s
   functions by the Lanczos method :math:`N_d`.

*  [double01]

   **Type :** Double

   **Description :** The value of the norm of the excited state.

*  [double02], [double03]

   **Type :** Double

   | **Description :** The value of the components of the tridiagonal
     matrix to recalculate dynamical Green’s functions by the Lanczos
     method :math:`\alpha_i, \beta_i (i =1,\cdots N_d)`.
     [double02] is :math:`\alpha_i` and
     [double03] is :math:`\beta_i`.

.. raw:: latex

   \newpage

