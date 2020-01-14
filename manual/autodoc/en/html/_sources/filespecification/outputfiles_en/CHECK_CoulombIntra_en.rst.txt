.. highlight:: none

CHECK_CoulombIntra.dat
----------------------

This file is outputted to check the input of the on-site interactions
:math:`U_i`,

.. math:: 

   \mathcal H+=\sum_{i} U_i n_{i\uparrow} n_{j \downarrow}.

An example of the file format is as follows.

::

    i=0 isite1=1 tmp_V=4.000000 
    i=1 isite1=2 tmp_V=4.000000 
    i=2 isite1=3 tmp_V=4.000000 
    i=3 isite1=4 tmp_V=4.000000 
    i=4 isite1=5 tmp_V=4.000000 
    i=5 isite1=6 tmp_V=4.000000 

.. _file_format_22:

File format
~~~~~~~~~~~

*  i=[int01] isite1=[int02] 
   tmp_V=[double01]

.. _parameters_22:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The counted number of inputting terms.

*  [int02]

   **Type :** Int

   **Description :** An integer for showing the site index of
   :math:`U_i`.

*  [double01]

   **Type :** Double

   **Description :** A value for :math:`U_i`.

.. raw:: latex

   \newpage