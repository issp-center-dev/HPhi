.. highlight:: none

.. _Subsec:checkchemi:

CHECK_Chemi.dat
---------------

This file is outputted to check the input of chemical potential
:math:`\mu_{i\sigma}`,

.. math:: \mathcal H+=\sum_{i,\sigma} \mu_{i\sigma} c_{i\sigma}^{\dagger}c_{i\sigma}.

An example of the file format is as follows.

::

    i=0 spin=0 isite1=1 tmp_V=0.000000 
    i=1 spin=0 isite1=2 tmp_V=0.000000 
    i=2 spin=0 isite1=3 tmp_V=0.000000 
    i=3 spin=0 isite1=4 tmp_V=0.000000 
    i=4 spin=0 isite1=5 tmp_V=0.000000 
    i=5 spin=0 isite1=6 tmp_V=0.000000 
    ...

.. _file_format_20:

File format
~~~~~~~~~~~

*  i=[int01] spin=[int02] isite1=[int03] tmp_V=[double01]

.. _parameters_20:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The counted number of inputting terms.

*  [int02]

   **Type :** Int

   | **Description :** An integer for showing the spin index of
     :math:`\mu_{i\sigma}`:
   | 0: Up-spin,
   | 1: Down-spin.

*  [int03]

   **Type :** Int

   **Description :** An integer for showing the site index of
   :math:`\mu_{i\sigma}`.

*  [double01]

   **Type :** Double

   **Description :** A value for :math:`\mu_{i\sigma}`.

.. raw:: latex

   \newpage