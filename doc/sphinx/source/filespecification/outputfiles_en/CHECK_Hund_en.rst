.. highlight:: none

CHECK_Hund.dat
--------------

This file is outputted to check the input of the Hund couplings
:math:`J_{ij}^{\rm Hund}`,

.. math:: \mathcal H += -\sum_{i,j}J_{ij}^{\rm Hund} (n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}).

An example of the file format is as follows.

::

    i=0 isite1=1 isite2=2 tmp_V=0.250000 
    i=1 isite1=2 isite2=3 tmp_V=0.250000 
    i=2 isite1=3 isite2=4 tmp_V=0.250000 
    i=3 isite1=4 isite2=5 tmp_V=0.250000 
    i=4 isite1=5 isite2=6 tmp_V=0.250000 
    i=5 isite1=6 isite2=1 tmp_V=0.250000 

.. _file_format_23:

File format
~~~~~~~~~~~

*  i=[int01] isite1=[int02] 
   isite2=[int03] 
   tmp_V=[double01]

.. _parameters_23:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The counted number of inputting terms.

*  [int02], [int03]

   **Type :** Int

   | **Description :** An integer for showing the site index of
     :math:`J_{ij}^{\rm Hund}`.
   | [int02] and [int03]
     correspond to :math:`i` and :math:`j`, respectively.

*  [double01]

   **Type :** Double

   **Description :** A value for :math:`J_{ij}^{\rm Hund}`.

.. raw:: latex

   \newpage