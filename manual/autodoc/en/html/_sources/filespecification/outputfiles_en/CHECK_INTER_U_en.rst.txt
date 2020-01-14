.. highlight:: none

CHECK_INTER_U.dat
-----------------

This file is outputted to check the input of the diagonal components of
the off-site interactions :math:`V_{ij}`,

.. math:: \mathcal H+=\sum_{i} V_{ij} n_{i} n_{j}

An example of the file format is as follows.

::

    i=0 isite1=1 isite2=2 tmp_V=-0.125000 
    i=1 isite1=2 isite2=3 tmp_V=-0.125000 
    i=2 isite1=3 isite2=4 tmp_V=-0.125000 
    i=3 isite1=4 isite2=5 tmp_V=-0.125000 
    i=4 isite1=5 isite2=6 tmp_V=-0.125000 
    i=5 isite1=6 isite2=1 tmp_V=-0.125000 

.. _file_format_24:

File format
~~~~~~~~~~~

*  i=[int01] isite1=[int02] 
   isite2=[int03] 
   tmp_V=[double01]

.. _parameters_24:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The counted number of inputting terms.

*  [int02], [int03]

   **Type :** Int

   | **Description :** An integer giving the site index of
     :math:`V_{ij}`.
   | [int02] and [int03] 
     correspond to :math:`i` and :math:`j`, respectively.

*  [double01]

   **Type :** Double

   **Description :** A value for :math:`V_{ij}`.

.. raw:: latex

   \newpage