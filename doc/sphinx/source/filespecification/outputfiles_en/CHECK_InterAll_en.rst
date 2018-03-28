.. highlight:: none

CHECK_InterAll.dat
------------------

This file is outputted to check the input of the diagonal components of
general two-body interactions,

.. math:: 

   \mathcal H+=\sum_{i,j, \sigma} I_{iijj\sigma_1\sigma_1\sigma_2\sigma_2} c_{i\sigma_1}^{\dagger}c_{i\sigma_1}c_{i\sigma_2}^{\dagger}c_{i\sigma_2}.

An example of the file format is as follows.

::

    i=0 isite1=1 A_spin=0 isite2=2 B_spin=0 tmp_V=0.500000 
    i=1 isite1=1 A_spin=0 isite2=2 B_spin=1 tmp_V=-0.500000 
    i=2 isite1=1 A_spin=1 isite2=2 B_spin=0 tmp_V=-0.500000 
    i=3 isite1=1 A_spin=1 isite2=2 B_spin=1 tmp_V=0.500000 
    i=4 isite1=2 A_spin=0 isite2=3 B_spin=0 tmp_V=0.500000 
    i=5 isite1=2 A_spin=0 isite2=3 B_spin=1 tmp_V=-0.500000 
    ...

.. _file_format_21:

File format
~~~~~~~~~~~

*  i=[int01] isite1=[int02] 
   A_spin=[int03] 
   isite2=[int04] 
   B_spin=[int05] 
   tmp_V=[double01]

.. _parameters_21:

Parameters
~~~~~~~~~~

*  [int01]

   **Type :** Int

   **Description :** The counted number of inputting terms.

*  [int02], [int04]

   **Type :** Int

   | **Description :** An integer for showing the site index of
     :math:`I_{iijj\sigma_1\sigma_1\sigma_2\sigma_2}`.
   | [int02] and [int04]
     correspond to :math:`i` and :math:`j`, respectively.

*  [int03], [int05]

   **Type :** Int

   | **Description :** An integer for showing the spin index of
     :math:`I_{iijj\sigma_1\sigma_1\sigma_2\sigma_2}`:
   | 0: Up-spin
   | 1: Down-spin.
   | [int03] and [int05]
     correspond to :math:`\sigma_1` and :math:`\sigma_2`, respectively.

*  [double01]

   **Type :** Double

   **Description :** A value for
   :math:`I_{iijj\sigma_1\sigma_1\sigma_2\sigma_2}`.

.. raw:: latex

   \newpage