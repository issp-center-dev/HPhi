.. highlight:: none

.. _Subsec:residual:

residual.dat
~~~~~~~~~~~~~

This file is outputted to show the calculation process information for the dynamical Green's function by using the LOBCG method (\ ``CalcType``\ =3 in ``CalcMod``\  file). 
An example of the file format is as follows.

File name
^^^^^^^^^^

-  residual.dat

File format
^^^^^^^^^^^^

-  Lines 1-\ ``NOmega``\ :\ :math:`[`\ int01\ :math:`]` 
   :math:`[`\ double01\ :math:`]`  :math:`[`\ double02\ :math:`]`
   :math:`[`\ double03\ :math:`]`  :math:`[`\ double04\ :math:`]`

-  Line \ ``NOmega``\ +1 :empty

-  Repeat Lines 1-\ ``NOmega``\ +1.

Parameters
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **Type :** Int

   **Description :** Number of Iteration Steps. Only multiples of 10 are printed.

-  :math:`[`\ double01\ :math:`]`

   **Type :** Double

   **Description :** The value of frequencies.

-  :math:`[`\ double02\ :math:`]`, :math:`[`\ double03\ :math:`]`

   **Type :** Double

   | **Description :** The value of dynamical Green’s functions
     :math:`G(z)` at \ :math:`[`\ int01\ :math:`]`\  iteration step.
   | [double02] and [double03]
     are a real and an imaginary part of :math:`G(z)`, respectively.

-  :math:`[`\ double04\ :math:`]`

   **Type :** Double

   | **Description :** The norm of the residual vector at \ :math:`[`\ int01\ :math:`]`\  iteration step.

.. raw:: latex

   \newpage
