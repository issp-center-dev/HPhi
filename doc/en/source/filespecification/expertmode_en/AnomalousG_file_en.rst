.. highlight:: none

.. _Subsec:anomalousg:

AnomalousG file
---------------

This file determines the target components of anomalous pair Green's
functions for HubbardGC calculations. A line with type ``0`` represents
:math:`\langle c_{i\sigma_1}c_{j\sigma_2}\rangle`, and a line with type
``1`` represents
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger}\rangle`. An
example of the file format is as follows.

::

    ========================
    NAnomalousG 2
    ========================
    ========AnomalousG=======
    ========================
    0 0 0 0 1
    1 0 1 0 0

File format
~~~~~~~~~~~

*  Line 1: Header

*  Line 2: [string01] [int01]

*  Lines 3-5: Header

*  Lines 6-: [int02] [int03] [int04] [int05] [int06].

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String (a blank parameter is not allowed)

   **Description :** A keyword for the total number of anomalous pair
   Green's functions. You can freely give a name to the keyword.

*  [int01]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving the total number of anomalous pair
   Green's functions.

*  [int02]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** The anomalous pair type:
   | 0: :math:`\langle c_{i\sigma_1}c_{j\sigma_2}\rangle`
   | 1: :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger}\rangle`.

*  [int03], [int05]

   **Type :** Int (a blank parameter is not allowed)

   **Description :** An integer giving a site index
   (:math:`0<=` [int03], [int05] :math:`<` ``Nsite``).

*  [int04], [int06]

   **Type :** Int (a blank parameter is not allowed)

   | **Description :** An integer giving a spin index:
   | 0: Up-spin,
   | 1: Down-spin.

Use rules
~~~~~~~~~

*  Headers cannot be omitted.

*  This file is currently supported only for HubbardGC calculations.

*  MPI execution is not supported for ``AnomalousG``.

*  ``AnomalousG`` is not supported for ``CalcSpec`` calculations.

*  The same fermion operator cannot be used twice in one anomalous pair.

*  A program is terminated when [int01] is different from the total number
   of anomalous pair Green's functions defined in this file.

*  A program is terminated when [int02]-[int06] are outside the range of the
   defined values.

.. raw:: latex

   \newpage
