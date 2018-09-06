.. highlight:: none

.. _Subsec:InputFileList:

List file for the input files
-----------------------------

This file determines the input filenames, which are correlated with the keywords. The file format is as follows.

::

    CalcMod  calcmod.def
    ModPara  modpara.def
    LocSpin  zlocspn.def
    Trans    ztransfer.def
    InterAll zinterall.def
    OneBodyG zcisajs.def
    TwoBodyG    zcisajscktaltdc.def

File format
~~~~~~~~~~~

[string01][string02]

Parameters
~~~~~~~~~~

*  [string01]

   **Type :** String

   **Description :** Select a word from keywords.

*  [string02]

   **Type :** String

   **Description :** An input filename that is correlated with the
   keywords.

Use rules
~~~~~~~~~

*  After setting keywords at [string01], the half-width state is needed
   for writing a filename. You can set the filename freely.

*  Keywords for input files are shown in :numref:`Table 4.2`

*  Essential keywords are “CalcMod", “ModPara", and “LocSpin".

*  Keywords can be set in random order.

*  If the keywords or filenames are incorrect, the program is
   terminated.

*  When the head of a line is \"#", the line is skipped.

.. _Table 4.2:
.. csv-table:: List of the definition files
    :header: "Keywords", "Details of corresponding files"
    :widths: 4, 20

    "CalcMod", "Parameters for modes of calculation"
    "ModPara", "Parameters for calculation"
    "LocSpin", "Configurations of the local spins for Hamiltonian"
    "Trans", "Transfer and chemical potential for Hamiltonian"
    "InterAll", "Two-body interactions for Hamiltonian"
    "CoulombIntra", "CoulombIntra interactions"
    "CoulombInter", "CoulombInter interactions"
    "Hund", "Hund couplings"
    "PairHop", "Pair hopping couplings"
    "Exchange", "Exchange couplings"
    "Ising", "Ising interactions"
    "PairLift", "Pair lift couplings."
    "OneBodyG", "Output components for one-body Green’s functions :math:`\langle c_{i\sigma}^{\dagger}c_{j\sigma}\rangle`"
    "TwoBodyG", "Output components for two-body Green’s functions :math:`\langle c_{i\sigma}^{\dagger}c_{j\sigma}c_{k\tau}^{\dagger}c_{l\tau}\rangle`"
    "SingleExcitation", "Operators for generating a single excited state"
    "PairExcitation", "Operators for generating a pair excited state"
    "SpectrumVec", "An input vector to calculate a restart vector"

.. raw:: latex

   \newpage