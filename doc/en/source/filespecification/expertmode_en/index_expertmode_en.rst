.. highlight:: none

.. _Ch:HowToExpert:

Input files for *Expert* mode
=============================

This section explains the input files for the expert mode.

Quick Reference
---------------

The following table summarizes all input files for expert mode.

**Basic Settings**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - File
     - Required
     - Description
   * - List
     - Yes
     - List of input file names with keywords
   * - CalcMod
     - Yes
     - Calculation mode settings
   * - ModPara
     - Yes
     - Basic parameters (site number, electron number, Lanczos steps, etc.)
   * - LocSpin
     - Kondo only
     - Location of local spins

**Hamiltonian Definition**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - File
     - Required
     - Description
   * - Trans
     - No
     - One-body terms: :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`
   * - InterAll
     - No
     - General two-body interactions: :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}`
   * - CoulombIntra
     - No
     - On-site Coulomb: :math:`n_{i\uparrow}n_{i\downarrow}`
   * - CoulombInter
     - No
     - Off-site Coulomb: :math:`n_i n_j`
   * - Hund
     - No
     - Hund coupling: :math:`n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}`
   * - PairHop
     - No
     - Pair hopping: :math:`c_{i\uparrow}^{\dagger}c_{j\uparrow}c_{i\downarrow}^{\dagger}c_{j\downarrow}`
   * - Exchange
     - No
     - Exchange coupling: :math:`c_{i\uparrow}^{\dagger}c_{j\uparrow}c_{j\downarrow}^{\dagger}c_{i\downarrow}`
   * - Ising
     - No
     - Ising interaction: :math:`S_i^z S_j^z`
   * - PairLift
     - No
     - Pair lift: :math:`c_{i\uparrow}^{\dagger}c_{i\downarrow}c_{j\uparrow}^{\dagger}c_{j\downarrow}`

**Output Specifications**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - File
     - Required
     - Description
   * - OneBodyG
     - No
     - One-body Green's functions: :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}\rangle`
   * - TwoBodyG
     - No
     - Two-body Green's functions: :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}c^{\dagger}_{k\sigma_3}c_{l\sigma_4}\rangle`

**Spectrum and Time Evolution**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - File
     - Required
     - Description
   * - SingleExcitation
     - Spectrum
     - Single excitation operator for dynamical Green's functions
   * - PairExcitation
     - Spectrum
     - Pair excitation operator for dynamical Green's functions
   * - SpectrumVec
     - Spectrum
     - Input vector for spectrum calculations
   * - OneBodyTE
     - Time evolution
     - Time-dependent one-body terms
   * - TwoBodyTE
     - Time evolution
     - Time-dependent two-body interactions

----

Detailed Specifications
-----------------------

Basic Settings
^^^^^^^^^^^^^^

These files define the fundamental parameters for the calculation.

.. toctree::
   :maxdepth: 1

   List_file_for_the_input_files_en
   CalcMod_file_en
   ModPara_file_en
   LocSpin_file_en

Hamiltonian Definition
^^^^^^^^^^^^^^^^^^^^^^

These files specify the terms in the Hamiltonian.

**One-body terms:**

.. toctree::
   :maxdepth: 1

   Trans_file_en

**Two-body interactions:**

.. toctree::
   :maxdepth: 1

   InterAll_file_en
   CoulombIntra_file_en
   CoulombInter_file_en
   Hund_file_en
   PairHop_file_en
   Exchange_file_en
   Ising_file_en
   PairLift_file_en

Output Specifications
^^^^^^^^^^^^^^^^^^^^^

These files specify which physical quantities to calculate and output.

.. toctree::
   :maxdepth: 1

   OneBodyG_file_en
   TwoBodyG_file_en

Spectrum and Time Evolution
^^^^^^^^^^^^^^^^^^^^^^^^^^^

These files are used for dynamical Green's function calculations and time evolution.

.. toctree::
   :maxdepth: 1

   SingleExcitation_file_en
   PairExcitation_file_en
   SpectrumVec_File_en
   OneBodyTE_File_en
   TwoBodyTE_File_en
