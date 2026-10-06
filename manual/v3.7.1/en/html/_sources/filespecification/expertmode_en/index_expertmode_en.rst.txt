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
   :align: left

   * - File
     - Required
     - Description
   * - :doc:`List <List_file_for_the_input_files_en>`
     - Yes
     - List of input file names with keywords
   * - :doc:`CalcMod <CalcMod_file_en>`
     - Yes
     - Calculation mode settings
   * - :doc:`ModPara <ModPara_file_en>`
     - Yes
     - Basic parameters (site number, electron number, Lanczos steps, etc.)
   * - :doc:`LocSpin <LocSpin_file_en>`
     - Kondo only
     - Location of local spins

**Hamiltonian Definition**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65
   :align: left

   * - File
     - Required
     - Description
   * - :doc:`Trans <Trans_file_en>`
     - No
     - One-body terms: :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`
   * - :doc:`AnomalousTerm <AnomalousTerm_file_en>`
     - No
     - HubbardGC anomalous pair terms: :math:`c_{i\sigma_1}c_{j\sigma_2}` and :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger}`
   * - :doc:`InterAll <InterAll_file_en>`
     - No
     - General two-body interactions: :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}`
   * - :doc:`NBodyInterAll <NBodyInterAll_file_en>`
     - No
     - Generic N-body interactions: :math:`\prod_p c_{i_p\sigma'_p}^{\dagger}c_{j_p\sigma_p}`
   * - :doc:`CoulombIntra <CoulombIntra_file_en>`
     - No
     - On-site Coulomb: :math:`n_{i\uparrow}n_{i\downarrow}`
   * - :doc:`CoulombInter <CoulombInter_file_en>`
     - No
     - Off-site Coulomb: :math:`n_i n_j`
   * - :doc:`Hund <Hund_file_en>`
     - No
     - Hund coupling: :math:`n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}`
   * - :doc:`PairHop <PairHop_file_en>`
     - No
     - Pair hopping: :math:`c_{i\uparrow}^{\dagger}c_{j\uparrow}c_{i\downarrow}^{\dagger}c_{j\downarrow}`
   * - :doc:`Exchange <Exchange_file_en>`
     - No
     - Exchange coupling: :math:`c_{i\uparrow}^{\dagger}c_{j\uparrow}c_{j\downarrow}^{\dagger}c_{i\downarrow}`
   * - :doc:`Ising <Ising_file_en>`
     - No
     - Ising interaction: :math:`S_i^z S_j^z`
   * - :doc:`PairLift <PairLift_file_en>`
     - No
     - Pair lift: :math:`c_{i\uparrow}^{\dagger}c_{i\downarrow}c_{j\uparrow}^{\dagger}c_{j\downarrow}`

**Output Specifications**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65
   :align: left

   * - File
     - Required
     - Description
   * - :doc:`OneBodyG <OneBodyG_file_en>`
     - No
     - One-body Green's functions: :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}\rangle`
   * - :doc:`AnomalousG <AnomalousG_file_en>`
     - No
     - HubbardGC anomalous pair Green's functions: :math:`\langle c_{i\sigma_1}c_{j\sigma_2}\rangle` and :math:`\langle c^{\dagger}_{i\sigma_1}c^{\dagger}_{j\sigma_2}\rangle`
   * - :doc:`TwoBodyG <TwoBodyG_file_en>`
     - No
     - Two-body Green's functions: :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}c^{\dagger}_{k\sigma_3}c_{l\sigma_4}\rangle`
   * - :doc:`NBodyG <NBodyG_file_en>`
     - No
     - Generic N-body Green's functions: :math:`\left\langle \prod_p c^{\dagger}_{i_p\sigma'_p}c_{j_p\sigma_p}\right\rangle`

**Spectrum and Time Evolution**

.. list-table::
   :header-rows: 1
   :widths: 20 15 65
   :align: left

   * - File
     - Required
     - Description
   * - :doc:`SingleExcitation <SingleExcitation_file_en>`
     - Spectrum
     - Single excitation operator for dynamical Green's functions
   * - :doc:`PairExcitation <PairExcitation_file_en>`
     - Spectrum
     - Pair excitation operator for dynamical Green's functions
   * - :doc:`SingleExcitationBra <SingleExcitationBra_file_en>`
     - Optional
     - Bra single excitation operator for off-diagonal Green's functions
   * - :doc:`PairExcitationBra <PairExcitationBra_file_en>`
     - Optional
     - Bra pair excitation operator for off-diagonal Green's functions
   * - :doc:`SpectrumVec <SpectrumVec_File_en>`
     - Spectrum
     - Input vector for spectrum calculations
   * - :doc:`OneBodyTE <OneBodyTE_File_en>`
     - Time evolution
     - Time-dependent one-body terms
   * - :doc:`TwoBodyTE <TwoBodyTE_File_en>`
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

**Anomalous pair terms:**

.. toctree::
   :maxdepth: 1

   AnomalousTerm_file_en

**Two-body interactions:**

.. toctree::
   :maxdepth: 1

   InterAll_file_en
   NBodyInterAll_file_en
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
   AnomalousG_file_en
   TwoBodyG_file_en
   NBodyG_file_en

Spectrum and Time Evolution
^^^^^^^^^^^^^^^^^^^^^^^^^^^

These files are used for dynamical Green's function calculations and time evolution.

.. toctree::
   :maxdepth: 1

   SingleExcitation_file_en
   PairExcitation_file_en
   SingleExcitationBra_file_en
   PairExcitationBra_file_en
   SpectrumVec_File_en
   OneBodyTE_File_en
   TwoBodyTE_File_en
