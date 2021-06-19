.. highlight:: none

.. _Ch:HowToExpert:

Input files for *Expert* mode
=============================

In this section, the details of the input files for the expert mode are
explained. The input files are categorized according to the following
four parts.

(1) List:
    This file is a list of the input file names with the keywords. Each
    keyword is fixed, but the file names can be determined freely.

(2) Basic parameters:
    | The following input files give the basic parameters. The types of
      input files are determined by the keywords.  
    | **CalcMod**: Set the parameters for the calculation modes.  
    | **ModPara**: Set the parameters for the basic parameters, such as
      site number, electron number, and Lanczos step.  
    | **LocSpin**: Set the location of the local spin (used only in the
      Kondo model).

(3) Hamiltonian:
    | The Hamiltonian for :math:`{\mathcal H}\Phi` is denoted by the format of the
      interactions for the electron system. The types of interaction are
      determined by the following keywords.  
    | **Trans**: The one-body part,
      :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`  
    | **InterAll**: The general two-body interactions,
      :math:`c_ {i \sigma_1}^{\dagger}c_{j\sigma_2}c_{k \sigma_3}^{\dagger}c_{l \sigma_4}`.
       
    | We can set the interactions that are frequently used by the
      following keywords.  
    | **CoulombIntra**: On-site Coulomb interactions,
      :math:`n_ {i \uparrow}n_{i \downarrow}`
      (:math:`n_{i \sigma}=c_{i\sigma}^{\dagger}c_{i\sigma}`)  
    | **CoulombInter**: Off-site Coulomb interactions,
      :math:`n_ {i}n_{j}` (:math:`n_i=n_{i\uparrow}+n_{i\downarrow}`)  
    | **Hund**: Hund couplings,
      :math:`n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}`
       
    | **PairHop**: Pair hopping couplings,
      :math:`c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{i \downarrow}^{\dagger}c_{j  \downarrow}`
       
    | **Exchange**: Exchange couplings,
      :math:`c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{j \downarrow}^{\dagger}c_{i  \downarrow}`
       
    | **Ising**: Ising interactions, :math:`S_i^z S_j^z`  
    | **PairLift**: PairLift couplings,
      :math:`c_ {i \uparrow}^{\dagger}c_{i\downarrow}c_{j \uparrow}^{\dagger}c_{j \downarrow}`.

(4) Output:
    | The target for the output is determined.  
    | **OneBodyG** : One-body Green’s functions,
      :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}\rangle`  
    | **TwoBodyG** : Two-body Green’s functions,
      :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}c^{\dagger}_{k \sigma_3}c_{l\sigma_4}\rangle`.

      List file for the input files

.. toctree::
   :maxdepth: 1

   List_file_for_the_input_files_en
   CalcMod_file_en
   ModPara_file_en
   LocSpin_file_en
   Trans_file_en
   InterAll_file_en
   CoulombIntra_file_en
   CoulombInter_file_en
   Hund_file_en
   PairHop_file_en
   Exchange_file_en
   Ising_file_en
   PairLift_file_en
   OneBodyG_file_en
   TwoBodyG_file_en
   SingleExcitation_file_en   
   PairExcitation_file_en
   SpectrumVec_File_en
   OneBodyTE_File_en
   TwoBodyTE_File_en
