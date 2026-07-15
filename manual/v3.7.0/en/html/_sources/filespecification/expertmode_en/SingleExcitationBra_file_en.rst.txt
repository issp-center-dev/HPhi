.. highlight:: none

.. _Subsec:singleexcitationbra:

SingleExcitationBra file
------------------------

The operators :math:`B=c_{i\sigma_1}(c_{i\sigma_1}^{\dagger})` to generate the
**bra** single excited state :math:`B|\phi\rangle`, used to compute the
*off-diagonal* dynamical Green's function

.. math::

   G_{BA}(z) = \langle \phi | B^{\dagger} \frac{1}{z-H} A | \phi \rangle ,

where the ket operator :math:`A` is given in
:doc:`SingleExcitation <SingleExcitation_file_en>`. If this file is omitted, the
ordinary diagonal Green's function (:math:`B=A`) is calculated.

The file format and parameters are the same as the
:doc:`SingleExcitation <SingleExcitation_file_en>` file; only the keyword for the
number of operators is ``NSingleExcitationBra``. An example::

    ===============================
    NSingleExcitationBra    1
    ===============================
    ====== Single Excitation Bra ===
    ===============================
        0     0     0    1.0    0.0

.. _use_rules_singleexcitationbra:

Use rules
~~~~~~~~~

*  Headers cannot be omitted; the format and parameters are identical to the
   :doc:`SingleExcitation <SingleExcitation_file_en>` file.

*  The operator written here is :math:`B`, **not** :math:`B^{\dagger}`. The
   computed quantity is :math:`\langle\phi|B^{\dagger}(z-H)^{-1}A|\phi\rangle`.

*  The off-diagonal Green's function is available only with the shifted BiCG
   method (``method="CG"``) and ``CalcSpec="Normal"`` (restart/save spectrum
   modes are not supported together with a bra operator).

*  The ket (:doc:`SingleExcitation <SingleExcitation_file_en>`) and bra operators
   must map :math:`|\phi\rangle` into the **same** excited Hilbert sector, and all
   operators within one file must share the same sector change; otherwise the
   program is terminated. What "same sector" means depends on the model:

   *  Hubbard (canonical): the same change of particle number **and**
      :math:`S_z`.

   *  HubbardGC (grand canonical): a single excitation carries no sector shift,
      so the whole Fock space is one Hilbert space and the ket and bra (of either
      spin) always share it.

   *  HubbardNConserved (:math:`S_z` free): only the same change of particle
      number is required; the :math:`S_z` change may differ.

*  Supported models: Hubbard, HubbardGC (grand canonical), and
   HubbardNConserved (canonical, fixed electron number with free :math:`S_z`).
   (A single excitation is not defined for spin systems, and other models are
   not yet supported for the off-diagonal Green's function; the program is
   terminated for them.)

*  HubbardNConserved is selected from the Hubbard model (``CalcModel=0``) by
   giving ``Ncond`` (the electron number) without ``2Sz`` in the
   :ref:`ModPara <Subsec:modpara>` file. Because :math:`S_z` is then free, this
   enables the cross-spin (2Sz-free) off-diagonal route (used, for example, for
   spin-orbit-coupled impurities): the excited Hilbert
   sector is fixed by the electron number alone, so the ket and bra need only
   share the same change of particle number, and a cross-spin :math:`G_{BA}` is
   obtained from a single solve.

*  When ``SpectrumNumBra`` :math:`>1` is set in the
   :ref:`ModPara <Subsec:modpara>` file, the additional bra operator sets
   ``1`` ... :math:`(B-1)` are read from files named ``single_ex_bra_1.def`` ...
   ``single_ex_bra_``\ :math:`(B-1)`\ ``.def`` in the run directory, using the
   same format as this file; bra set ``0`` is this ``SingleExcitationBra`` file
   given in the namelist.

.. raw:: latex

   \newpage
