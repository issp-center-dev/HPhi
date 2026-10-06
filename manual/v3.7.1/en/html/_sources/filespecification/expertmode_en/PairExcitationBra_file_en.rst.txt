.. highlight:: none

.. _Subsec:pairexcitationbra:

PairExcitationBra file
----------------------

The operators :math:`B` to generate the **bra** pair excited state
:math:`B|\phi\rangle`, used to compute the *off-diagonal* dynamical Green's
function

.. math::

   G_{BA}(z) = \langle \phi | B^{\dagger} \frac{1}{z-H} A | \phi \rangle ,

where the ket operator :math:`A` is given in
:doc:`PairExcitation <PairExcitation_file_en>`. If this file is omitted, the
ordinary diagonal Green's function (:math:`B=A`) is calculated.

The file format and parameters are the same as the
:doc:`PairExcitation <PairExcitation_file_en>` file; only the keyword for the
number of operators is ``NPairExcitationBra``. An example
(a density operator :math:`n_{1\uparrow}` on site 1)::

    ===============================
    NPairExcitationBra    1
    ===============================
    ======= Pair Excitation Bra ====
    ===============================
        1     0     1     0     1    1.0    0.0

.. _use_rules_pairexcitationbra:

Use rules
~~~~~~~~~

*  Headers cannot be omitted; the format and parameters are identical to the
   :doc:`PairExcitation <PairExcitation_file_en>` file.

*  The operator written here is :math:`B`, **not** :math:`B^{\dagger}`. The
   computed quantity is :math:`\langle\phi|B^{\dagger}(z-H)^{-1}A|\phi\rangle`.

*  The off-diagonal Green's function is available only with the shifted BiCG
   method (``method="CG"``) and ``CalcSpec="Normal"`` (restart/save spectrum
   modes are not supported together with a bra operator).

*  The ket (:doc:`PairExcitation <PairExcitation_file_en>`) and bra operators
   must map :math:`|\phi\rangle` into the **same** excited Hilbert sector (the
   same change of particle number and :math:`S_z`). All operators within one
   file must also share the same sector change. Otherwise the program is
   terminated.

*  Supported models: Hubbard, Spin, and SpinGC. (Other models are not yet
   supported for the off-diagonal Green's function; the program is terminated for
   them.) ``SingleExcitationBra`` and ``PairExcitationBra`` cannot be specified
   at the same time.

.. raw:: latex

   \newpage
