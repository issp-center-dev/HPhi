.. highlight:: none

Parameters for conserved quantities
-----------------------------------

*  ``nelec``

   **Type :** Positive integer

   **Description :** The number of valence electrons is specified with
   this parameter. When model = ``"Fermion HubbardGC"``, ``"Spin"``, or
   ``"SpinGC"``, it should not be specified.

*  ``2Sz``

   **Type :** Integer

   **Description :** The :math:`z` component of the twofold total spin
   is specified with this parameter. When
   model = ``"Fermion HubbardGC"`` or ``"SpinGC"``, it should not be
   specified.

*  ``MomentumIndex``

   **Type :** Integer

   **Description :** The one-dimensional translation momentum sector is
   specified by an integer :math:`m`. The generated expert input contains
   ``qptransidx.def`` and the corresponding ``TransSym`` entry in
   ``namelist.def``. For a chain of length :math:`L`, the character of
   the translation by :math:`g` sites is
   :math:`\exp(-2\pi i m g/L)`, and ``MomentumIndex`` must satisfy
   :math:`0 \le m < L`.

   This parameter is currently supported only for HPhi Standard mode with
   ``model = "Spin"``, ``lattice = "chain"``, ``2S = 1``, fixed
   ``2Sz``, periodic boundary conditions without ``phase0``, and
   exchange-only spin couplings with ``Jz = 0`` and no field, general, or
   pair terms. The supported calculation methods are ``"Lanczos"`` and
   ``"CG"``.

.. raw:: latex

   \newpage
