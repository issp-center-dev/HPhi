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
   ``lattice = "chain"``, periodic boundary conditions without ``phase0``,
   and the ``"Lanczos"`` or ``"CG"`` calculation method. The supported
   models are:

   * ``model = "Spin"`` with ``2S = 1``, fixed ``2Sz``, and
     exchange-only spin couplings with ``Jz = 0`` and no field, general,
     or pair terms.
   * ``model = "SpinlessFermion"`` with fixed ``ncond`` and hopping-only
     terms. Density interactions (``V``/``CoulombInter``), general
     interactions, and pair terms are not supported with
     ``MomentumIndex`` in Standard mode.
   * ``model = "Fermion Hubbard"`` (or ``"Hubbard"``) with fixed
     ``nelec`` and ``2Sz``, nearest-neighbor hopping ``t``/``t0``, and
     onsite ``U``. Chemical-potential or field terms
     (``mu``, ``h``, ``Gamma``, ``Gamma_y``), longer-range hopping
     (``t'``, ``t''``), offsite Coulomb terms (``V``/``CoulombInter``),
     and general or pair terms are not supported with ``MomentumIndex``.

   To obtain the lowest energy at every momentum point, run one Standard-mode
   calculation for each :math:`m=0,\ldots,L-1`, changing only
   ``MomentumIndex``. For example, a six-site spin chain in the :math:`m=1`
   sector can be specified as follows:

   .. code-block:: text

      L = 6
      model = Spin
      method = Lanczos
      lattice = chain
      Jx = 1.0
      Jy = 1.0
      Jz = 0.0
      2Sz = 0
      MomentumIndex = 1

   If this input is saved as ``stan.in``, run ``HPhi -s stan.in``. The
   corresponding momentum is :math:`k=2\pi m/L`, and the lowest energy is
   written to ``output/zvo_energy.dat``. Run each value of ``MomentumIndex``
   in a separate working directory, or save the output before the next run,
   so that the results are not overwritten.

.. raw:: latex

   \newpage
