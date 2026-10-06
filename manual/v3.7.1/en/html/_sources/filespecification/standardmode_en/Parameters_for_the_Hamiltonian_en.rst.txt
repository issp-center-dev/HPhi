.. highlight:: none

Parameters for the Hamiltonian
------------------------------

A default value is :math:`0` unless a specific value is written in the
description. \ :numref:`table_interactions`
shows the parameters for each models. In the case of a complex type, a
file format is “\ *a real part, an imaginary part* " while in the case
of a real type, only “\ *a real part* ".

Local terms
~~~~~~~~~~~

*  ``mu``

   **Type :** Real

   **Description :** It is available only for the Hubbard and Kondo
   lattice model. The chemical potential :math:`\mu` (including the site
   potential) is specified with this parameter.

*  ``U``

   **Type :** Real

   **Description :** It is available only for the Hubbard and Kondo
   lattice model. The onsite Coulomb integral :math:`U` is specified
   with this parameter.

*  ``Jx``, ``Jy``, ``Jz``, ``Jxy``, ``Jyx``, ``Jxz``, ``Jzx``, ``Jyz``,
   ``Jzy``

   **Type :** Real

   **Description :** It is available only for the Kondo lattice model.
   The spin-coupling constant between the valence and the local
   electrons is specified with this parameter. If the exchange coupling
   ``J`` is specified in the input file, instead of ``Jx, Jy, Jz``, the
   diagonal exchange couplings, ``Jx, Jy, Jz``, are set as
   ``Jx = Jy = Jz = J``. When both the set of exchange couplings
   (``Jx``, ``Jy``, ``Jz``) and the exchange coupling ``J`` are
   specified in the input file, :math:`{\mathcal H}\Phi` will stop.

*   ``h``

   **Type :** Real

   **Description :** The longitudinal magnetic field is specified with this parameter.

*   ``Gamma``, ``D``

   **Type :** Real

   **Description :** (Spin model) The transverse magnetic field, and the single-site anisotropy parameter
   are specified with these parameters. The single-site anisotropy
   parameter is not available for ``model=SpinGCCMA``.

The non-local terms described below should be specified differently
according to the lattice structure: For ``lattice=Ladder``, the
non-local terms are specified differently from those for the
``lattice=Chain Lattice``, ``Square Lattice``, ``Triangular Lattice``,
``Honeycomb Lattice``, ``Kagome``. Below, the available parameters for
each lattice are shown in :numref:`table_interactions` .

.. _table_interactions:
.. csv-table:: Interactions for each model defined in an input file. We can define spin couplings as a matrix format.
   :header: "Interactions", "1D chain", "2D square", "2D triangular", "Honeycomb", "Kagome", "Ladder"

   "J, t, V(simplified)", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", "\-"
   "J0, t0, V0", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`"
   "J1, t1, V1", "\-", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`"
   "J2, t2, V2", "\-", "\-", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`"
   "J', t', V'(simplified)", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", "\-"
   "J0', t0', V0'", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", "\-"
   "J1', t1', V1'", "\-", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`"
   "J2', t2', V2'", "\-", "\-", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`"
   "J'', t'', V''(simplified)", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", "\-", "\-"
   "J0'', t0'', V0''", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", "\-", "\-"
   "J1'', t1'', V1''", "\-", ":math:`{\circ}`", ":math:`{\circ}`", ":math:`{\circ}`", "\-", "\-"
   "J2'', t2'', V2''", "\-", "\-", ":math:`{\circ}`", ":math:`{\circ}`", "\-", "\-"


Non-local terms[ for Ladder ( :numref:`fig_ladder` )]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  ``t0``, ``t1``, ``t1'``, ``t2``, ``t2'``

   **Type :** Complex

   **Description :** (Hubbard and Kondo lattice model) Hopping integrals
   in the ladder lattice (see :numref:`fig_ladder` ) are
   specified with this parameter.

*  ``V0``, ``V1``, ``V1'``, ``V2``, ``V2'``

   **Type :** Real

   **Description :** (Hubbard and Kondo lattice model) Offsite Coulomb
   integrals on the ladder lattice ( :numref:`fig_chap04_1_honeycomb` ) 
   are specified with these parameters.

*  ``J0x``, ``J0y``, ``J0z``, ``J0xy``, ``J0yx``, ``J0xz``, ``J0zx``,
   ``J0yz``, ``J0zy``

*  ``J1x``, ``J1y``, ``J1z``, ``J1xy``, ``J1yx``, ``J1xz``, ``J1zx``,
   ``J1yz``, ``J1zy``

*  ``J1'x``, ``J1'y``, ``J1'z``, ``J1'xy``, ``J1'yx``, ``J1'xz``,
   ``J1'zx``, ``J1'yz``, ``J1'zy``

*  ``J2x``, ``J2y``, ``J2z``, ``J2xy``, ``J2yx``, ``J2xz``, ``J2zx``,
   ``J2yz``, ``J2zy``

*  ``J2'x``, ``J2'y``, ``J2'z``, ``J2'xy``, ``J2'yx``, ``J2'xz``,
   ``J2'zx``, ``J2'yz``, ``J2'zy``.

   **Type :** Real

   **Description :** (Spin model) Spin-coupling constants in the ladder
   lattice (see :numref:`fig_ladder` ) are specified with
   these parameters. If the simplified parameter ``J0`` is specified in
   the input file instead of the diagonal couplings, ``J0x, J0y, J0z``,
   these diagonal couplings are set as ``J0x = J0y = J0z = J0``. If both
   J0 and the set of the couplings (J0x, J0y, J0z) are
   specified, :math:`{\mathcal H}\Phi` will stop. The above rules are also valid
   for the simplified parameters, ``J1``, ``J1'``, ``J2``, and ``J2'``.

Non-local terms [other than Ladder ( :numref:`fig_chap04_1_lattice` , :numref:`fig_chap04_1_honeycomb` , :numref:`fig_kagome` )]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  ``t``, ``t0``, ``t1``, ``t2``

   **Type :** Complex

   **Description :** (Hubbard and Kondo lattice model) The nearest
   neighbor hoppings for each direction (see :numref:`fig_chap04_1_lattice` -
   :numref:`fig_kagome` )
   are specified with these parameters. If there is no bond dependence
   of the hoppings, the simplified parameter ``t`` is available to
   specify ``t0``, ``t1``, and ``t2`` as ``t0 = t1 = t2 = t``. If both
   ``t`` and the set of the hoppings (``t0``, ``t1``, ``t2``) are
   specified, :math:`{\mathcal H}\Phi` will stop.

*  ``t'``, ``t0'``, ``t1'``, ``t2'``

   **Type :** Complex

   **Description :** (Hubbard and Kondo lattice model) The second nearest
   neighbor hoppings for each direction (see :numref:`fig_chap04_1_lattice` - :numref:`fig_kagome` )
   are specified with these parameter.
   If there is no bond dependence
   of the hoppings, the simplified parameter ``t'`` is available to
   specify ``t0'``, ``t1'``, and ``t2'`` as ``t0' = t1' = t2' = t'``. If both
   ``t'`` and the set of the hoppings (``t0'``, ``t1'``, ``t2'``) are
   specified, :math:`{\mathcal H}\Phi` will stop.
   
*  ``t''``, ``t0''``, ``t1''``, ``t2''``

   **Type :** Complex

   **Description :** (Hubbard and Kondo lattice model) The third nearest
   neighbor hoppings for each direction (see :numref:`fig_chap04_1_lattice` - :numref:`fig_kagome` )
   are specified with these parameter.
   If there is no bond dependence
   of the hoppings, the simplified parameter ``t''`` is available to
   specify ``t0''``, ``t1''``, and ``t2''`` as ``t0'' = t1'' = t2'' = t''``. If both
   ``t''`` and the set of the hoppings (``t0''``, ``t1''``, ``t2''``) are
   specified, :math:`{\mathcal H}\Phi` will stop.
 
*  ``V``, ``V0``, ``V1``, ``V2``

   **Type :** Real

   **Description :** (Hubbard and Kondo lattice model) The nearest
   neighbor offsite Coulomb integrals :math:`V` for each direction
   (see :numref:`fig_chap04_1_lattice` - :numref:`fig_kagome` )
   are specified with these parameters. If there is no bond dependence
   of the offsite Coulomb integrals, the simplified parameter ``V`` is
   available to specify ``V0``, ``V1``, and ``V2`` as
   ``V0 = V1 = V2 = V``. If both ``V`` and the set of the Coulomb
   integrals (``V0``, ``V1``, ``V2``) are specified, :math:`{\mathcal H}\Phi` will
   stop.

*  ``V'``, ``V0'``, ``V1'``, ``V2'``

   **Type :** Real

   **Description :** (Hubbard and Kondo lattice model) The second nearest
   neighbor-offsite Coulomb integrals :math:`V'` for each direction (see :numref:`fig_chap04_1_lattice` - :numref:`fig_kagome` )
   are specified with this parameter.
   If there is no bond dependence
   of the offsite Coulomb integrals, the simplified parameter ``V'`` is
   available to specify ``V0'``, ``V1'``, and ``V2'`` as
   ``V0' = V1' = V2' = V'``. If both ``V'`` and the set of the Coulomb
   integrals (``V0'``, ``V1'``, ``V2'``) are specified, :math:`{\mathcal H}\Phi` will
   stop.

*  ``V''``, ``V0''``, ``V1''``, ``V2''``

   **Type :** Real

   **Description :** (Hubbard and Kondo lattice model) The third nearest
   neighbor-offsite Coulomb integrals :math:`V'` for each direction (see :numref:`fig_chap04_1_lattice` - :numref:`fig_kagome` )
   are specified with this parameter.
   If there is no bond dependence
   of the offsite Coulomb integrals, the simplified parameter ``V''`` is
   available to specify ``V0''``, ``V1''``, and ``V2''`` as
   ``V0'' = V1'' = V2'' = V''``. If both ``V''`` and the set of the Coulomb
   integrals (``V0''``, ``V1''``, ``V2''``) are specified, :math:`{\mathcal H}\Phi` will
   stop.

*  ``J0x``, ``J0y``, ``J0z``, ``J0xy``, ``J0yx``, ``J0xz``, ``J0zx``,
   ``J0yz``, ``J0zy``

*  ``J1x``, ``J1y``, ``J1z``, ``J1xy``, ``J1yx``, ``J1xz``, ``J1zx``,
   ``J1yz``, ``J1zy``

*  ``J2x``, ``J2y``, ``J2z``, ``J2xy``, ``J2yx``, ``J2xz``, ``J2zx``,
   ``J2yz``, ``J2zy``

   **Type :** Real

   **Description :** (Spin model) The nearest neighbor exchange
   couplings for each direction are specified with these parameters. If
   the simplified parameter ``J0`` is specified, instead of
   ``J0x, J0y, J0z``, the exchange couplings, ``J0x, J0y, J0z``, are set
   as ``J0x = J0y = J0z = J0``. If both ``J0`` and the set of the
   exchange couplings (``J0x, J0y, J0z``) are specified, :math:`{\mathcal H}\Phi`
   will stop. The above rules are valid for ``J1`` and ``J2``.

   If there is no bond dependence of the exchange couplings, the
   simplified parameters, ``Jx``, ``Jy``, ``Jz``, ``Jxy``, ``Jyx``,
   ``Jxz``, ``Jzx``, ``Jyz``, ``Jzy``, are available to specify the
   exchange couplings for every bond as ``J0x = J1x = J2x = Jx``. If any
   simplified parameter (``Jx``-``Jzy``) is specified in addition to its
   counterparts (``J0x``-``J2zy``), :math:`{\mathcal H}\Phi` will stop. Below,
   examples of parameter sets for nearest neighbor exchange couplings
   are shown.

   *  If there are no bond-dependent, and no anisotropic and offdiagonal
      exchange couplings (such as :math:`J_{x y}`), please specify ``J``
      in the input file.

   *  If there are no bond-dependent and offdiagonal exchange couplings
      but there are anisotropic couplings, please specify the non-zero
      couplings in the diagonal parameters, ``Jx, Jy, Jz``.

   *  If there are no bond-dependent exchange couplings but there are
      anisotropic and offdiagonal exchange couplings, please specify the
      non-zero couplings in the nine parameters,
      ``Jx, Jy, Jz, Jxy, Jyz, Jxz, Jyx, Jzy, Jzx``.

   *  If there are no anisotropic and offdiagonal exchange couplings,
      but there are bond-dependent couplings, please specify the
      non-zero couplings in the three parameters, ``J0, J1, J2``.

   *  If there are no anisotropic exchange couplings, but are
      bond-dependent and offdiagonal couplings, please specify the
      non-zero couplings in the nine parameters,
      ``J0x, J0y, J0z, J1x, J1y, J1z, J2x, J2y, J2z``.

   *  If there are bond-dependent, anisotropic, and offdiagonal exchange
      couplings, please specify the non-zero couplings in the
      twenty-seven parameters from ``J0x`` to ``J2zy``.

*  ``J'x``, ``J'y``, ``J'z``, ``J'xy``, ``J'yx``, ``J'xz``, ``J'zx``,
   ``J'yz``, ``J'zy``
*  ``J0'x``, ``J0'y``, ``J0'z``, ``J0'xy``, ``J0'yx``, ``J0'xz``, ``J0'zx``,
   ``J0'yz``, ``J0'zy``
*  ``J1'x``, ``J1'y``, ``J1'z``, ``J1'xy``, ``J1'yx``, ``J1'xz``, ``J1'zx``,
   ``J1'yz``, ``J1'zy``
*  ``J2'x``, ``J2'y``, ``J2'z``, ``J2'xy``, ``J2'yx``, ``J2'xz``, ``J2'zx``,
   ``J2'yz``, ``J2'zy``

   **Type :** Real

   **Description :** (Spin model) The second nearest neighbor exchange
   couplings are specified. However, for ``lattice = Honeycomb Lattice``
   and ``lattice = Kagome`` with ``model=SpinGCCMA``, the second nearest
   neighbor exchange couplings are not available in the :math:`Standard`
   mode. If the simplified parameter ``J'`` is specified, instead of
   ``J'x, J'y, J'z``, the exchange couplings are set as
   ``J'x = J'y = J'z = J'``. If both ``J'`` and the set of the couplings
   (``J'x, J'y, J'z``) are specified, :math:`{\mathcal H}\Phi` will stop.

*  ``J''x``, ``J''y``, ``J''z``, ``J''xy``, ``J''yx``, ``J''xz``, ``J''zx``,
   ``J''yz``, ``J''zy``
*  ``J0''x``, ``J0''y``, ``J0''z``, ``J0''xy``, ``J0''yx``, ``J0''xz``, ``J0''zx``,
   ``J0''yz``, ``J0''zy``
*  ``J1''x``, ``J1''y``, ``J1''z``, ``J1''xy``, ``J1''yx``, ``J1''xz``, ``J1''zx``,
   ``J1''yz``, ``J1''zy``
*  ``J2''x``, ``J2''y``, ``J2''z``, ``J2''xy``, ``J2''yx``, ``J2''xz``, ``J2''zx``,
   ``J2''yz``, ``J2''zy``

   **Type :** Real

   **Description :** (Spin model) The third nearest neighbor exchange
   couplings are specified. However, for ``lattice = Honeycomb Lattice``
   and ``lattice = Kagome`` with ``model=SpinGCCMA``, the third nearest
   neighbor exchange couplings are not available in the :math:`Standard`
   mode. If the simplified parameter ``J''`` is specified, instead of
   ``J''x, J''y, J''z``, the exchange couplings are set as
   ``J''x = J''y = J''z = J''``. If both ``J''`` and the set of the couplings
   (``J''x, J''y, J''z``) are specified, :math:`{\mathcal H}\Phi` will stop.

*  ``phase0``, ``phase1``

   **Type :** Double (``0.0`` as defaults)

   **Description :** We can specify the phase for the hopping through
   the cell boundary with these parameter (unit: degree). These factors
   for the :math:`\boldsymbol{a}_0` direction and the :math:`\boldsymbol{a}_1`
   direction can be specified independently. For the one-dimensional
   system, only ``phase0`` can be used. For example, a fopping from
   :math:`i`-th site to :math:`j`-th site through the cell boundary with
   the positive direction becomes as

   .. math::

          \exp(i \times {\rm phase0}\times\pi/180) \times t {c}_{j \sigma}^\dagger {c}_{i \sigma}
          + \exp(-i \times {\rm phase0}\times\pi/180) \times t^* {c}_{i \sigma}^\dagger {c}_{j \sigma}

.. raw:: latex

   \newpage
