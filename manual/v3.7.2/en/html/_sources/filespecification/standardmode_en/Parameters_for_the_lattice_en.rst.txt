.. highlight:: none

Parameters for the lattice
--------------------------

Chain [ :numref:`fig_chap04_1_lattice` (a)]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  ``L``

   **Type :** Integer

   **Description :** The length of the chain is specified with this
   parameter.

   .. figure:: ../../../../figs/chap04_1_lattice.png
      :name: fig_chap04_1_lattice
      :alt: Schematic illustration of (a) one-dimensional chain lattice,
            (b) two-dimensional square lattice, and (c) two-dimensional
            triangular lattice. They have :math:`t`, :math:`V`, and :math:`J`
            as the nearest neighbor hopping, an offsite Coulomb integral, and
            a spin-coupling constant, respectively (magenta solid lines); they
            also have :math:`t'`, :math:`V'`, and :math:`J'` as the next
            nearest neighbor hopping, offsite Coulomb integral, and
            spin-coupling constant, respectively (green dashed line).
            
      Schematic illustration of (a) one-dimensional chain lattice, (b)
      two-dimensional square lattice, and (c) two-dimensional triangular
      lattice. They have :math:`t`, :math:`V`, and :math:`J` as the
      nearest neighbor hopping, an offsite Coulomb integral, and a
      spin-coupling constant, respectively (magenta solid lines); they
      also have :math:`t'`, :math:`V'`, and :math:`J'` as the next
      nearest neighbor hopping, offsite Coulomb integral, and
      spin-coupling constant, respectively (green dashed line). 

   .. figure:: ../../../../figs/chap04_1_honeycomb.png
      :name: fig_chap04_1_honeycomb
      :alt: Schematic illustration of the anisotropic honeycomb lattice.
            The first/second/third nearest neighbor hopping integral,
            spin coupling, and offsite
            Coulomb integral depend on the bond direction.
            
      Schematic illustration of the anisotropic honeycomb lattice.
      The first/second/third nearest neighbor hopping integral,
      spin coupling, and offsite
      Coulomb integral depend on the bond direction.

   .. figure:: ../../../../figs/kagome.png
      :name: fig_kagome
      :alt: Schematic illustration of the Kagome lattice.
      
      Schematic illustration of the Kagome lattice. 

   .. figure:: ../../../../figs/ladder.png
      :name: fig_ladder
      :alt: Schematic illustration of the ladder lattice.
      
      Schematic illustration of the ladder lattice. 

Ladder ( :numref:`fig_ladder` )
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  ``L``

   **Type :** Integer

   **Description :** The length of the ladder is specified with this
   parameter.

*  ``W``

   **Type :** Integer

   **Description :** The number of the ladder is specified with this
   parameter.

  .. figure:: ../../../../figs/chap04_1_unitlattice.png
     :name: fig_chap04_1_unitlattice
     :scale: 100%
     :alt: Shape of the numerical cell when
           :math:`{\boldsymbol a}_0 = (6, 2), {\boldsymbol a}_1 = (2, 4)` in the triangular
           lattice. The region surrounded by :math:`{\boldsymbol a}_0` (magenta dashed
           arrow) and :math:`{\boldsymbol a}_1` (green dashed arrow) becomes the cell
           to be calculated (20 sites).
     
     Shape of the numerical cell when
     :math:`{\boldsymbol a}_0 = (6, 2), {\boldsymbol a}_1 = (2, 4)` in the triangular
     lattice. The region surrounded by :math:`{\boldsymbol a}_0` (magenta dashed
     arrow) and :math:`{\boldsymbol a}_1` (green dashed arrow) becomes the cell
     to be calculated (20 sites). 

**Tetragonal lattice** [ :numref:`fig_chap04_1_lattice` (b)], triangular lattice [ :numref:`fig_chap04_1_lattice` (c)], 
honeycomb lattice [ :numref:`fig_chap04_1_honeycomb` ], Kagome lattice [ :numref:`fig_kagome` ]

In these lattices, we can specify the shape of the numerical cell by
using the following two methods.

*   ``W``, ``L``

   **Type :** Integer

   **Description :** The alignment of the original unit cells (dashed
   black lines in :numref:`fig_chap04_1_lattice`  - :numref:`fig_kagome` ) is specified with this parameter.

*  ``a0W``, ``a0L``, ``a1W``, ``a1L``

   **Type :** Integer

   **Description :** We can specify two vectors
   (:math:`{\boldsymbol a}_0, {\boldsymbol a}_1`) that surround the numerical cell
   (:numref:`fig_chap04_1_unitlattice` ).
   These vectors should be specified in the fractional coordinate.

If we use both these methods, :math:`{\mathcal H}\Phi` stops. When
``model=SpinGCCMA``, we can use only the former.

We can check the shape of the numerical cell by using a file
``lattice.gp`` which is written in Standard mode. This file can be read
by ``gnuplot`` as follows:

::

    $ gnuplot lattice.gp

.. raw:: latex

   \newpage
