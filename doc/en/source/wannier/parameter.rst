Input parameters for Standard mode
==================================

We show the following example of the input file.

:download:`stan.in <../../../../samples/Wannier/Sr2VO4/stan.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/stan.in

The input parameters for the Standard mode to perform calculation
of the downfolded model are as follows: 

- lattice

  * ``lattice = "wannier90"``

- Parameters related to the lattice

  *  ``W``, ``L``, ``Height``

     **Type :** int

     **Description :** The alignment of original unit cells is specified.
     
  *  ``a0W``, ``a0L``, ``a0H``, ``a1W``, ``a1L``, ``a1H``, ``a2W``, ``a2L``, ``a2H``

     **Type :** int

     **Description :**  Three vectors (:math:`{\vec a}_0, {\vec a}_1, {\vec a}_2`) that specify the lattice . These vectors should be written in the Fractional coordinates of the original transrational vectors.

     
  *  ``Wsub``, ``Lsub``, ``Hsub``

     **Type :** int (positive). In the default setting, ``Wsub=W``, ``Lsub=L``, ``Hsub=Height`` . Namely, there is no sublattice.

     **Description :**
     They are available only in mVMC.
     By using these parameters, we can force the pair-orbital symmetrical to the translation of the sublattice. If the sublattice is incommensurate with the original lattice, ``vmcdry.out`` stops.

     
  *  ``a0Wsub``, ``a0Lsub``, ``a0Hsub``, ``a1Wsub``, ``a1Lsub``, ``a1Hsub``, ``a2Wsub``, ``a2Lsub``, ``a2Hsub``

     **Type :** int (positive). In the default setting, ``a0Wsub=a0W``, ``a0Lsub=a0L``, ``a0Hsub=a0H``,
     ``a1Wsub=a1W``, ``a1Lsub=a1L``, ``a1Hsub=a1H``, ``a2Wsub=a2W``, ``a2Lsub=a2L``, ``a2Hsub=a2H``. Namely, there is no sublattice.

     **Description :**
     The manner to set these aparameters is same as that for ``a0W``, ``a0L``, ``a0H``,
     ``a1W``, ``a1L``, ``a1H``,   ``a2W``, ``a2L``, ``a2H``.
     If the sublattice is incommensurate with the original lattice, ``vmcdry.out`` stops.


- Parameters related to interactions
     
  * ``lambda_u``

    **Type :** float (greater than or equal to 0)

    **Default :** ``1.0``

     **Description :** A parameter to tune the strength of Coulomb interactions by multiplying :math: `lambda_u` by them.

  * ``lambda_j``

    **Type :** float (greater than or equal to 0)

    **Default :** ``1.0``

     **Description :** A parameter to tune the strength of exchange Coulomb interactions by multiplying :math: `lambda_j` by them.
    
  
  * ``lambda``

    **Type :** float (greater than or equal to 0)

    **Default :** ``1.0``

     **Description :** A parameter to tune the strength of Coulomb and exchange interactions by multiplying :math: `lambda` by them. When  :math:`\lambda_U` , :math:`\lambda_J` are specified, these settings are used.

  * ``cutoff_t``, ``cutoff_u``, ``cutoff_j``

    **Type :** float

    **Default :** ``1.0e-8``

    **Description :** The cutoff parameters for the hopping, Coulomb, exchange integrals.
    We ignore these integrals smaller than cutoff values.

  * ``cutoff_tW``, ``cutoff_tL``, ``cutoff_tH``
  * ``cutoff_UW``, ``cutoff_UL``, ``cutoff_UH``
  * ``cutoff_JW``, ``cutoff_JL``, ``cutoff_JH``

    **Type :** 

    **Default :**  ``cutoff_tW = int((W-1)/2)``, ``cutoff_tL=int((L-1)/2)``, ``cutoff_tH=int((Height-1)/2)`` (when ``W`` , ``L`` and ``Height`` are not defined, the values are set to 0) and others are set to ``0``. 

    **Description :** The cutoff parameters for the hopping, Coulomb, exchange integrals.
    We ignore these integrals that have lattice vector :math:`{\bf R}` larger than these values.

  * ``cutoff_length_t``, ``cutoff_length_U``, ``cutoff_length_J``

    **Type :** float

    **Default :** ``cutoff_length_t = -1.0`` (include all terms), others are set to ``0.3``.

    **Description**

    The cutoff parameters for the hopping, Coulomb, exchange integrals.
    We ignore these integrals whose distances are longer than these values.
    The distances are computed from the position of the Wannier center and  unit lattice vectors.

- Parameters for one body correction

  To avoid double countings in analyzing the lattice model,
  one body correction is done by subtracting the following terms from one body terms: 
  
       .. math::
	  \begin{aligned}
	  t_{mm}^{\rm DC}({\bf 0}) &\equiv \alpha U_{mm}({\bf 0}) D_{mm}({\bf 0})
	  + \sum_{({\bf R}, n) \neq ({\bf 0}, m)} U_{m n} ({\bf R})D_{nn}({\bf 0})\\
	  & - (1-\alpha) \sum_{({\bf R}, n) \neq ({\bf 0}, 0)} J_{m n}({\bf R}) D_{nn}({\bf R}),\\
	  t_{mn}^{\rm DC}({\bf R}_{ij}) &\equiv \frac{1}{2} J_{mn}({\bf R}_{ij}) \left(D_{nm}({\bf R}_{ji}) + 2 {\rm Re} [D_{nm}({\bf R}_{ji})]\right)\\
	  &-\frac{1}{2}  U_{mn}({\bf R}_{ij}) D_{nm}({\bf R}_{ji}),
	  \quad ({\bf R}_{ij}, m) \neq ({\bf 0}, n),
	  \\
	  D_{mn}({\bf R}_{ij}) &\equiv \sum_{\sigma}
	  \left\langle c_{im \sigma}^{\dagger} c_{jn \sigma}\right\rangle_{\rm KS},
	  \end{aligned}

  where, the first and second terms correspond to the Hartree and Fock corrections, respectively.
  :math:`\alpha` is a tuning parameter for one body correction from the on-site Coulomb interactions.

  * ``doublecounting``

     **Type :** char

     **Default :** ``none``

     **Description :** 
     
     ``none``: One body correction is not considered.  ``Hartree_U``: Hartree correction only considered the contribution from Coulomb interactions :math:`U_{Rii}` . ``Hartree``: Hartree correction. ``full``: One body correction including Fock correction. 
     Charge densities :math:`D_{Rij}` are obtained by ``[CDataFileHead]_dr.dat`` which is automatically outputted by RESPACK.
     It is noted that the charge densities are assumed not to depend on spin components.
     
  * ``alpha``

     **Type :** float

     **Default :** ``0.5``

     **Description :**

     A tuning parameter for one body correction from the on-site Coulomb interactions (:math:`0\le \alpha \le 1`).
