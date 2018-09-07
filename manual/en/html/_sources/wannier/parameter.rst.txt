Input parameter for Standard mode
=================================

We show the following example of the input file.

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/stan.in

The input parameters for the Standard mode to perform calculation
of the downfolded model are as follows: 
                    
* ``lattice = "wannier90"``

* ``cutoff_t``, ``cutoff_u``, ``cutoff_j``

   **Type :** float

   **Default :** ``1.0e-8``

   The cutoff for the hopping, Coulomb, exchange integrals.
   We ignore these integrals smaller than cutoffs.

*  ``W``, ``L``, ``Height``
*  ``a0W``, ``a0L``, ``a0H``, ``a1W``, ``a1L``, ``a1H``, ``a2W``, ``a2L``, ``a2H``
*  ``Wsub``, ``Lsub``, ``Hsub``
*  ``a0Wsub``, ``a0Lsub``, ``a0Hsub``, ``a1Wsub``, ``a1Lsub``, ``a1Hsub``, ``a2Wsub``, ``a2Lsub``, ``a2Hsub``

   The third dimension appears.
