Hubbard Trimer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let's solve the following the Hubbard trimer model 
(Hubbard model on a triangle).

.. math::

 H = -t \sum_{\sigma}(c_{0\sigma}^{\dagger}c_{1\sigma}+c_{1\sigma}^{\dagger}c_{2\sigma}
   +c_{2\sigma}^{\dagger}c_{0\sigma}+{\rm H.c.})
   +U\sum_{i}(n_{i\uparrow}n_{i\downarrow})

The input file (``samples/tutorial_1.3/stan1.in``) is as follows::

 model = "Hubbard" 
 method = "FullDiag" 
 lattice = "chain" 
 L = 3
 t = -1
 U = 4
 2Sz = 0
 nelec = 2

Note that the filling is not half filling.

You can execute HPhi as follows ::

 HPhi -s stan.in

Ferromagnetic ground state
"""""""""""""""""""""""""""""""
For the Hubbard model on a triangle with one hole, 
it is known that the **perfect ferromagnetism** becomes ground state.
Please check that. 

If you want know the mechanism of the 
ferromagnetism, please see 
**Hal Tasaki, Kotai Butsuri, Vol. 31, 173 (1996)**.
This is one of the simplest example of the 
Nagaoka's ferromagnetism.


Effects of transfer integrals
"""""""""""""""""""""""""""""""
Please the effects of the sign of
the transfer integrals. 
**For example, what happens if you take t = 1 ?**

Another interesting example is by changing 
the transfer integrals between site 0 and site 2.
Following an example of the **trans.def** ::

  ======================== 
  NTransfer      12  
  ======================== 
  ========i_j_s_tijs====== 
  ======================== 
    1     0     0     0         -1.000000000000000         0.000000000000000
    0     0     1     0         -1.000000000000000         0.000000000000000
    1     1     0     1         -1.000000000000000         0.000000000000000
    0     1     1     1         -1.000000000000000         0.000000000000000
    2     0     0     0         -2.000000000000000         0.000000000000000
    0     0     2     0         -2.000000000000000         0.000000000000000
    2     1     0     1         -2.000000000000000         0.000000000000000
    0     1     2     1         -2.000000000000000         0.000000000000000
    2     0     1     0         -1.000000000000000         0.000000000000000
    1     0     2     0         -1.000000000000000         0.000000000000000
    2     1     1     1         -1.000000000000000         0.000000000000000
    1     1     2     1         -1.000000000000000         0.000000000000000

Using this transfer integrals, please examine the
U dependence of the ground state.
Is there phase transition between singlet ground state and
the perfect ferromagnetism ?
