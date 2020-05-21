Spin 1/2 Dimer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let's solve the following spin 1/2 dimer model (2-site Heisenberg model).

.. math::

 H = J {\bf S}_{0}\cdot{\bf S}_{1}

The input file (stan.in) is as follows::

 L=2
 model = "Spin" 
 method = "FullDiag" 
 lattice = "chain"
 J = 0.5
 2Sz = 0
 2S  = 1

It should be noted that the ``L=2`` chain has two sites connected by *two* bonds with each other, and so we let ``J=0.5`` to simulate the dimer model with :math:`J=1`.

You can execute HPhi as follows ::

 HPhi -s stan.in

Check the energy
"""""""""""""""""""""""""""""""
Please check whether the energies are given as follows.

 :math:`E_{\rm min}=-3/4` (singlet) 

 :math:`E_{\rm max}=1/4` (triplet) 

Check S dependence
"""""""""""""""""""""""""""""""
By changing 2S=1 in stan.in, you can treat spin-S dimer (eg. 2S=2 means S=1).
Please check whether the energies are given as follows.

 :math:`E_{\rm min}=-S(S+1)` 

 :math:`E_{\rm max}=S^2` 

Add magnetic field H
"""""""""""""""""""""""""""""""
By adding  H in stan.in, selecting model as "SpinGC", and deleting "2Sz=0" 
you can examine the 
effects of the magnetic field in the dimer model.
An example of the input file (stan.in) is as follows::

 L=2
 model = "SpinGC" 
 method = "FullDiag" 
 lattice = "chain"
 J = 0.5
 2S  = 1
 H   = 2

Please check whether the ground state becomes polarized state (Sz=1).


Try to use Lanczos method
"""""""""""""""""""""""""""""""
**This is just a pedagogical example.**

By selecting method as "Lanczos"
you can perform the Lanczos calculations.
An example of the input file (stan.in) is as follows::

 L=2
 model = "SpinGC" 
 method = "Lanczos" 
 lattice = "chain"
 J = 0.5
 2S  = 1
 H   = 2

This calculation will **fail**! 
Please think why the Lanczos method fails for the dimer. 


Try to use LOBCG method
"""""""""""""""""""""""""""""""
LOBCG is locally optimal block conjugate gradient method.
By selecting method as "CG",
you can perform the LOBCG calculations.
An example of the input file (stan.in) is as follows::

 L=2
 model = "SpinGC" 
 method = "CG" 
 lattice = "chain"
 J = 0.5
 2S  = 1
 H   = 2


In contrast to the Lanczos method, 
this calculation will work well ! 
Please think why the CG method works well for the dimer. 

Please also check the excited states can
be correctly obtained by using LOBCG method.
An example of the input file (stan.in) is as follows::

 L=2
 model = "SpinGC" 
 method = "CG" 
 lattice = "chain"
 J = 0.5
 2S  = 1
 H   = 2
 exct = 4

Here, exct represents the number of excited states, which are
obtained by the LOBCG method.
