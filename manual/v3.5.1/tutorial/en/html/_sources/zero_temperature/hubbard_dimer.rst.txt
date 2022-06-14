Hubbard Dimer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let's solve the following the Hubbard dimer model.

.. math::

 H = -t \sum_{\sigma}(c_{0\sigma}^{\dagger}c_{1\sigma}+{\rm H.c.})
   +U(n_{0\uparrow}n_{0\downarrow}+n_{1\uparrow}n_{1\downarrow})

The input file (``samples/tutorial_1.2/stan1.in``) is as follows::

 model = "Hubbard" 
 method = "FullDiag" 
 lattice = "chain" 
 L=2
 t = -0.5 
 U = 4
 2Sz = 0
 nelec = 2

You can execute HPhi as follows ::

 HPhi -s stan.in

Check the energy
"""""""""""""""""""""""""""""""
For the Hubbard dimer at half filling with total Sz=0, 
energies are given as follows:

 :math:`E=0,U,\frac{U}{2}\times(1\pm\sqrt{(1+(4t/U)^2)})` 

For example, by taking :math:`U=4,t=-1`, the 
energies are  given as follows:

 :math:`E=-0.828427, 0, 4, 4.828427` 

It is note that simple mathematical calculations 
can be done using:: 

 bc -l 

on the terminal.

Try to use LOBCG method
"""""""""""""""""""""""""""""""
The input file (``samples/tutorial_1.2/stan2.in``) is as follows::

 model = "Hubbard" 
 method = "CG" 
 lattice = "chain" 
 L=2
 t = -0.5 
 U = 4
 2Sz = 0
 nelec = 2
 exct = 4

Please check whether LOBCG method correctly 
reproduces the energies including the excited states.
