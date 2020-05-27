U quench in Hubbard model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let's solve the following U quench in 2D Hubbard model at half filling.

.. math::

 H(\tau) = -t \sum_{\langle i,j\rangle , \sigma}(c_{i\sigma}^{\dagger}c_{j\sigma}+{\rm H.c.})
   +U \sum_{i} n_{i\uparrow}n_{i\downarrow}
   +U_{\rm quench} h(\tau) \sum_{i} n_{i\uparrow}n_{i\downarrow}

:math:`h(\tau)` means the step function.
The input files (``samples/tutorial_3.1/stan1.in`` and ``samples/tutorial_3.1/stan2.in``) are as follows ::

 stan1.in

 model = "Hubbard" 
 method = "CG" 
 lattice = "square"
 a0W =  2
 a0L =  2
 a1W =  2
 a1L = -2
 2Sz = 0
 t = 1.0
 U = 4.0
 nelec = 8
 EigenvecIO = "out"

:: 

 stan2.in

 model = "Hubbard"
 method = "Time-Evolution"
 lattice = "square"
 a0W =  2
 a0L =  2
 a1W =  2
 a1L = -2
 2Sz = 0
 t = 1.0
 U = 4.0
 nelec = 8
 PumpType = "Quench"
 Uquench = -8
 EigenvecIO = "in"
 dt = 0.01
 lanczos_max = 1000
 
You can execute HPhi as follows ::

 HPhi -s stan1.in
 HPhi -s stan2.in

Check norm and energy
"""""""""""""""""""""""""""""""
Unitary dynamics of the norm of a wavefunction should be conserved during the real-time evolution.
Using gnuplot, check the dynamics of the norm for this problem ::
  
  plot "output/Norm.dat" u 1:2 w l

In sudden U-quench simulations, the total energy should be conserved for :math:`\tau>0`.
Using gnuplot, check whether the energy is conserved during the real-time evolution ::
  
  plot "output/SS.dat" u 1:2 w l

Dynamics of double occupation
""""""""""""""""""""""""""""""""""
The double occupation :math:`D=\sum_i \langle n_{i\uparrow}n_{i\downarrow} \rangle` is not conserved because :math:`[D, H] \neq 0`.
You can check the dynamics of :math:`D` by executing the following command on gnuplot ::
  
  plot "output/SS.dat" u 1:4 w l
