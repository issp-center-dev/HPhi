Tutorials for real-time evolution
=================================

U quench in Hubbard model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let's solve the following U quench in 2D Hubbard model at half filling.

.. math::

 H(\tau) = -t \sum_{\langle i,j\rangle , \sigma}(c_{i\sigma}^{\dagger}c_{j\sigma}+{\rm H.c.})
   +U \sum_{i} n_{i\uparrow}n_{i\downarrow}
   +U_{\rm quench} h(\tau) \sum_{i} n_{i\uparrow}n_{i\downarrow}

:math:`h(\tau)` means the step function.
The input files (stan1.in and stan2.in) are as follows ::

 stan1.in

 model = "Hubbard" 
 method = "CG" 
 lattice = "square"
 a0W =  2
 a0L =  2
 a1W =  2
 a1L = -2
 J = 0.5
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


Dynamical phase transition in 1D transverse-field Ising model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Dynamical phase transition (DPT) is a phase transition due to nonequilibrium process such as sudden quench.
In this section, we introduce the DPT in 1D transverse Ising model, which is one of the famous examples for the DPT.
The details are found in this paper (M. Heyl, A. Polkovnikov, and S. Kehrein, Phys. Rev. Lett. **110**, 135704 (2013)).

.. math::

 H = J \sum_{\langle i,j\rangle} S^z_{i} S^z_{j} + \Gamma \sum_{i} S^x_i 

Ground state
"""""""""""""
First, you need to obtain the initial state for simulations of the real-time dynamics.
Please make the following input file (stan1.in) ::

 model = "SpinGC"
 method = "CG"
 lattice = "chain"
 L = 12
 Jz = -1.0
 Gamma = 0.1
 h = 1e-5
 EigenvecIO = "out"

and run the following command ::
 
  HPhi -s stan1.in

Check the total magnetization along :math:`z`-direction :math:`M_z = \sum_i \langle S^z_i \rangle` in "output/zvo_energy.dat". 

Question: If a longitudinal magnetic field :math:`h` in stan1.in becomes 0, what happens?

:math:`\Gamma` quench
""""""""""""""""""""""
After obtaining the ground state, you can perform simulations for :math:`\Gamma` quench in the 1D transverse-field Ising model.
Make the following input file (stan2.in) ::

 model = "SpinGC" 
 method = "Time-Evolution"
 lattice = "chain"
 L = 12
 Jz = -1.0
 Gamma = 0.2 
 h = 1e-5
 EigenvecIO = "in"
 dt = 0.01
 lanczos_max = 1000

In this case, :math:`\Gamma` quench from 0.1 to 0.2 will be performed by executing the following command ::

  HPhi -s stan2.in

Now you can check the real-time evolution of :math:`M_z` in "output/Flct.dat".
The result is plotted by executing the following command on gnuplot ::

 gnuplot
 gnuplot> set xlabel "time t"
 gnuplot> set ylabel "Magnetization per site |M_z/L|"
 gnuplot> p "output/Flct.dat" u 1:(abs($6/12)) w l tit "Gamma=0.1 -> 0.2"

One of the features of DPT in this model is that cusp structures in the dynamics of :math:`M_z` appears at the same intervals.
The following figure is an example for several results for :math:`\Gamma` quench.
You can see the cusp structure for :math:`\Gamma > 0.5`.

.. image:: ../../figs/Mz_DPT_1D_TIM.pdf
   :height: 500px
   :width: 700px
   :align: center



