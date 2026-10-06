Heisenberg chain (zero temperature)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's solve the following spin 1/2 Heisenberg model on the chain.

.. math::

 H = J \sum_{\langle i,j\rangle}{\bf S}_{i}\cdot{\bf S}_{j}

The input file (``samples/tutorial_1.4/stan1.in``) for 16-site Heisenberg model is as follows::

 L       = 16
 model   = "Spin" 
 method  = "CG" 
 lattice = "chain"
 J = 1
 2Sz = 0
 2S  = 1

You can execute HPhi as follows ::

 HPhi -s stan.in

Check the energy
"""""""""""""""""""""""""""""""
Please check whether the energies are given as follows.

.. math::

 E_{0}= -7.142296 

Obtaining the excited state
"""""""""""""""""""""""""""""""
By adding **exct=2**, you can obtain the 2 low-energy states (``samples/tutorial_1.4/stan2.in``).
Please check the energies.

.. math::
  
 E_{0}= -7.142296

 E_{1}= -6.872107 

Size dependence of the spin gap
"""""""""""""""""""""""""""""""
The spin gap at finite system size is defined
as :math:`\Delta=E_{1}-E_{0}`. For 16-site,
we obtain :math:`\Delta\sim 0.2701`.

Please examine how :math:`\Delta` behaves
as a function of system size L (``samples/tutorial_1.4/stan3.in`` for L=20).
(available system size on PC may be L=24)

Haldane gap
"""""""""""""""""""""""""""""""
By performing a similar calculations for S=1 system,
please examine  how :math:`\Delta` behaves
as a function of system size L (``samples/tutorial_1.4/stan4.in``).
It is known that the finite spin gap exists even
in the thermodynamic limit (:math:`L=\infty`).
This spin gap is often called Haldane gap.
