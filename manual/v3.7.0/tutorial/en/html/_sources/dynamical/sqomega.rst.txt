Spectrum calculation for 12-site one-dimeinsional Heisenberg chain model.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Hamiltonian is given by

.. math::

 H = J \sum_{i=0}^{11}{\bf S}_{i}\cdot{\bf S}_{i+1},

where

.. math::
   
 {\bf S}_0 = {\bf S}_{12}.
 

The spectrum function can be calculated by following steps.

1. Calculate the ground state.
2. Define excitation operators in the ``pair.def`` file.
3. Calculate spectrum function.

See the manual for details. To simply do above steps, we prepare the
script file ``samples/tutorial_4.3/spinchain_example.py`` 
.. in https://github.com/issp-center-dev/HPhi-gallery/tree/master/Spin/HeisenbergSpectrum. 
In the following, we show the
procedure to obtain the specrum function by using the script file.

1. Execute the script file (``spinchain_example.py``)

   ::

      $ python spinchain_example.py


2. Plot ``spectrum.dat`` by gnuplot.

   ::

      $ gnuplot
      $ set yrange [0:5]
      $ set pm3d map
      $ splot "./spectrum.dat" using 1:2:3

   You can see the following figure, where horizontal and vertical
   axises correspond to the index of wave vector and frequency,
   respectively.

.. image:: ../../../figs/sqomega.*
   :height: 500px
   :width: 700px
   :align: center

