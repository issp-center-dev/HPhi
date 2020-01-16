.. _tutorialwannier:

Tutorial
========

In this tutorial, we downfold Sr\ :sub:`2`\ CuO\ :sub:`3`
into single orbital 1D Hubbard model,
and simulate that model with HPhi/mVMC.
We use QuantumESPRESSO for the DFT calculation.
Input files are served in ``samples/Wannier/Sr2CuO3`` directory.


In actual studies, the input files etc. of each solver should be modified for more high accuracy calculation.
Please refer to the manuals of each solver for the details of the input files.

SCF calculation of charge density
---------------------------------

First, we perform the SCF calculation of the charge density.
The input file is as follows:

:download:`scf.in <../../../../samples/Wannier/Sr2CuO3/scf.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2CuO3/scf.in

The pseudopotential (UPF file) are downloaded from
`The SG15 Optimized Norm-Conserving Vanderbilt (ONCV) pseudopotentials <www.quantum-simulation.org/potentials/sg15_oncv/>`_.
Put the pseudopotential files into ``../pseudo`` directory.

http://www.quantum-simulation.org/potentials/sg15_oncv/sg15_oncv_upf_2015-10-07.tar.gz

We use the program ``pw.x`` in QuantumESPRESSO as follows.

.. code-block:: bash

   $ pw.x -in scf.in

(Optional) Band structure
-------------------------

:download:`band.in <../../../../samples/Wannier/Sr2CuO3/band.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2CuO3/band.in

We use ``pw.x``.
                    
.. code-block:: bash

   $ pw.x -in band.in

:download:`bands.in <../../../../samples/Wannier/Sr2CuO3/bands.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2CuO3/bands.in

We use ``bands.x`` QuantumESPRESSO.
                    
.. code-block:: bash

   $ bands.x -in bands.in

We can plot the band structure by reading output ``bands.out.gnu`` from
GnuPlot etc.
   
Kohn-Sham orbitals for Wannier
------------------------------

:download:`nscf.in <../../../../samples/Wannier/Sr2CuO3/nscf.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2CuO3/nscf.in

We use ``pw.x`` as
                    
.. code-block:: bash

   $ pw.x -in nscf.in

Then, we use the utility ``qe2respack.py`` which is included in the RESPACK package.
The command-line argument is the name of ``[prefix].save`` directory.

.. code-block:: bash

   $ qe2respack.py sr2cuo3.save
                
Wannier function, dielectric function, effective interaction
------------------------------------------------------------

:download:`respack.in <../../../../samples/Wannier/Sr2CuO3/respack.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2CuO3/respack.in

We use ``calc_wannier``, ``calc_chiqw``, ``calc_j3d``,
``calc_w3d`` in RESPACK.
                    
.. code-block:: bash

   $ calc_wannier < respack.in
   $ calc_chiqw < respack.in
   $ calc_w3d < respack.in
   $ calc_j3d < respack.in

After finishing calculations, the files are outputted in ``dir-model`` folder. 
The format of these files is Wannier90 format and the data such as the hopping integrals are written.
(If you use the old version of RESPACK (20190226), the folder name is  ``dir-mvmc`` .)

   
Quantum lattice mode for HPhi/mVMC
----------------------------------

Using standard mode of HPhi/mVMC, the calculation will be done by reading the files in ``dir-model`` folder.
First, the files in ``dir-model`` directory should be moved to the current directry.
Then, the calculation will be started by using standard mode.
For example, in HPhi, the calculation will be dobe by typing the following command:
                    
:download:`stan.in <../../../../samples/Wannier/Sr2CuO3/stan.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2CuO3/stan.in

.. code-block:: bash

   $ cp ./dir-model/* .
   $ HPhi -s stan.in
   
