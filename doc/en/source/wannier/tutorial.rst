.. _tutorialwannier:

Tutorial
========

In this tutorial, we downfold Sr\ :sub:`2`\ VO\ :sub:`4`
into three-orbitals 2D Hubbard model,
and simulate that model with HPhi/mVMC.
We employ QuantumESPRESSO for the DFT calculation.

SCF calculation of charge density
---------------------------------

First, we perform the SCF calculation of the charge density.
The input file is as follows:

:download:`scf.in <../../../../samples/Wannier/Sr2VO4/scf.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/scf.in

The pseudopotential (UPF file) are downloaded from
`The SG15 Optimized Norm-Conserving Vanderbilt (ONCV) pseudopotentials <www.quantum-simulation.org/potentials/sg15_oncv/>`_.

http://www.quantum-simulation.org/potentials/sg15_oncv/sg15_oncv_upf_2015-10-07.tar.gz

We use the program ``pw.x`` in QuantumESPRESSO as follows.

.. code-block:: bash

   $ pw.x -in scf.in

(Optional) Band structure
-------------------------

:download:`band.in <../../../../samples/Wannier/Sr2VO4/band.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/band.in

We use ``pw.x``.
                    
.. code-block:: bash

   $ pw.x -in band.in

:download:`bands.in <../../../../samples/Wannier/Sr2VO4/bands.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/bands.in

We use ``bands.x`` QuantumESPRESSO.
                    
.. code-block:: bash

   $ bands.x -in bands.in

We can plot the band structure by reading output ``bands.out.gnu`` from
GnuPlot etc.
   
Kohn-Sham orbitals for Wannier
------------------------------

:download:`nscf.in <../../../../samples/Wannier/Sr2VO4/nscf.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/nscf.in

We use ``pw.x`` as
                    
.. code-block:: bash

   $ pw.x -in nscf.in

Then, we use the utility ``qe2respack.sh`` which is included in the RESPACK package.
The command-line argument is the name of ``[prefix].save`` directory.

.. code-block:: bash

   $ qe2respack.sh sr2cuo3.save
                
Wannier function, dielectric function, effective interaction
------------------------------------------------------------

:download:`respack.in <../../../../samples/Wannier/Sr2VO4/respack.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/respack.in

We use ``calc_wannier``, ``calc_chiqw``, ``calc_j3d``,
``calc_w3d`` in RESPACK.
                    
.. code-block:: bash

   $ calc_wannier < respack.in
   $ calc_chiqw < respack.in
   $ calc_w3d < respack.in
   $ calc_j3d < respack.in

Quantum lattice mode for HPhi/mVMC
----------------------------------

First, we translate the hopping file etc. created by RESPACK
into the Wannier90 format.
For this purpose, we use the utility ``respack2wan90.py`` included in HPhi/mVMC.
The command-line argument should be the same as the input parameter ``CDataFileHead``
for Standard mode of HPhi/mVMC.
If we do not specify that argument, the default value ``zvo`` is used.

.. code-block:: bash

   $ respack2wan90.py zvo

Then we can run HPhi/mVMC with the standard mode.
                    
:download:`respack.in <../../../../samples/Wannier/Sr2VO4/stan.in>`

.. literalinclude:: ../../../../samples/Wannier/Sr2VO4/stan.in

.. code-block:: bash

   $ vmc.out -s stan.in
