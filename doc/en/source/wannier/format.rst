File format
===========

Standard mode of HPhi/mVMC reads the following files.
We can obtain these files by executing RESPACK.
The files are outputted in the ``dir-model`` directory.

Geometry
--------

The file name is ``[CDataFileHead]_geom.dat``.
By editing this file, we can modify the number of orbitals treated in HPhi/mVMC.

::

   -1.917800 1.917800 6.280100
   1.917800 -1.917800 6.280100
   1.917800 1.917800 -6.280100
   3
   0.000000 -0.000000 -0.000000
   -0.000000 -0.000000 -0.000000
   0.000000 0.000000 0.000000

* Lines 1 - 3

  Unit lattice vectors in the Cartesian coordinate (arbitrary unit).

* Line 4

  The number of orbitals par unit cell treated by mVMC/HPhi.
  When we reduce the number by editing this file,
  the model including the same number of orbitals from the top.

* Line 5 - end

  Wannier centers in the fractional coordinate. They are used by the Fourier utility.
  The order in which wannier functions are defined corresponds to the index of wannier functions 
  in the following four files.
  
Hopping, Coulomb, exchange integrals, charge densities
-------------------------------------------------------

The file names are ``[CDataFileHead]_hr.dat``, ``[CDataFileHead]_ur.dat``,
, ``[CDataFileHead]_jr.dat`` and ``[CDataFileHead]_dr.dat``, respectively.
They are formatted as the hopping-integral file of Wannier90 is used.
See details in the ``8.19 seedname_hr.dat`` in the user_guide for wannier90.
	
	::
	
	  wannier90 format for mvmcdry
	  8
	  343
	  1    1    1    1    1    1    1    1    1    1    1    1    1
	  ...
	   -3   -3   -3    1    1   0.0004104251  -0.0000000000
	   -3   -3   -3    1    2   0.0001515941  -0.0000000006
	   -3   -3   -3    1    3  -0.0001515941   0.0000000002
	
	* Line 1
	
	  File Header
	
	* Line 2
	
	  Total number of wannier functions.
	
	* Line 3
	
	  Total number of super cells ``nrpts`` .
	
	* Line 4- Line 5 + int(``nrpts``/15) 

	  The degeneracy at each super cell (basically, the number is set as 1).
	
	* Line 6 + int(``nrpts``/15) -
	
	  1-3rd columns: The supercell lattice vectors.
	  
	  4-th column: The index of wannier functions at the original cell.

	  5-th column: The index of wannier functions at the supercell.

	  6(7)-th column: The real (imaginary) value.
