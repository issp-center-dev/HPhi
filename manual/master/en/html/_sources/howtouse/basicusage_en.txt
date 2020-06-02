.. highlight:: none

Basic usage
===========

:math:`{\mathcal H}\Phi` has two modes: Standard mode and Expert mode. Here, the basic flows of calculations of the Standard and expert modes are shown.

*Standard* mode
----------------

The procedure of calculation through the standard mode is as follows:

 1. Create a directory for a calculation scenario
 
  First, you create a working directory for the calculation.
  
 2. Create input files for Standard mode
 
  In Standard mode, you can choose a model (the Heisenberg model, Hubbard model, etc.) and a lattice (the square lattice, triangular lattice, etc.) from those provided; you can specify some parameters (such as the first/second nearest neighbor hopping integrals and the on-site Coulomb integral) for them. Finally, you have to specify the numerical method (such as the Lanczos method) employed in this calculation. The input file format is described in :ref:`How to use HPhi <Ch:Prerequisite>`.
  
 3. Run
 
  Run an executable ``HPhi`` in the terminal by setting option \"``-s``\" (or \"``--standard``\") and the name of the input file written in the previous step.
  
   * Serial/OpenMP parallel
   
   ``$ Path /HPhi -s Input_file_name``
   
   * MPI parallel/ Hybrid parallel
   
   ``$ mpiexec -np number_of_processes Path /HPhi -s Input_file_name``
   
  When you use a queuing system in workstations or super computers, sometimes the number of processes is specified as an argument for the job-submitting command. If you need more information, please refer to your system manuals. The number of processes depends on the target system of the models. The details of setting the number of processes are shown in :ref:`Subsec:CreatingExpert` .
  
 4. Watch calculation logs
 
  Log files are outputted in the \"output\" folder, which is automatically created in the directory for a calculation scenario. The details of the output files are shown in :ref:`Sec:outputfile` .
  
 5. Results 
 
  If the calculation is completed normally, the result files are outputted in  the \"output\" folder. The details of the output files are shown in :ref:`Sec:outputfile` 
  
.. tip::

 | **The number of threads for OpenMP**
 | If you specify the number of OpenMP threads for :math:`{\mathcal H}\Phi`, you should set it as follows (in the case of 16 threads) before running:
 | ``export OMP_NUM_THREADS=16``
  
*Expert* mode
-------------

The calculation procedure for Expert mode is as follows. 
 1. Create a directory for a calculation scenario
 
  First, you create a directory and give it the name of a calculation scenario (you can attach an arbitrary name to a directory).
  
 2. Create input files for Expert mode
 
  For Expert mode, you should create input files for constructing Hamiltonian operators, calculation conditions, and a list file for the filenames of the input files (see the file formats shown in :ref:`Ch:HowToExpert`).
  
 .. note::

  | A list file can be easily created by using Standard mode.
  
 3. Run
 
  Run :math:`{\mathcal H}\Phi` in the terminal by setting option \"``-e``\" (or \"``--expert``\") and the file name of a list file.
  
   * Serial/OpenMP
   
    ``$ Path/HPhi -e Input_List_file\_name``
   
   * MPI/Hybrid
   
    | ``$ mpiexec -np number_of_processes Path/HPhi -e Input_List_file_name``
    | A number of processes depend on a target of system for models. The details of setting a number of processes are shown in :ref:`Subsec:CreatingExpert`.
   
 4. While running
 
  Log files are outputted in the \"output\" folder which is automatically created in the directory for a calculation scenario. The details of the output files are shown in :ref:`Sec:outputfile`.
  
 5. Results
 
  If the calculation is finished normally, the result files are outputted in the \"output\" folder. The details of the output files are shown in :ref:`Sec:outputfile`. 
  
.. _Subsec:CreatingExpert:
  
Creating input files for *Expert* mode
--------------------------------------
  
This mode is for creating input files for *Expert* mode.
A set of input files created using this mode gives a model provided in *Standard* mode.
The usage is shown as follows.

 1. Create an input file for *Standard* mode.
 
 2. Setting an option \"-sdry\" and an input file (in this example, StdFace.def), run :math:`{\mathcal H}\Phi`.
    ::

     $ Path/HPhi -sdry StdFace.def
     
    In this case, you should not use MPI parallelization (mpirun, mpiexec, etc.).

 3. The following files are created as the input files for *Expert* mode in the current working directory.

    ::

     calcmod.def   greentwo.def  namelist.def  zTrans.def
     greenone.def  modpara.def   zInterAll.def zlocspn.def
  
.. _Subsec:process:
  
Setting the process number for MPI/hybrid parallelization
---------------------------------------------------------

For using MPI/hybrid parallelization, the process number must be set as follows.

 1. Standard mode
 
  * Hubbard/Kondo model
  
   When ``model`` in the input file for Standard mode is set as ``"Fermion Hubbard"``, ``"Kondo Lattice"``, or ``"Fermion HubbardGC"``, the process number must be equal to :math:`4^n`.
   
  * Spin model
  
   When ``model`` in the input file for Standard mode is set as ``"Spin"`` or ``"SpinGC"``, the process number must be equal to :math:`(2S+1)^n`, where ``2S`` is set in the input file (the default value is :math:`1`).
   
 2. Expert mode
 
  * Hubbard/Kondo model
  
   When the model is selected as the Fermion Hubbard model or Kondo model by setting ``CalcModel`` in a **CalcMod** file, the process number must be equal to :math:`4^n`. See :ref:`Subsec:calcmod` for details of the ``CalcModel`` file. 
   
  * Spin model
  
   When the model is selected as the spin model by setting ``CalcModel`` in a **CalcMod** file, the process number is fixed by a **LocSpin** file. The process number must be equal to the number calculated by multiplying the state number of the localized spin (``2S`` +1) in descending order by the site number. See :ref:`Subsec:locspn` for details of the **LocSpin** file.
   
   For example, when a **LocSpin** file is given as follows, the process number must be equal to :math:`2=1+1,~6=2\times(2+1),~24=6\times(3+1)`. 

  ::
  
   ================================ 
   NlocalSpin     3
   ================================
   ========i_0IteElc_2S ======
   ================================
       0      3
       1      2
       2      1

Printing version ID
-------------------

By using the ``-v`` option as follows, you can check which version of :math:`{\mathcal H}\Phi` you are using.

 ``$ PATH/HPhi -v``
