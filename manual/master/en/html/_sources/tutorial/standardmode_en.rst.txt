Quick guide to *Standard* mode
==============================

Heisenberg model
----------------

This tutorial should be performed in ::

 samples/CG/Heisenberg/
 
The input file is provided as follows::

 samples/CG/Heisenberg/stan.in
 
In this case, we treat the two-dimensional antiferromagnetic Heisenberg model that has a nearest neighbor spin coupling.

.. math::

   \hat{\mathcal H}=J \sum_{i,j=1}^{4} (\hat{ S }_{i j} \cdot \hat{ S }_{i+1 j} + \hat{ S }_{i j} \cdot \hat{ S }_{i j+1},)

where we use the periodic boundary condition :math:`(S_{15}=S_{51}= S_{11})`.

The input file is as follows::

 model = "Spin"
 method = "CG"
 lattice = "square"
 W = 4
 L = 4
 J = 1.0
 2Sz = 0
 
In this tutorial, J and the number of sites are set to 1 (arbitrary unit) and 16, respectively.

**Log output**
^^^^^^^^^^^^^^

Log messages are outputted to the standard output.
Log files for the calculation procedure are created in the \"output\" directory which is automatically created.
In this example, the following files are outputted. ::

 CHECK_InterAll.dat     Time_CG_EigenVector.dat  zvo_Lanczos_Step.dat  
 CHECK_Memory.dat       WarningOnTransfer.dat    zvo_sz_TimeKeeper.dat
 CHECK_Sdim.dat         zvo_TimeKeeper.dat
 
The details of the outputted files are shown in
:ref:`Subsec:checkchemi` etc.
We execute ::

 $ Path/HPhi -s stan.in
 
and obtain the following standard outputs (the compilation mode is MPI parallel/hybrid parallel) ::


       ,ammmmmmmmmmmmmmb,,       Welcome to the
     ,@@` dm          mb  ===m
   ,@@` d@@@@@@@@@@@@@@@@b Pm,   @@          @@       @@
  d@  d@@@ @@@ @@@@@@ @@@@b ~@a  @@          @@    @@@@@@@@
 d@   @@@@ ^^^ @@@@ m m @@@   @, @@          @@  @@@  @@  @@@
 @    @@@@_@@@_@@@@mm mm@@@   @| @@mmmmmmmmmm@@ @@    @@    @@
 P@    9@@@@@@@@@@@@@@@@@P    @~ @@@@@@@@@@@@@@ @@    @@    @@
  @@      ~~9@@@@@@PPP~      @P  @@          @@  @@@  @@  @@@
   ~@@b      @@@@@@@      ,@@~   @@          @@    @@@@@@@@
     ~@@@m,,@@@@@@@@@  ,m@~`     @@          @@       @@
         ~~9@@@@@@@@@  ~
            9@P~~~9@P            Version 2.0.3


 #####  Parallelization Info.  #####

   OpenMP threads : 1
   MPI PEs : 1


 ######  Standard Intarface Mode STARTS  ######

   Open Standard-Mode Inputfile stan.in

   KEYWORD : model                | VALUE : Spin
   KEYWORD : method               | VALUE : CG
   KEYWORD : lattice              | VALUE : square
   KEYWORD : w                    | VALUE : 4
   KEYWORD : l                    | VALUE : 4
   KEYWORD : j                    | VALUE : 1.0
   KEYWORD : 2sz                  | VALUE : 0

 #######  Parameter Summary  #######

   @ Lattice Size & Shape

                 a = 1.00000     ######  DEFAULT VALUE IS USED  ######
           Wlength = 1.00000     ######  DEFAULT VALUE IS USED  ######
           Llength = 1.00000     ######  DEFAULT VALUE IS USED  ######
                Wx = 1.00000     ######  DEFAULT VALUE IS USED  ######
                Wy = 0.00000     ######  DEFAULT VALUE IS USED  ######
                Lx = 0.00000     ######  DEFAULT VALUE IS USED  ######
                Ly = 1.00000     ######  DEFAULT VALUE IS USED  ######
            phase0 = 0.00000     ######  DEFAULT VALUE IS USED  ######
            phase1 = 0.00000     ######  DEFAULT VALUE IS USED  ######

   @ Super-Lattice setting

                 L = 4
                 W = 4
            Height = 1           ######  DEFAULT VALUE IS USED  ######
          Number of Cell = 16

   @ Hamiltonian

                 h = 0.00000     ######  DEFAULT VALUE IS USED  ######
             Gamma = 0.00000     ######  DEFAULT VALUE IS USED  ######
                2S = 1           ######  DEFAULT VALUE IS USED  ######
                 D = 0.00000     ######  DEFAULT VALUE IS USED  ######
               J0x = 1.00000
               J0y = 1.00000
               J0z = 1.00000
               J1x = 1.00000
               J1y = 1.00000
               J1z = 1.00000

   @ Numerical conditions

        LargeValue = 4.50000     ######  DEFAULT VALUE IS USED  ######

 ######  Print Expert input files  ######

     locspn.def is written.
     coulombinter.def is written.
     hund.def is written.
     exchange.def is written.
     CDataFileHead = zvo         ######  DEFAULT VALUE IS USED  ######
       Lanczos_max = 2000        ######  DEFAULT VALUE IS USED  ######
        initial_iv = -1          ######  DEFAULT VALUE IS USED  ######
              exct = 1           ######  DEFAULT VALUE IS USED  ######
        LanczosEps = 14          ######  DEFAULT VALUE IS USED  ######
     LanczosTarget = 2           ######  DEFAULT VALUE IS USED  ######
            NumAve = 5           ######  DEFAULT VALUE IS USED  ######
     ExpecInterval = 20          ######  DEFAULT VALUE IS USED  ######
            NOmega = 200         ######  DEFAULT VALUE IS USED  ######
          OmegaMax = 72.00000    ######  DEFAULT VALUE IS USED  ######
          OmegaMin = -72.00000   ######  DEFAULT VALUE IS USED  ######
           OmegaIm = 0.04000     ######  DEFAULT VALUE IS USED  ######
               2Sz = 0
      modpara.def is written.

   @ Spectrum

        SpectrumQW = 0.00000     ######  DEFAULT VALUE IS USED  ######
        SpectrumQL = 0.00000     ######  DEFAULT VALUE IS USED  ######
        SpectrumQH = 0.00000     ######  DEFAULT VALUE IS USED  ######
      SpectrumType = szsz        ######  DEFAULT VALUE IS USED  ######
         pair.def is written.


   @ CalcMod

           Restart = none        ######  DEFAULT VALUE IS USED  ######
    InitialVecType = c           ######  DEFAULT VALUE IS USED  ######
        EigenVecIO = none        ######  DEFAULT VALUE IS USED  ######
          CalcSpec = none        ######  DEFAULT VALUE IS USED  ######
      calcmod.def is written.

       ioutputmode = 1           ######  DEFAULT VALUE IS USED  ######
     greenone.def is written.
     greentwo.def is written.
     namelist.def is written.

 ######  Input files are generated.  ######

   Read File 'namelist.def'.
   Read File 'calcmod.def' for CalcMod.
   Read File 'modpara.def' for ModPara.
   Read File 'locspn.def' for LocSpin.
   Read File 'coulombinter.def' for CoulombInter.
   Read File 'hund.def' for Hund.
   Read File 'exchange.def' for Exchange.
   Read File 'greenone.def' for OneBodyG.
   Read File 'greentwo.def' for TwoBodyG.
   Read File 'pair.def' for PairExcitation.

 ######  Definition files are correct.  ######

   Read File 'locspn.def'.
   Read File 'coulombinter.def'.
   Read File 'hund.def'.
   Read File 'exchange.def'.
   Read File 'greenone.def'.
   Read File 'greentwo.def'.
   Read File 'pair.def'.

 ######  Indices and Parameters of Definition files(*.def) are complete.  ######

   MAX DIMENSION idim_max=12870
   APPROXIMATE REQUIRED MEMORY  max_mem=0.001647 GB


 ######  MPI site separation summary  ######

   INTRA process site
     Site    Bit
        0       2
        1       2
        2       2
        3       2
        4       2
        5       2
        6       2
        7       2
        8       2
        9       2
       10       2
       11       2
       12       2
       13       2
       14       2
       15       2

   INTER process site
     Site    Bit

   Process element info
     Process       Dimension   Nup  Ndown  Nelec  Total2Sz   State
           0           12870     8      8      8         0

    Total dimension : 12870


 ######  LARGE ALLOCATE FINISH !  ######

   Start: Calculate HilbertNum for fixed Sz.
   End  : Calculate HilbertNum for fixed Sz.

   Start: Calculate diagaonal components of Hamiltonian.
   End  : Calculate diagaonal components of Hamiltonian.

 ######  Eigenvalue with LOBPCG  #######

   initial_mode=1 (random): iv = -1 i_max=12870 k_exct =1

     Step   Residual-2-norm     Threshold      Energy
         1     2.44343e+00     1.00000e-07          -5.27456e-01
         2     2.76604e+00     1.87217e-07          -1.87217e+00
         3     2.61923e+00     4.19088e-07          -4.19088e+00
         4     2.57106e+00     5.97098e-07          -5.97098e+00

 ( snip )

        40     7.39431e-06     1.12285e-06          -1.12285e+01
        41     4.15948e-06     1.12285e-06          -1.12285e+01
        42     2.04898e-06     1.12285e-06          -1.12285e+01
        43     9.92048e-07     1.12285e-06          -1.12285e+01

 ######  End  : Calculate Lanczos EigenValue.  ######


 ######  End  : Calculate Lanczos EigenVec.  ######

 i=    0 Energy=-11.228483 N= 16.000000 Sz=  0.000000 Doublon=  0.000000
 
In the beginning of this run,
files describing the details of the considered Hamiltonian 
(``locspin.def``, ``trans.def``, ``exchange.def``, ``coulombintra.def``, ``hund.def``, ``namelist.def``, ``calcmod.def``, ``modpara.def``) and files specifying the elements of the correlation functions that will be calculated(``greenone.def``, ``greentwo.def``) are generated.

**Outputs for calculation results**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Locally Optimal Block Conjugate Gradient (LOBCG) method**
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When a calculation by the LOBCG method is finished normally, eigenenergies, one-body Green's functions, and two-body Green's functions are calculated and outputted to the files, respectively. In this sample, the following files are outputted. ::
 
 zvo_energy.dat
 zvo_cisajscktalt_eigen_xx.dat  zvo_phys_Nup4_Ndown4.dat
  
where xx is the number of the eigenstate counting from 0.

**Lanczos method**
""""""""""""""""""

When a calculation by the Lanczos method is completed normally, eigenenergies, one-body Green's functions, and two-body Green's functions are calculated and outputted to the files, respectively. In this sample, the following files are outputted. ::
 
 zvo_energy.dat zvo_cisajs.dat
 zvo_cisajscktalt.dat
 
For Standard mode, all pairs of :math:`\langle n_{i\sigma} \rangle` are calculated as one-body Green's functions and those of :math:`\langle n_{i\sigma} n_{j\sigma'} \rangle` are calculated as two-body Green's functions on the basis of the definition files, ``greenone.def`` and ``greentwo.def``. When the accuracy of the Lanczos vectors is sufficient, one-body and two-body Green's functions are calculated by the eigenvectors obtained by the Lanczos method. When the accuracy of the Lanczos vectors is *not* sufficient, a message \"Accuracy of Lanczos vector is not enough\" is outputted to the standard output and the one-body and two-body Green's functions are calculated by the eigenvectors obtained by CG method. The details of output files are shown in :ref:`Subsec:energy.dat` , :ref:`Subsec:cgcisajs` , :ref:`Subsec:cisajscktalt`.

**TPQ method**
""""""""""""""

When ``method="TPQ"`` is selected in an input file, a calculation by the TPQ method is started. After the calculation is completed normally, the following files are outputted, where \%\% is the number of runs and \&\& is the number of steps for the TPQ method. ::
  
 Norm_rand%%.dat SS_rand%%.dat
 zvo_cisajs_set%%step&&.dat
 zvo_cisajscktalt_set%%step&&.dat
 
In Norm\_rand\%\%.dat, basic information such as the inverse of temperature and the norm of the wave function before normalization is outputted with a TPQ step for each number of runs. In SS\_rand\%\%.dat, physical quantities such as the inverse of temperature, energy, and expected value of the square of the Hamiltonian are outputted with a TPQ step for each number of runs. In zvo\_cisajs\_set\%\%step\&\&.dat and zvo\_cisajscktalt\_set\%\%step\&\&.dat, one-body and two-body Green's functions are outputted for each number of a TPQ steps and runs. The details of these files are shown in :ref:`Subsec:normrand`, :ref:`Subsec:ssrand`, :ref:`Subsec:cgcisajs`, :ref:`Subsec:cisajscktalt`.

**Full diagonalization method**
"""""""""""""""""""""""""""""""

When ``method = "fulldiag"`` is selected in an input file, a calculation by the full diagonalization method is started. After the calculation is completed normally, the following files are outputted, where xx is the number of the eigenstate counting from 0. ::
 
 Eigenvalue.dat zvo_cisajs_eigen_xx.dat
 zvo_cisajscktalt_eigen_xx.dat  zvo_phys_Nup4_Ndown4.dat
   
In Eigenvalue.dat, an eigennumber and an eigenvalue are outputted for each line.In zvo\_cisajs\_eigen\_xx.dat and zvo\_cisajscktalt\_eigen\_xx.dat,one-body Green's functions and two-body Green's functions are outputted for each eigennumber. In zvo\_phys\_Nup4\_Ndown4.dat, physical quantities, such as the expected values of energy and the doublon are outputted. The details of these files are shown in :ref:`Subsec:eigenvalue` - :ref:`Subsec:cisajscktalt`.

Other tutorials
---------------

There are many tutorials in ``samples/Standard/``. For more details, please see ``README.md`` at each directory.
