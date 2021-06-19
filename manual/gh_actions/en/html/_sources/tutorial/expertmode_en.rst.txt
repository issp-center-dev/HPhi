.. highlight:: none

Quick guide to *Expert* mode
============================

For Expert mode, the following input files are needed.

1. A file list for input files
2. Files for basic parameters
3. Files for constructing Hamiltonian
4. Files for setting output components.

The process after calculation is the same as in Standard mode.
In this section, we demonstrate Expert mode in the directory where
the tutorial at the previous section was performed.

File list for input files
-------------------------

In namelist.def, the types of input files and filenames are defined as shown below. By writing the keyword and filenames at each line, the types of files are distinguished. The details of namelist.def are shown in :ref:`Subsec:InputFileList`. ::

        ModPara  modpara.def
        LocSpin  locspn.def
   CoulombInter  coulombinter.def
           Hund  hund.def
       Exchange  exchange.def
       OneBodyG  greenone.def
       TwoBodyG  greentwo.def
        CalcMod  calcmod.def
 PairExcitation  pair.def
    SpectrumVec  zvo_eigenvec_0
    
Files for basic parameters
--------------------------

In this subsection, we show how to set a calculation mode, the parameters for the calculation, and the positions of the localized spins.

**Setting a calculation mode**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The calculation mode is set in a CalcMod file (in this sample file, calcmod.def). The contents of the files are as follows. ::

 #CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, ...
 #CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, ...
 #Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart
 #CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, ...
 CalcType   3
 CalcModel   1
 ReStart   0
 CalcSpec   0
 CalcEigenVec   0
 InitialVecType   0
 InputEigenVec   0
  
We select a calculation method in CalcType and a target model in CalcModel. In this sample, we set the Lanczos method as a calculation method and the target model as the spin system (canonical ensemble). The details of a CalcMod file are shown in :ref:`Subsec:calcmod`.

**Setting parameters for calculation**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters for the calculation are set in a ModPara file(in this sample, modpara.def). The contents of this file are as follows. ::

 --------------------
 Model_Parameters   0
 --------------------
 HPhi_Cal_Parameters
 --------------------
 CDataFileHead  zvo
 CParaFileHead  zqp
 --------------------
 Nsite          16   
 2Sz            0    
 Lanczos_max    2000 
 initial_iv     -1   
 exct           1    
 LanczosEps     14   
 LanczosTarget  2    
 LargeValue     4.500000000000000e+00    
 NumAve         5    
 ExpecInterval  20   
 NOmega         200  
 OmegaMax       7.200000000000000e+01     4.000000000000000e-02    
 OmegaMin       -7.200000000000000e+01    4.000000000000000e-02    
 OmegaOrg       0.0 0.0
  
In this file, we set the parameters for the calculation, such as the site number, the total number of conduction electrons, the total :math:`S_z` and the number of Lanczos steps. The details of the ModPara file are shown in :ref:`Subsec:modpara`.
  
**Setting positions of localized spins**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The positions :math:`S` of the localized spins are defined by a LocSpin file (in this sample, locspn.def). The contents of the files are as follows. ::

 ================================
 NlocalSpin    16  
 ================================ 
 ========i_0IteElc_1LocSpn ====== 
 ================================ 
     0      1
     1      1
     2      1
     3      1
     4      1
     5      1
 ...
 
When CalcModel in a CalcMod file is set as the spin system, all the sites are automatically treated as localized spins. The details of a LocSpin file are shown in :ref:`Subsec:locspn`.

Files for constructing Hamiltonian
----------------------------------

After setting the basic parameters, we create input files for constructing the Hamiltonian. Since the calculations are performed by using the representation of the fermion operators in :math:`{\mathcal H}\Phi`, we must rewrite the spin operator. For example,  in the case of :math:`S = 1/2`, we rewrite the equation by using the relation

.. math::

    S^z_{i}&=(c_{i\uparrow}^{\dagger}c_{i\uparrow}-c_{i\downarrow}^{\dagger}c_{i\downarrow})/2,\\
    S^+_{i}&=c_{i\uparrow}^{\dagger}c_{i\downarrow},\\
    S^-_{i}&=c_{i\downarrow}^{\dagger}c_{i\uparrow}.

**Setting transfer integrals**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a Trans file (in this sample, zTrans.def), we set the transfer part of the Hamiltonian,

.. math::

   \mathcal{H} +=-\sum_{ij\sigma_1\sigma2}
   t_{ij\sigma_1\sigma2}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}.
   
The contents of the files are as follows. ::

 ======================== 
 NTransfer       0  
 ======================== 
 ========i_j_s_tijs====== 
 ======================== 
   
We can use this term when an electric magnetic field is added in the spin system. For example, when a magnetic field is added at a site 1 such as :math:`-0.5 S_z^{(1)}` for :math:`S=1/2`, this term can be rewritten as :math:`-0.5/2(c_{1\uparrow}^{\dagger}c_{1\uparrow}-c_{1\downarrow}^{\dagger}c_{1\downarrow})`. Thus, the input file becomes as follows. ::

 ======================== 
 NTransfer      1   
 ======================== 
 ========i_j_s_tijs====== 
 ======================== 
 1 0 1 0 -0.25 0
 1 1 1 1 0.25 0
 
The details for a Trans file are shown in :ref:`Subsec:Trans`.

**Setting general two-body interactions**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In an InterAll file (in this sample, zInterall.def), we set the general two-body interaction part of the Hamiltonian,

.. math::

   \mathcal{H}+=\sum_{i,j,k,l}\sum_{\sigma_1,\sigma_2, \sigma_3, \sigma_4}I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}.

The contents of the files are as follows. ::

 ====================== 
 NInterAll      96  
 ====================== 
 ========zInterAll===== 
 ====================== 
     0     0     0     0     1     0     1     0   0.500000  0.000000
     0     0     0     0     1     1     1     1  -0.500000  0.000000
     0     1     0     1     1     0     1     0  -0.500000  0.000000
     0     1     0     1     1     1     1     1   0.500000  0.000000
     0     0     0     1     1     1     1     0   1.000000  0.000000
     0     1     0     0     1     0     1     1   1.000000  0.000000
 ...
 
Here, we explain the interaction between site :math:`i` and site :math:`j` in the case of :math:`S = 1/2`, for simplicity. Using fermion operators, the interaction terms for the spin operators can be rewritten as

.. math::
   \mathcal{H}_{i,i+1}&=J(S^x_{i}S^x_{i+1}+S^y_{i}S^y_{i+1}+S^z_{i}S^z_{i+1}) \nonumber\\
   &=J \left( \frac{1}{2}S^+_{i}S^-_{i+1}+\frac{1}{2}S^-_{i}S^+_{i+1}+S^z_{i}S^z_{i+1} \right) \nonumber\\
   &=J \left[ \frac{1}{2}c_{i\uparrow}^{\dagger}c_{i\downarrow}c_{i+1\downarrow}^{\dagger}c_{i+1\uparrow}+\frac{1}{2}c_{i\downarrow}^{\dagger}c_{i\uparrow}c_{i+1\uparrow}^{\dagger}c_{i+1\downarrow}+\frac{1}{4}(c_{i\uparrow}^{\dagger}c_{i\uparrow}-c_{i\downarrow}^{\dagger}c_{i\downarrow})(c_{i+1\uparrow}^{\dagger}c_{i+1\uparrow}-c_{i+1\downarrow}^{\dagger}c_{i+1\downarrow}) \right]. \nonumber 

Thus, the interaction :math:`S^z_{i}S^z_{i+1}` for :math:`J=2` can be written as ::

    i     0     i     0    i+1     0    i+1     0   0.500000  0.000000
    i     0     i     0    i+1     1    i+1     1  -0.500000  0.000000
    i     1     i     1    i+1     0    i+1     0  -0.500000  0.000000
    i     1     i     1    i+1     1    i+1     1   0.500000  0.000000
  
in the format of an InterAll file. The other terms can be written as follows. ::

    i     0     i     1    i+1     1    i+1     0   1.000000  0.000000
    i     1     i     0    i+1     0    i+1     1   1.000000  0.000000
  
There are other file formats for constructing the Hamiltonian. The details of the input formats of two-body interactions are shown in :ref:`Subsec:interall` - :ref:`Subsec:pairlift`.

Setting output components
-------------------------

In OneBodyG and TwoBodyG files, the indices of one-body and two-body Green's functions are defined, respectively. 

**Setting indices of one-body Green's functions**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a OneBodyG file (in this sample, greenone.def), the indices of :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2} \rangle` are defined. The contents of files are as follows. ::

 ===============================
 NCisAjs         32
 ===============================
 ======== Green functions ======
 ===============================
    0     0     0     0
    0     1     0     1
    1     0     1     0
    1     1     1     1
    2     0     2     0
 ...
 
The details of the input formats of a OneBodyG file are shown in :ref:`Subsec:onebodyg`.

**Setting indices of two-body Green's functions**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the TwoBodyG file (in this sample, greentwo.def), the indices of :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4} \rangle` are defined. The contents of this file are as follows. ::

 =============================================
 NCisAjsCktAltDC       1024
 =============================================
 ======== Green functions for Sq AND Nq ======
 =============================================
    0     0     0     0     0     0     0     0
    0     0     0     0     0     1     0     1
    0     0     0     0     1     0     1     0
    0     0     0     0     1     1     1     1
    0     0     0     0     2     0     2     0
 ...

The details of the input formats of the TwoBodyG file are shown in :ref:`Subsec:twobodyg`.

Running
-------

After creating all the input files above, we are ready to run a program. For Expert mode, we must set an option \"-e\" and a file name list (in this sample, namelist.def) as arguments to run :math:`{\mathcal H}\Phi`. ::

 $ Path/HPhi -e namelist.def

The process after the calculation is the same as that of Standard mode.
