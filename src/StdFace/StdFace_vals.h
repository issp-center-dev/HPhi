/*
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**@file
@brief Variables used in the Standard mode.
These variables are passed as a pointer of the structure(StdIntList).
*/
#include <complex.h>

struct StdIntList {
  /*
   Initial (undefined)
  */
  int NaN_i;/**<@brief It is used for initializing input parameter. 
            This means that a parameter wich is not specified in input file.
            Set in StdFace_ResetVals().*/
  double pi;/**<@brief @f$\pi=3.14...@f$*/
  /*
  Parameters for LATTICE
  */
  char lattice[256];/**<@brief Name of lattice. Input parameter.*/
  double a; /**<@brief The lattice constant. Input parameter.*/
  double length[3];/**<@brief Anisotropic lattice constant, 
                   input parameter wlength, llength, hlength.*/
  int W;/**<@brief Number of sites along the 1st axis, input parameter.*/
  int L;/**<@brief Number of sites along the 2nd axis, input parameter.*/
  int Height;/**<@brief Number of sites along the 3rd axis, input parameter.*/
  double direct[3][3];/**<@brief The unit direct lattice vector. 
                      Set in StdFace_InitSite().*/
  int box[3][3];/**<@brief The shape of the super-cell. Input parameter
                a0W, a0L, a0H, etc. or defined from StdIntList::W, etc. in 
                StdFace_InitSite().*/
  int rbox[3][3];/**<@brief The inversion of StdIntList::box.
                 Set in StdFace_InitSite().*/
  int NCell;/**<@brief The number of the unit cell in the super-cell
            (determinant of StdIntList::box). Set in StdFace_InitSite().*/
  int **Cell;/**<@brief [StdIntList][3] The cell position in the fractional 
             coordinate. Malloc and Set in StdFace_InitSite().*/
  int NsiteUC;/**<@brief Number of sites in the unit cell. Defined in the
              beginning of each lattice function*/
  double **tau;/**<@brief Cell-internal site position in the fractional 
               coordinate. Defined in the beginning of each lattice function*/
  /*
  Parameters for MODEL
  */
  char model[256];/**<@brief Name of model, input parameter*/
  double mu;/**<@brief Chemical potential, input parameter*/
  double complex t;/**<@brief Nearest-neighbor hopping, input parameter*/
  double complex tp;/**<@brief 2nd-nearest hopping, input parameter*/
  double complex t0;/**<@brief Anisotropic hopping (1st), input parameter*/
  double complex t0p;/**<@brief Anisotropic hopping (2nd), input parameter*/
  double complex t0pp;/**<@brief Anisotropic hopping (3rd), input parameter*/
  double complex t1;/**<@brief Anisotropic hopping (1st), input parameter*/
  double complex t1p;/**<@brief Anisotropic hopping (2nd), input parameter*/
  double complex t1pp;/**<@brief Anisotropic hopping (3rd), input parameter*/
  double complex t2;/**<@brief Anisotropic hopping (1st), input parameter*/
  double complex t2p;/**<@brief Anisotropic hopping (2nd), input parameter*/
  double complex t2pp;/**<@brief Anisotropic hopping (3rd), input parameter*/
  double complex tpp;/**<@brief 3rd-nearest hopping, input parameter*/
  double U;/**<@brief On-site Coulomb potential, input parameter*/
  double V;/**<@brief Off-site Coulomb potential (1st), input parameter*/
  double Vp;/**<@brief Off-site Coulomb potential (2nd), input parameter*/
  double V0;/**<@brief Anisotropic Coulomb potential (1st), input parameter*/
  double V0p;/**<@brief Anisotropic Coulomb potential (2nd), input parameter*/
  double V0pp;/**<@brief Anisotropic Coulomb potential (3rd), input parameter*/
  double V1;/**<@brief Anisotropic Coulomb potential (1st), input parameter*/
  double V1p;/**<@brief Anisotropic Coulomb potential (2nd), input parameter*/
  double V1pp;/**<@brief Anisotropic Coulomb potential (3rd), input parameter*/
  double V2;/**<@brief Anisotropic Coulomb potential (1st), input parameter*/
  double V2p;/**<@brief Anisotropic Coulomb potential (2nd), input parameter*/
  double V2pp;/**<@brief Anisotropic Coulomb potential (3rd), input parameter*/
  double Vpp;/**<@brief Off-site Coulomb potential (3rd), input parameter*/
  /**/
  double JAll;/**<@brief Isotropic, diagonal spin coupling (1st Near.), 
              input parameter J.*/
  double JpAll;/**<@brief Isotropic, diagonal spin coupling (2nd Near), 
               input parameter Jp.*/
  double J0All;/**<@brief Anisotropic, diagonal spin coupling (1st Near), 
               input parameter J0.*/
  double J0pAll;/**<@brief Anisotropic, diagonal spin coupling (2nd Near), 
               input parameter J0'.*/
  double J0ppAll;/**<@brief Anisotropic, diagonal spin coupling (3rd Near),
               input parameter J0''.*/
  double J1All;/**<@brief Anisotropic, diagonal spin coupling (1st Near),
               input parameter J1.*/
  double J1pAll;/**<@brief Anisotropic, diagonal spin coupling (2nd Near), 
               input parameter J1'.*/
  double J1ppAll;/**<@brief Anisotropic, diagonal spin coupling (3rd Near),
               input parameter J1''.*/
  double J2All;/**<@brief Anisotropic, diagonal spin coupling (1st Near),
               input parameter J2.*/
  double J2pAll;/**<@brief Anisotropic, diagonal spin coupling (2nd Near), 
               input parameter J2'.*/
  double J2ppAll;/**<@brief Anisotropic, diagonal spin coupling (3rd Near),
               input parameter J2''.*/
  double JppAll;/**<@brief Isotropic, diagonal spin coupling (3rd Near),
               input parameter J''.*/
  double J[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                 (1st Near.), input parameter Jx, Jy, Jz, Jxy, etc.*/
  double Jp[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J'x, J'y, J'z, J'xy, etc.*/
  double J0[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                 (1st Near.), input parameter J0x, J0y, J0z, J0xy, etc. 
                 or set in StdFace_InputSpinNN().*/
  double J0p[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J0'x, J0'y, J0'z, J0'xy, etc. 
                   or set in StdFace_InputSpin().*/
  double J0pp[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J0''x, J0''y, J0''z, J0''xy, etc.
                   or set in StdFace_InputSpin().*/
  double J1[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                  (1st Near.), input parameter J1x, J1y, J1z, J1xy, etc. 
                  or set in StdFace_InputSpinNN().*/
  double J1p[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J1'x, J1'y, J1'z, J1'xy, etc. 
                   or set in StdFace_InputSpin().*/
  double J1pp[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J1''x, J1''y, J1''z, J1''xy, etc.
                   or set in StdFace_InputSpin().*/
  double J2[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                  (1st Near.), input parameter J2x, J2y, J2z, J2xy, etc. 
                  or set in StdFace_InputSpinNN().*/
  double J2p[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (2nd Near.), input parameter J2'x, J2'y, J2'z, J2'xy, etc. 
                   or set in StdFace_InputSpin().*/
  double J2pp[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J2''x, J2''y, J2''z, J2''xy, etc.
                   or set in StdFace_InputSpin().*/
  double Jpp[3][3];/**<@brief Isotropic, diagonal/off-diagonal spin coupling
                   (3rd Near.), input parameter J''x, J''y, J''z, J''xy, etc.*/
  double D[3][3];/**<@brief Coefficient for @f${\hat S}_{i z} {\hat S}_{i z}@f$
                 input parameter D. Only D[2][2] is used.*/
  double h;/**<@brief Longitudinal magnetic field, input parameter.*/
  double Gamma;/**<@brief Transvars magnetic field, input parameter.*/
  double K;/**<@brief 4-spin term. Not used.*/
  /*
   Phase for the boundary
  */
  double pi180;/**<@brief @f$\pi/180@f$, set in StdFace_ResetVals().*/
  double phase[3];/**<@brief Boundary phase, input parameter phase0, etc.*/
  double complex ExpPhase[3];/**<@brief @f$\exp(i \pi {\rm phase}/180)@f$.*/
  int AntiPeriod[3];/**<@brief If corresponding StdIntList::phase = 180,
                    it becomes 1.*/
  /*
   Transfer, Interaction, Locspin
  */
  int nsite;/**<@brief Number of sites, set in the each lattice file.*/
  int *locspinflag;/**<@brief [StdIntList::nsite] LocSpin in Expert mode, 
                   malloc and set in each lattice file.*/
  int ntrans;/**<@brief Number of transfer, counted in each lattice file.*/
  int **transindx;/**<@brief [StdIntList::ntrans][4] Site/spin indices of 
                  one-body term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_trans().*/
  double complex *trans;/**<@brief [StdIntList::ntrans] Coefficient of 
                  one-body term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_trans().*/
  int nintr;/**<@brief Number of InterAll, counted in each lattice file.*/
  int Lintr;/**<@brief Print interall.def or not, set in PrintInteractions().*/
  int **intrindx;/**<@brief [StdIntList::nintr][8] Site/spin indices of 
                  two-body term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  double complex *intr;/**<@brief [StdIntList::nintr] Coefficient of general
                  two-body term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  int NCintra;/**<@brief Number of intra-site Coulomb interaction,
              counted in each lattice file.*/
  int LCintra;/**<@brief Print coulombintra.def or not, set in PrintInteractions().*/
  int **CintraIndx;/**<@brief [StdIntList::NCintra][1] Site indices of 
                  intra-site Coulomb term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  double *Cintra;/**<@brief [StdIntList::NCintra] Coefficient of intra-site
                  Coulomb term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  int NCinter;/**<@brief Number of inter-site Coulomb interaction,
              counted in each lattice file.*/
  int LCinter;/**<@brief Print coulombinter.def or not, set in PrintInteractions().*/
  int **CinterIndx;/**<@brief [StdIntList::NCinter][2] Site indices of 
                  inter-site Coulomb term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  double *Cinter;/**<@brief [StdIntList::NCinter] Coefficient of inter-site
                  Coulomb term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  int NHund;/**<@brief Number of Hund term, counted in each lattice file.*/
  int LHund;/**<@brief Print hund.def or not, set in PrintInteractions().*/
  int **HundIndx;/**<@brief [StdIntList::NHund][2] Site indices of 
                  Hund term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  double *Hund;/**<@brief [StdIntList::NHund] Coefficient of Hund term, 
               malloc in StdFace_MallocInteractions()
                   and set in StdFace_intr().*/
  int NEx;/**<@brief Number of exchange term, counted in each lattice file.*/
  int LEx;/**<@brief Print exchange.def or not, set in PrintInteractions().*/
  int **ExIndx;/**<@brief [StdIntList::NEx][2] Site indices of 
                  exchange term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  double *Ex;/**<@brief [StdIntList::NEx] Coefficient of exchange term, 
               malloc in StdFace_MallocInteractions()
                   and set in StdFace_intr().*/
  int NPairLift;/**<@brief Number of pair-lift term, counted in each lattice file.*/
  int LPairLift;/**<@brief Print pairlift.def or not, set in PrintInteractions().*/
  int **PLIndx;/**<@brief [StdIntList::NPairLift][2] Site indices of 
                  pair-lift term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_intr().*/
  double *PairLift;/**<@brief [StdIntList::NPairLift] Coefficient of 
                   pair-lift term, malloc in StdFace_MallocInteractions()
                   and set in StdFace_intr().*/
  int NPairHopp;/**<@brief Number of pair-hopping term, counted in each lattice file.*/
  int LPairHopp;/**<@brief Print pairhopp.def or not, set in PrintInteractions().*/
  int **PHIndx;/**<@brief [StdIntList::NPairLift][2] Site indices of
               pair-hopping term, malloc in StdFace_MallocInteractions()
               and set in StdFace_intr().*/
  double *PairHopp;/**<@brief [StdIntList::NPairLift] Coefficient of
                   pair-hopping term, malloc in StdFace_MallocInteractions()
                   and set in StdFace_intr().*/
  int lBoost;
  /*
   Calculation conditions
  */
  int lGC;/**<@brief Switch for computing Grandcanonical ensemble(== 1).
          Setted in StdFace_main() after all keywords are read.*/
  int nelec;/**<@brief Number of electrons, input from file.*/
  int S2;/**<@brief Total spin |S| of a local spin, input from file.*/
  char outputmode[256];/**<@brief Select amount of correlation function,
                       input from file.*/
  char CDataFileHead[256];/**<@brief Header of the output files.
                          Input from file*/
  int Sz2;/**<@brief Total Sz, input from file.*/
  int ioutputmode;/**<@brief Switch associated to StdIntList::outputmode*/
  /*
   Wannier90 mode
  */
  double cutoff_t;/**<@brief Cutoof for the hopping in wannier90, input from file*/
  double cutoff_u;/**<@brief Cutoof for the Coulomb in wannier90, input from file*/
  double cutoff_j;/**<@brief Cutoof for the Hund in wannier90, input from file*/
  double cutoff_length_t; /**<@brief Cutoof for R in wannier90, input from file.*/
  double cutoff_length_U; /**<@brief Cutoof for R in wannier90, input from file.*/
  double cutoff_length_J; /**<@brief Cutoof for R in wannier90, input from file.*/
  int cutoff_tR[3];
  int cutoff_UR[3];
  int cutoff_JR[3];
  double lambda; /**<@brief Tuning parameter of U and J in wannier90, input from file.*/
  double lambda_U; /**<@brief Tuning parameter of U in wannier90, input from file.*/
  double lambda_J; /**<@brief Tuning parameter of J in wannier90, input from file.*/
  char double_counting_mode[256];/**<@brief Select mode of double counting, input from file.*/
  double alpha;/**<@brief Tuning parameter of chemical potential correction in wannier90, input from file.*/

#if defined(_HPhi)
  /*
  HPhi modpara
  */
  char method[256];/**<@brief The name of method, input from file.*/
  char Restart[256];/**<@brief The name of restart mode, input from file.*/
  char InitialVecType[256];/**<@brief The name of initialguess-type, input from file.*/
  char EigenVecIO[256];/**<@brief The name of I/O mode for eigenvector, input from file*/
  char HamIO[256];/**<@brief The name of I/O mode for Hamiltonian, input from file*/
  int FlgTemp;/**<@brief */
  int Lanczos_max;/**<@brief The maxixmum number of iterations, input from file*/
  int initial_iv; /**<@brief the number for generating random number, input from file.*/
  int nvec;/**<@brief */
  int exct;/**<@brief The number of eigenvectors to be computed. input from file*/
  int LanczosEps;/**<@brief Convergence threshold for the Lanczos method.*/
  int LanczosTarget;/**<@brief Which eigenvector is used for the convergence check.*/
  int NumAve;/**<@brief Number of trials for TPQ calculation.*/
  int ExpecInterval;/**<@brief Interval for the iteration when the expectation 
                    value is computed.*/
  double LargeValue;/**<@brief The shift parameter for the TPQ calculation.*/
  /*
  Boost
  */
  int ***list_6spin_pair;/**<@brief */
  int **list_6spin_star;/**<@brief */
  int num_pivot;/**<@brief */
  int ishift_nspin;/**<@brief */
  /*
  Spectrum
  */
  char CalcSpec[256];/**<@brief The name of mode for spectrum, input from file.*/
  char SpectrumType[256];/**<@brief The type of mode for spectrum, input from file.*/
  int Nomega;/**<@brief Number of frequencies, input from file.*/
  double OmegaMax;/**<@brief Maximum of frequency for spectrum, input from file.*/
  double OmegaMin;/**<@brief Minimum of frequency for spectrum, input from file.*/
  double OmegaOrg;/**<@brief Origin of frequency for spectrum, input from file.*/
  double OmegaIm;/**<@brief Imaginary part of frequency.*/
  double SpectrumQ[3];/**<@brief wavenumver (q-vector) in fractional coordinate*/
  int SpectrumBody;/**<@brief one- or two-body excitation, defined from
                   StdIntList::SpectrumType*/
  char OutputExVec[256];/**<@brief The name of output mode for the excited vector, input from file.*/
    /*
  Time evolution
  */
  double dt;/**<@brief Time step*/
  double tshift;/**<@brief Shift of time-step of laser*/
  double tdump;/**<@brief Time scale of dumping*/
  double freq;/**<@brief Frequency of laser*/
  double Uquench;/**<@brief Quenched on-site potential*/
  double VecPot[3];/**<@brief Vector potential*/
  char PumpType[256];/**<@brief The type of pump*/
  int PumpBody;/**<@brief one- or two-body pumping, defined from
                   StdIntList::PumpType*/
  int *npump;/**<@brief [StdIntList::nt] Number of transfer, counted in each lattice file.*/
  int ***pumpindx;/**<@brief [StdIntList::nt][StdIntList::npump][4] Site/spin indices of
                  one-body term, malloc in StdFace_MallocInteractions()
                  and set in StdFace_trans().*/
  double complex **pump;/**<@brief [StdIntList::nt][StdIntList::npump] Coefficient of
                        one-body term, malloc in StdFace_MallocInteractions()
                        and set in StdFace_trans().*/
  double **At;/**<@brief [StdIntList::nt][3] Vector potential.*/
  int ExpandCoef;/**<@brief The number of Hamiltonian-vector operation for the time-evolution*/
#elif defined(_mVMC)
  /*mVMC modpara*/
  char CParaFileHead[256];/**<@brief Header of the optimized wavefunction,
                          input from file*/
  int NVMCCalMode;/**<@brief Optimization(=0) or compute correlation
                  function(=1), input from file.*/
  int NLanczosMode;/**<@brief Power Lanczos(=1), input from file*/
  int NDataIdxStart;/**<@brief Start index of trials, input from file.*/
  int NDataQtySmp;/**<@brief Number of trials, input from file.*/
  int NSPGaussLeg;/**<@brief Number of Gauss-Legendre points for spin projection,
                  input from file.*/
  int NMPTrans;/**<@brief Number of translation symmetry*/
  int NSROptItrStep;/**<@brief Number of iterations for stocastic reconfiguration*/
  int NSROptItrSmp;/**<@brief Number of steps for sampling*/
  int NSROptFixSmp;/**<@brief */
  double DSROptRedCut;/**<@brief Stocastic reconfiguration parameter, input from file.*/
  double DSROptStaDel;/**<@brief Stocastic reconfiguration parameter, input from file.*/
  double DSROptStepDt;/**<@brief Stocastic reconfiguration parameter, input from file.*/
  int NVMCWarmUp;/**<@brief */
  int NVMCInterval;/**<@brief */
  int NVMCSample;/**<@brief */
  int NExUpdatePath;/**<@brief */
  int RndSeed;/**<@brief */
  int NSplitSize;/**<@brief */
  int NSPStot;/**<@brief */
  int NStore;/**<@brief */
  int NSRCG;/**<@brief */
  int ComplexType;/**<@brief */
  /*
   Sub-lattice
  */
  int Lsub;/**<@brief Sublattice*/
  int Wsub;/**<@brief Sublattice*/
  int Hsub;/**<@brief Sublattice*/
  int NCellsub;/**<@brief Number of cells in a sublattice*/
  int boxsub[3][3];/**<@brief Sublattice*/
  int rboxsub[3][3];/**<@brief Sublattice*/
  /*
   2-body part of the trial wavefunction
  */
  int **Orb;/**<@brief [StdIntList::nsite][StdIntList::nsite] Orbital index*/
  int **AntiOrb;/**<@brief [StdIntList::nsite][StdIntList::nsite] Anti-periodic switch*/
  int NOrb;/**<@brief Number of independent orbital index*/
  int NSym;/**<@brief Number of translation symmetries, 
           Defined from the number of cells in the sub-lattice.*/
#endif
};
