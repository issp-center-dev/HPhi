/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*-------------------------------------------------------------
 *[ver.2009.3.31]
 * Exact Diagonalization with Lanczos Method
 *-------------------------------------------------------------
 * Copyright (C) 2006- Takahiro MISAWA. All rights reserved.
 *-------------------------------------------------------------*/
/**@file
@brief Binded struct
*/
#ifndef HPHI_STRUCT_H
#define HPHI_STRUCT_H

#include "Common.h"

/*=================================================================================================*/
//For TEM
struct ParamList {
    //For Time Evolution
    double Tinit;
    double TimeSlice;
    int OutputInterval;
    int ExpandCoef;
    int ExpecInterval;
};
/**
@brief Definision of system (Hamiltonian) etc.
*/
struct DefineList {
  char *CDataFileHead;/**<@brief Read from Calcmod in readdef.h. 
                      Header of output file such as Green's function*/
  char *CParaFileHead;/**<@brief Read from Calcmod in readdef.h.
                      It is not used. Just for the compatibility to mVMC*/
  unsigned int nvec;/**<@brief Read from Calcmod in readdef.h*/
  unsigned int k_exct;/**<@brief Read from Calcmod in readdef.h*/
  int LanczosEps;/**<@brief log(10 base) of the convergence threshold.
                 Read from Calcmod in readdef.h*/
  int LanczosTarget;/**<@brief Which eigenstate is used to check convergence.
                    Read from Calcmod in readdef.h.*/
  int read_hacker;/**<@brief Whether use an efficient method (=1) in sz.c or not (=0)*/
  int READ;/**<@brief It is ALWAYS 0 ???*/
  int WRITE;/**<@brief It is ALWAYS 0 ???*/

  unsigned int Nsite;/**<@brief Number of sites in the INTRA process region*/
  unsigned int NsiteMPI;/**<@brief Total number of sites, differ from DefineList::Nsite*/
  unsigned int Nup;/**<@brief Number of spin-up electrons in this process. */
  unsigned int Ndown;/**<@brief Number of spin-down electrons in this process. */
  unsigned int NupMPI;/**<@brief Total number of spin-up electrons across processes.
                      Deffer from DefineList::Nup. Read from modpara in readdef.h*/
  unsigned int NdownMPI;/**<@brief Total number of spin-down electrons across processes.
                      Deffer from DefineList::Ndown. Read from modpara in readdef.h*/
  unsigned int NupOrg;/**<@brief Number of spin-up electrons before exitation. Used only in
                      the spectrum calculation. Read from modpara in readdef.h*/
  unsigned int NdownOrg;/**<@brief Number of spin-down electrons before exitation. Used only in
                      the spectrum calculation. Read from modpara in readdef.h*/

  int Total2Sz;/**<@brief Total @f$2S_z@f$ in this process.*/
  int Total2SzMPI;/**<@brief Total @f$2S_z@f$ across processes.*/
  unsigned int Ne;/**<@brief Number of electrons in this process.*/
  unsigned int NeMPI;/**<@brief Total number of electrons across process.
                     Differ from DefineList::Ne .*/
  unsigned int Lanczos_max;/**<@brief Maximum number of iterations.*/
  int Lanczos_restart;/**<@brief Number of iterations performed in the restart computation.*/
  long int initial_iv;/**<@brief Seed of random number for initial guesss of wavefunctions.*/

  int istep;/**<@brief Index of TPQ step ???*/
  int irand;/**<@brief Input keyword TargetTPQRand ???*/
  int St;/**<@brief 0 or 1, but it affects nothing.*/

  int *LocSpn;/**<@brief [DefineList::NLocSpn] Flag (and size) of the local spin.
              malloc in setmem_def().*/
  unsigned int NLocSpn;/**<@brief Number of local spins*/
  unsigned int NCond;/**<@brief Number of itinerant electrons*/
  int iFlgGeneralSpin;/**<@brief Flag for the general (Sz/=1/2) spin*/
  int iFlgSzConserved;/**<@brief Flag whether Sz is conserved.*/

  int fidx;/**<@brief Always 0, it is not used ???*/
  long unsigned int *Tpow;/**<@brief [2 * DefineList::NsiteMPI] @f$2^n@f$
                          malloc in setmem_def().*/
  long unsigned int *OrgTpow;/**<@brief [2 * DefineList::NsiteMPI] @f$2^n@f$
                             malloc in setmem_def().*/
  long int *SiteToBit;/**<@brief [DefineList::NsiteMPI] Similar to DefineList::Tpow.
                      For general spin.*/

  unsigned int EDNChemi;/**<@brief Number of on-site term.*/
  int *EDChemi;/**<@brief [DefineList::Nsite] Chemical potential. malloc in setmem_def().*/
  int *EDSpinChemi;/**<@brief [DefineList::Nsite]*/
  double *EDParaChemi;/**<@brief [DefineList::Nsite] On-site potential parameter.
                      malloc in setmem_def().*/

    //[s] Transfer
  unsigned int NTransfer;/**<@brief Number of transfer integrals obtained by a def file.*/
  unsigned int EDNTransfer;/**<@brief Number of transfer integrals for calculation. */
  int **GeneralTransfer;/**<@brief Index of transfer integrals obtained by a def file. 
                        malloc in setmem_def().\n
                        Data Format [DefineList::NTransfer][4]: 
                        0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  int **EDGeneralTransfer;/**<@brief Index of transfer integrals for calculation. 
                          malloc in setmem_def().\n
                          Data Format [DefineList::NTransfer][4]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  double complex *ParaGeneralTransfer;/**<@brief Value of general transfer integrals by a def file. 
                                      malloc in setmem_def().\n
                                      Data Format [DefineList::NTransfer].*/
  double complex *EDParaGeneralTransfer;/**<@brief Value of general transfer integrals  by a def file. 
                                        malloc in setmem_def().\n
                                        Data Format [DefineList::NTransfer].*/
    //[e] Transfer

  unsigned int NCoulombIntra;/**< Number of on-site Coulomb interaction*/
  int **CoulombIntra;/**< [DefineList::NCoulombIntra][1] Index of on-site coulomb interaction.
                     malloc in setmem_def().*/
  double *ParaCoulombIntra;/**< [DefineList::NCoulombIntra] Coupling constant of on-site 
                           Coulomb interaction. malloc in setmem_def().*/

  unsigned int NCoulombInter;/**<@brief Number of off-site Coulomb interaction*/
  int **CoulombInter;/**< [DefineList::NCoulombInter][2] Index of off-site coulomb interaction.
                     malloc in setmem_def().*/
  double *ParaCoulombInter;/**<@brief [DefineList::NCoulombInter]Coupling constant of off-site 
                           Coulomb interaction. malloc in setmem_def().*/

  unsigned int NHundCoupling;/**<@brief Number of Hund coupling*/
  int **HundCoupling;/**<@brief [DefineList::NHundCoupling][2] Index of Hund coupling.
                     malloc in setmem_def().*/
  double *ParaHundCoupling;/**<@brief [DefineList::NHundCoupling] Hund coupling constant.
                           malloc in setmem_def().*/

  unsigned int NPairHopping;/**<@brief Number of pair-hopping term*/
  int **PairHopping;/**<@brief [DefineList::NPairHopping][2] Index of pair-hopping.
                    malloc in setmem_def().*/
  double *ParaPairHopping;/**<@brief [DefineList::NPairHopping] Coupling constant of
                          pair-hopping term. malloc in setmem_def().*/

  unsigned int NExchangeCoupling;/**<@brief Number of exchange term*/
  int **ExchangeCoupling;/**<@brief [DefineList::NExchangeCoupling][2] Index of exchange term.
                         malloc in setmem_def().*/
  double *ParaExchangeCoupling;/**<@brief [DefineList::NExchangeCoupling] Coupling constant of
                               exchange term. malloc in setmem_def().*/

  unsigned int NIsingCoupling;/**<@brief Number of Ising term.*/

  unsigned int NPairLiftCoupling;/**<@brief Number of pair-lift term*/
  int **PairLiftCoupling;/**<@brief [DefineList::NPairHopping][2] Index of pair-lift term.
                         malloc in setmem_def().*/
  double *ParaPairLiftCoupling;/**<@brief [DefineList::NPairHopping] Coupling constant of
                               pair-lift term. malloc in setmem_def().*/

    //[s] For InterAll
  int **InterAll;/**<@brief [DefineList::NinterAll][8] Interacted quartet*/
  int **InterAll_OffDiagonal;/**<@brief [DefineList::NinterAll_OffDiagonal][8] Interacted quartet*/
  int **InterAll_Diagonal;/**<@brief [DefineList::NinterAll_Diagonal][4] Interacted quartet*/
  unsigned int NInterAll;/**<@brief Total Number of Interacted quartet*/
  unsigned int NInterAll_Diagonal;/**<@brief Number of interall term (diagonal)*/
  unsigned int NInterAll_OffDiagonal;/**<@brief Number of interall term (off-diagonal)*/
  double complex *ParaInterAll;/**<@brief [DefineList::NInterAll] Coupling constant of
                               inter-all term. malloc in setmem_def().*/
  double *ParaInterAll_Diagonal;/**<@brief [DefineList::NInterAll_Diagonal] Coupling constant of
                               diagonal inter-all term. malloc in setmem_def().*/
  double complex *ParaInterAll_OffDiagonal;/**<@brief [DefineList::NInterAll_OffDiagonal] Coupling constant of
                               off-diagonal inter-all term. malloc in setmem_def().*/
    //[e] For InterAll

  int **CisAjt;/**<@brief [DefineList::NCisAjt][4] Indices of one-body correlation function. malloc in setmem_def().*/
  unsigned int NCisAjt;/**<@brief Number of indices of two-body correlation function.*/

  int **CisAjtCkuAlvDC;/**<@brief [DefineList::NCisAjtCkuAlvDC][8] Indices of two-body correlation function. malloc in setmem_def().*/
  unsigned int NCisAjtCkuAlvDC;/**<@brief Number of indices of two-body correlation function.*/

  int **TBody;/**<@brief [DefineList::TBody][12] Indices of three-body correlation function. malloc in setmem_def().*/
  unsigned int NTBody;/**<@brief Number of indices of three-body correlation function.*/

  int **FBody;/**<@brief [DefineList::FBody][16] Indices of four-body correlation function. malloc in setmem_def().*/
  unsigned int NFBody;/**<@brief Number of indices of four-body correlation function.*/

  int **SBody;/**<@brief [DefineList::SBody][24] Indices of six-body correlation function. malloc in setmem_def().*/
  unsigned int NSBody;/**<@brief Number of indices of six-body correlation function.*/

  int **SingleExcitationOperator;/**<@brief [DefineList::NSingleExcitationOperator][3] 
                                 Indices of single excitaion operator for spectrum. malloc in setmem_def().*/
  unsigned int NSingleExcitationOperator;/**<@brief Number of single excitaion operator for spectrum.*/
  double complex *ParaSingleExcitationOperator;/**<@brief [DefineList::NSingleExcitationOperator] 
              Coefficient of single excitaion operator for spectrum. malloc in setmem_def().*/

  int **PairExcitationOperator;/**<@brief [DefineList::NPairExcitationOperator][5] 
                               Indices of pair excitaion operator for spectrum. malloc in setmem_def().*/
  unsigned int NPairExcitationOperator;/**<@brief Number of pair excitaion operator for spectrum.*/
  double complex *ParaPairExcitationOperator;/**<@brief [DefineList::NPairExcitationOperator]
                           Coefficient of pair excitaion operator for spectrum. malloc in setmem_def().*/
  
  int iCalcType;/**<@brief Switch for calculation type. 0:Lanczos, 1:TPQCalc, 2:FullDiag.*/
  int iCalcEigenVec;/**<@brief Switch for method to calculate eigenvectors. 
                    0:Lanczos+CG, 1: Lanczos. default value is set as 0 in readdef.c*/
  int iInitialVecType;/**<@brief Switch for type of inital vectors. 
                      0:complex type, 1: real type. default value is set as 0 in readdef.c*/
  int iFlgFiniteTemperature;/**<@brief ???*/
  int iCalcModel;/**<@brief Switch for model. 0:Hubbard, 1:Spin, 2:Kondo, 
                 3:HubbardGC, 4:SpinGC, 5:KondoGC, 6:HubbardNConserved*/
  int iOutputMode;/**<@brief Switch for output mode. 0: OneBodyG and TwoBodyG. 
                  1: OneBodyG and TwoBodyG and correlations for charge and spin.*/
  int iOutputEigenVec;/**<@brief ASwitch for outputting an eigenvector. 0: no output, 1:output.*/
  int iInputEigenVec;/**<@brief Switch for reading an eigenvector. 0: no input, 1:input*/
  int iOutputHam;/**<brief Switch for outputting a Hamiltonian. 0: no output, 1:output*/
  int iInputHam;/**<brief Switch for reading a Hamiltonian. 0: no input, 1:input*/
  int iOutputExVec; /**<brief Switch for outputting an excited vector. 0: no output, 1:output*/

    //[s] For Spectrum
  double complex dcOmegaMax;/**<@brief Upper limit of the frequency for the spectrum.*/
  double complex dcOmegaMin;/**<@brief Lower limit of the frequency for the spectrum.*/
  double complex dcOmegaOrg;/**<@brief Origin limit of the frequency for the spectrum.*/
  int iNOmega;/**<@brief Number of frequencies for spectrum.*/
  int iFlgSpecOmegaMax;/**<@brief Whether DefineList::dcOmegaMax is input or not.*/
  int iFlgSpecOmegaMin;/**<@brief Whether DefineList::dcOmegaMin is input or not.*/
  int iFlgSpecOmegaOrg;/**<@brief Whether DefineList::dcOmegaOrg is input or not.*/
  int iFlgCalcSpec;/**<@brief Input parameter CalcSpec in teh CalcMod file.*/
  int iFlagListModified;/**<@brief When the Hilbert space of excited state differs from the original one.*/

    //[e] For Spectrum

  int iReStart;/**< An integer for restarting output a Hamiltonian.
     - 0: not restart
     - 1:restart (output restart vector),
     - 2: restart (input and output restart vector) */
  int iFlgMPI;/**<@brief MPI mode
    - 0: butterfly
    - 1: Parallel Interaction [to be supported]
    */

    int iNGPU;/**<@brief GPU mode ( only for FullDiag )
    - 0: Use lapack
    - >0: Use GPU
    */

    int iFlgScaLAPACK;/**<@brief ScaLAPACK mode ( only for FullDiag )
    - 0: Use lapack
    - 1: Use ScaLAPACK
    */


    struct ParamList Param;

    //[s] For Time Evolution

    //Information of Time
    unsigned int NTETimeSteps;
    double *TETime;

    //[s]For Ido-san version
    unsigned int NLaser;
    double *ParaLaser;
    //[e]For Ido-san version

    //Information of Transfer integrals
    unsigned int NTETransferMax;
    unsigned int *NTETransfer;        /**< Number of time-dependent transfer integrals for Time Evolution.\n
               Data Format [NTE]*/
    unsigned int *NTETransferDiagonal;        /**< Number of time-dependent transfer integrals for Time Evolution.\n
               Data Format [NTE]*/
    int ***TETransfer;      /**< Index of time-dependent transfer integrals for Time Evolution. \n
               Data Format [NTE][Ntransfer][4]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
    int ***TETransferDiagonal;      /**< Index of time-dependent transfer integrals for Time Evolution. \n
               Data Format [NTE][Ntransfer][2]: 0->site number i, 1-> spin index on i. */
    double complex **ParaTETransfer;  /**< Value of time-dependent transfer integrals for Time Evolution. \n
               Data Format [NTE][Ntransfer]. */
    double **ParaTETransferDiagonal;  /**< Value of time-dependent transfer integrals for Time Evolution. \n
               Data Format [NTE][Ntransfer]. */

    //Two-body part
    unsigned int NTEInterAllMax;
    unsigned int *NTEInterAll;        /**< Number of time-dependent InterAll for Time Evolution.\n
               Data Format [NTE]*/
    unsigned int *NTEInterAllOffDiagonal;        /**< Number of off-diagonal part of time-dependent InterAll for Time Evolution.\n
               Data Format [NTE]*/

    unsigned int *NTEInterAllDiagonal;        /**< Number of diagonal part of time-dependent InterAll for Time Evolution.\n
               Data Format [NTE]*/
    int ***TEInterAll;      /**< Index of time-dependent InterAll for Time Evolution. \n
               Data Format [NTE][NTEInterAll][8]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j.
               4->site number k, 5-> spin index on k, 6-> site number l, 7-> spin index on l.*/
    int ***TEInterAllOffDiagonal;      /**< Index of off-diagonal part of time-dependent InterAll for Time Evolution. \n
               Data Format [NTE][NTEInterAll][8]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j.
               4->site number k, 5-> spin index on k, 6-> site number l, 7-> spin index on l.*/
    int ***TEInterAllDiagonal;      /**< Index of diagonal part of time-dependent InterAll for Time Evolution. \n
               Data Format [NTE][NTEInterAll][4]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
    double complex **ParaTEInterAll;  /**< Value of time-dependent InterAll for Time Evolution. \n
               Data Format [NTE][Ntransfer]. */
    double complex **ParaTEInterAllOffDiagonal;  /**< Value of off-diagonal part of time-dependent InterAll for Time Evolution. \n
               Data Format [NTE][Ntransfer]. */

    double **ParaTEInterAllDiagonal;  /**< Value of diagonal part of time-dependent InterAll for Time Evolution. \n
               Data Format [NTE][Ntransfer]. */
    int **TEChemi;    /**< [NTE][Nsite] */
    unsigned int *NTEChemi;   /**< [NTE] */
    int **SpinTEChemi;  /**< [NTE][Nsite] */
    double **ParaTEChemi;  /**< [NTE][Nsite] */
    //[e] For Time Evolution
    int PreCG;/**< 1: Use preconditioning in LOBPCG */
};/*struct DefineList*/
/**
@brief Size of the Hilbert space
*/
struct CheckList {
  unsigned long int idim_max;/**<@brief The dimension of the Hilbert space of this process.*/
  unsigned long int idim_maxMPI;/**<@brief The total dimension across process.*/
  unsigned long int idim_maxOrg;/**<@brief The local Hilbert-space dimention of original state for the spectrum.*/
  unsigned long int idim_maxMPIOrg;/**<@brief The global Hilbert-space dimention of original state for the spectrum.*/
  unsigned long int sdim;/**<@brief Dimension for Ogata-Lin ???*/
  double max_mem;/**<@brief Estimated memory size.*/
};/*struct CheckList*/
/**
@brief For Matrix-Vector product
*/
struct LargeList {
  double complex prdct;/**<@brief The expectation value of the energy.*/
  int itr;/**<@brief Iteration number.*/
  long int iv;/**<@brief Used for initializing vector.*/
  long int i_max;/**<@brief Length of eigenvector*/
  long int SizeOflist_2_1;/**<@brief Size of ::list_2_1*/
  long int SizeOflist_2_2;/**<@brief Size of ::list_2_2*/
  long int SizeOflistjb;/**<@brief Used for computing Sz.*/

  double complex tmp_trans;/**<@brief Hopping parameter.*/
  double complex tmp_J;/**<@brief Coupling constant*/

  long unsigned int is1_up;/**<@brief Mask used in the bit oeration.*/
  long unsigned int is1_down;/**<@brief Mask used in the bit oeration.*/
  long unsigned int is2_up;/**<@brief Mask used in the bit oeration.*/
  long unsigned int is2_down;/**<@brief Mask used in the bit oeration.*/

  int mode;/**<@brief multiply or expectation value.*/
  double sgn;/**<@brief Not used ???*/
  long unsigned int is1_spin;/**<@brief Mask used in the bit oeration.*/
  long unsigned int is2_spin;/**<@brief Mask used in the bit oeration.*/
  long unsigned int is3_spin;/**<@brief Mask used in the bit oeration.*/
  long unsigned int is4_spin;/**<@brief Mask used in the bit oeration.*/
  int isite1;/**<@brief Is it realy used ???*/
  int isite2;/**<@brief Is it realy used ???*/
  int isite3;/**<@brief Is it realy used ???*/
  int isite4;/**<@brief Is it realy used ???*/

  long unsigned int A_spin;/**<@brief Mask used in the bit oeration.*/
  long unsigned int B_spin;/**<@brief Mask used in the bit oeration.*/
  long unsigned int irght;/**<@brief Used for Ogata-Lin ???*/
  long unsigned int ilft;/**<@brief Used for Ogata-Lin ???*/
  long unsigned int ihfbit;/**<@brief Used for Ogata-Lin ???*/
  long unsigned int isA_spin;/**<@brief Mask used in the bit oeration.*/
  long unsigned int isB_spin;/**<@brief Mask used in the bit oeration.*/
  double complex tmp_V;/**<@brief Coupling constant*/
};/*struct LargeList*/
/**
@brief Physical quantities (Expectation value)
*/
struct PhysList {
  //double energy,doublon;
  double energy;/**<@brief Expectation value of the total energy.*/
  double doublon;/**<@brief Expectation value of the Doublon*/
  double doublon2;/**<@brief Expectation value of the Square of doublon*/
  double num;/**<@brief Expectation value of the Number of electrons*/
  double num2;/**<@brief Expectation value of the quare of the number of electrons*/
  double Sz;/**<@brief Expectation value of the Total Sz*/
  double Sz2;/**<@brief Expectation value of the Square of total Sz*/
    /*[s] For TPQ*/
  double var;/**<@brief Expectation value of the Energy variance.*/
    /*[e] For TPQ*/

    /*[s] For Full Diagonalization*/
  int eigen_num;/**<@brief Index of eigenstate used for the file name of the correlation function.*/
  double num_up;/**<@brief Expectation value of the number of up-spin electtrons*/
  double num_down;/**<@brief Expectation value of the number of down-spin electtrons*/
  double s2;/**<@brief Expectation value of the square of the total S.*/
  double *all_energy;/**<@brief [CheckList::idim_max+1] Energy for FullDiag and LOBPCG.
                     malloc in setmem_large().*/
  double *all_doublon;/**<@brief [CheckList::idim_max+1] Doublon for FullDiag and LOBPCG.
                      malloc in setmem_large().*/
  double *all_sz;/**<@brief [CheckList::idim_max+1] @f$S_z@f$ for FullDiag and LOBPCG.
                 malloc in setmem_large().*/
  double *all_s2;/**<@brief [CheckList::idim_max+1] @f$S_z^2@f$ for FullDiag and LOBPCG.
                 malloc in setmem_large().*/
  double *all_num_up;/**<@brief [CheckList::idim_max+1] Number of spin-up electrons
                     for FullDiag and LOBPCG. malloc in setmem_large().*/
  double *all_num_down;/**<@brief[CheckList::idim_max+1] Number of spin-down electrons 
                       for FullDiag and LOBPCG. malloc in setmem_large().*/
    /*[e] For Full Diagonalization*/

  double *spin_real_cor;/**<@brief Malloc, but Not used ???*/
  double *charge_real_cor;/**<@brief Malloc, but Not used ???*/
  double *loc_spin_z;/**<@brief Malloc, but Not used ???*/
  double Target_energy;/**<@brief Is it really used ???*/
  double Target_CG_energy;/**<@brief Taget energy of CG-inversed iteration (NOT LOBCG) method.*/
};/*struct PhysList*/
/**
@brief For Boost
*/
struct BoostList {
  int flgBoost;/**<@brief Flag whether use CMA algorithm.*/
  long unsigned int R0;
  long unsigned int W0;
  long unsigned int num_pivot;
  long unsigned int ishift_nspin;
  unsigned int NumarrayJ;/**<@brief */
  double complex ***arrayJ;/**<@brief */
  double complex vecB[3];/**<@brief */
  int **list_6spin_star;/**<@brief */
  int ***list_6spin_pair;/**<@brief */
};/*struct BoostList*/
/**
@brief Bind
*/
struct BindStruct {
  struct DefineList Def;/**<@brief Definision of system (Hamiltonian) etc.*/
  struct CheckList Check;/**<@brief Size of the Hilbert space*/
  struct LargeList Large;/**<@brief Variables for Matrix-Vector product*/
  struct PhysList Phys;/**<@brief Physical quantities*/
  struct BoostList Boost;/**<@brief For Boost*/
};/*struct BindStruct*/
/**
@brief Bind
*/
struct EDMainCalStruct {
  struct BindStruct Bind;/**<@brief Binded struct*/
};/*struct EDMainCalStruct*/
/**
@brief Time keeping
*/
struct TimeKeepStruct {
  time_t tstart;/**<@brief Start time*/
  time_t tnow;/**<@brief Current time*/
  time_t tend;/**<@brief End time*/
};/*struct TimeKeepStruct*/

/*global variables---------------------------------------------*/
extern struct EDMainCalStruct X;
/*-------------------------------------------------------------*/

#endif /* HPHI_STRUCT_H */
