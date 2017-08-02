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

/*=================================================================================================*/
/**
@brief What is defined ?
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
  int READ;/**<@brief It is ALWAYS 0.*/
  int WRITE;/**<@brief It is ALWAYS 0.*/

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
  unsigned int Nsize;/**<@brief Defined but NEVER used ???*/
  unsigned int Lanczos_max;/**<@brief Maximum number of iterations.*/
  int Lanczos_restart;/**<@brief Number of iterations performed in the restart computation.*/
  long int initial_iv;/**<@brief Seed of random number for initial guesss of wavefunctions.*/

  int istep;/**<@brief Index of TPQ step ???*/
  int irand;/**<@brief */
  int St;/**<@brief 0 or 1, but it affects nothing.*/

  int *LocSpn;/**<@brief [DefineList::NLocSpn] Flag (and size) of the local spin*/
  unsigned int NLocSpn;/**<@brief Number of local spins*/
  unsigned int NCond;/**<@brief Number of itinerant electrons*/
  int iFlgGeneralSpin;/**<@brief Flag for the general (Sz/=1/2) spin*/
  int iFlgSzConserved;/**<@brief Flag whether Sz is conserved.*/

  int fidx;/**<@brief Always 0, it is not used ???*/
  long unsigned int *Tpow;/**<@brief [2 * DefineList::NsiteMPI] @f$2^n@f$*/
  long unsigned int *OrgTpow;/**<@brief [2 * DefineList::NsiteMPI] @f$2^n@f$*/
  long int *SiteToBit;/**<@brief [DefineList::NsiteMPI] Similar to DefineList::Tpow.
                      For general spin.*/

  int *EDChemi;/**<@brief [DefineList::Nsite] Chemical potential*/
  unsigned int EDNChemi;/**<@brief*/
  int *EDSpinChemi;/**<@brief [DefineList::Nsite]*/

  double *EDParaChemi;/**<@brief [DefineList::Nsite]*/

    //[s] Transfer
  unsigned int NTransfer;/**<@brief Number of transfer integrals obtained by a def file.*/
  unsigned int EDNTransfer;/**<@brief Number of transfer integrals for calculation. */
  int **GeneralTransfer;/**<@brief Index of transfer integrals obtained by a def file. \n
						   Data Format [DefineList::Ntransfer][4]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  int **EDGeneralTransfer;/**<@brief Index of transfer integrals for calculation. \n
						   Data Format [DefineList::Ntransfer][4]: 0->site number i, 1-> spin index on i, 2-> site number j, 3-> spin index on j. */
  double complex *ParaGeneralTransfer;/**<@brief Value of general transfer integrals  by a def file. \n
						   Data Format [DefineList::Ntransfer].*/
  double complex *EDParaGeneralTransfer;/**<@brief Value of general transfer integrals  by a def file. \n
						   Data Format [DefineList::Ntransfer].*/
    //[e] Transfer

  int **CoulombIntra;/**< [DefineList::NCoulombIntra][1]*/
  unsigned int NCoulombIntra;/**< */
  double *ParaCoulombIntra;/**< [DefineList::NCoulombIntra]*/

  int **CoulombInter;/**< [DefineList::NCoulombInter][2]*/
  unsigned int NCoulombInter;/**<@brief */
  double *ParaCoulombInter;/**<@brief [DefineList::NCoulombInter]*/

  int **HundCoupling;/**<@brief */
  unsigned int NHundCoupling;/**<@brief */
  double *ParaHundCoupling;/**<@brief */

  int **PairHopping;/**<@brief */
  unsigned int NPairHopping;/**<@brief */
  double *ParaPairHopping;/**<@brief */

  int **ExchangeCoupling;/**<@brief */
  unsigned int NExchangeCoupling;/**<@brief */
  double *ParaExchangeCoupling;/**<@brief */

  unsigned int NIsingCoupling;/**<@brief */

  int **PairLiftCoupling;/**<@brief */
  unsigned int NPairLiftCoupling;/**<@brief */
  double *ParaPairLiftCoupling;/**<@brief */


    //[s] For InterAll
  int **InterAll;/**<@brief [DefineList::NinterAll][8] Interacted quartet */
  int **InterAll_OffDiagonal;/**<@brief [DefineList::NinterAll_OffDiagonal][8] Interacted quartet */
  int **InterAll_Diagonal;/**<@brief [DefineList::NinterAll_Diagonal][4] Interacted quartet */
  unsigned int NInterAll;/**<@brief */
  unsigned int NInterAll_Diagonal;/**<@brief */
  unsigned int NInterAll_OffDiagonal;/**<@brief */
  double complex *ParaInterAll;/**<@brief */
  double *ParaInterAll_Diagonal;/**<@brief */
  double complex *ParaInterAll_OffDiagonal;
    //[e] For InterAll

  int **CisAjt;/**<@brief */
  unsigned int NCisAjt;/**<@brief */

  int **CisAjtCkuAlvDC;/**<@brief */
  unsigned int NCisAjtCkuAlvDC;/**<@brief */

  int **SingleExcitationOperator;/**<@brief */
  unsigned int NSingleExcitationOperator;/**<@brief */
  double complex *ParaSingleExcitationOperator;/**<@brief */

  int **PairExcitationOperator;/**<@brief */
  unsigned int NPairExcitationOperator;/**<@brief */
  double complex *ParaPairExcitationOperator;/**<@brief */
  
  int iCalcType;/**<@brief An integer for selecting calculation type. 0:Lanczos, 1:TPQCalc, 2:FullDiag.*/
  int iCalcEigenVec;/**<@brief An integer for selecting method to calculate eigenvectors. 0:Lanczos+CG, 1: Lanczos. default value is set as 0 in readdef.c*/
  int iInitialVecType;/**<@brief An integer for setting a type of inital vectors. 0:complex type, 1: real type. default value is set as 0 in readdef.c*/
  int iFlgFiniteTemperature;/**<@brief */
  int iCalcModel;/**<@brief An integer for selecting calculation model. 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC, 6:HubbardNConserved*/
  int iOutputMode;/**<@brief An integer for selecting output mode. 0: OneBodyG and TwoBodyG. 1: OneBodyG and TwoBodyG and correlations for charge and spin.*/
  int iOutputEigenVec;/**<@brief An integer for selecting output an eigenvector. 0: no output, 1:output.*/
  int iInputEigenVec;/**<@brief An integer for selecting output an eigenvector. 0: no input, 1:input*/
  int iOutputHam;/**< An integer for selecting output a Hamiltonian. 0: no output, 1:output*/
  int iInputHam;/**< An integer for selecting input a Hamiltonian. 0: no input, 1:input*/


    //[s] For Spectrum
  double complex dcOmegaMax;/**<@brief */
  double complex dcOmegaMin;/**<@brief */
  double complex dcOmegaOrg;/**<@brief */
  int iNOmega;/**<@brief */
  int iFlgSpecOmegaMax;/**<@brief */
  int iFlgSpecOmegaMin;/**<@brief */
  int iFlgSpecOmegaOrg;/**<@brief */
  int iFlgCalcSpec;/**<@brief */
  int iFlagListModified;/**<@brief */
    /**<@brief An integer for selecting calculation type. 0:Lanczos, 1:TPQCalc, 2:FullDiag.*/

    //[e] For Spectrum

  int iReStart;/**< An integer for restarting output a Hamiltonian.
     - 0: not restart
     - 1:restart (output restart vector),
     - 2: restart (input and output restart vector) */
  int iFlgMPI;/**<@brief MPI mode
    - 0: butterfly
    - 1: Parallel Interaction [to be supported]
    */
};/*struct DefineList*/

struct CheckList {
  unsigned long int idim_max;/**<@brief */
  unsigned long int idim_maxMPI;/**<@brief */
  unsigned long int idim_maxOrg;/**<@brief calcspectrum*/
  unsigned long int idim_maxMPIOrg;/**<@brief */
  unsigned long int sdim;/**<@brief */
  double max_mem;/**<@brief */
};/*struct CheckList*/

struct LargeList {
  double complex prdct;/**<@brief */
  int itr;/**<@brief */
  long int iv;/**<@brief */
  long int i_max;/**<@brief */
  long int SizeOflist_2_1;/**<@brief */
  long int SizeOflist_2_2;/**<@brief */
  long int SizeOflistjb;/**<@brief */

  double complex tmp_trans;/**<@brief */
  double complex tmp_J;/**<@brief */

  long unsigned int is1_up;/**<@brief */
  long unsigned int is1_down;/**<@brief */
  long unsigned int is2_up;/**<@brief */
  long unsigned int is2_down;/**<@brief */

  int mode;/**<@brief */
  double sgn;/**<@brief */
  long unsigned int is1_spin;/**<@brief */
  long unsigned int is2_spin;/**<@brief */
  long unsigned int is3_spin;/**<@brief */
  long unsigned int is4_spin;/**<@brief */
  int isite1;/**<@brief */
  int isite2;/**<@brief */
  int isite3;/**<@brief */
  int isite4;/**<@brief */

  long unsigned int A_spin;/**<@brief */
  long unsigned int B_spin;/**<@brief */
  long unsigned int irght;/**<@brief */
  long unsigned int ilft;/**<@brief */
  long unsigned int ihfbit;/**<@brief */
  long unsigned int isA_spin;/**<@brief */
  long unsigned int isB_spin;/**<@brief */
  double complex tmp_V;/**<@brief */
};/*struct LargeList*/

struct PhysList {
  //double energy,doublon;
  double energy;/**<@brief */
  double doublon, doublon2;/**<@brief */
  double num, num2;/**<@brief */
  double Sz, Sz2;/**<@brief */
    /*[s] For TPQ*/
  double var;/**<@brief */
    /*[e] For TPQ*/

    /*[s] For Full Diagonalization*/
  int eigen_num;/**<@brief */
  double num_up;/**<@brief */
  double num_down;/**<@brief */
  double s2;/**<@brief */
  double sz;/**<@brief */
  double *all_energy;/**<@brief */
  double *all_doublon;/**<@brief */
  double *all_sz;/**<@brief */
  double *all_s2;/**<@brief */
  double *all_num_up;/**<@brief */
  double *all_num_down;/**<@brief */
    /*[e] For Full Diagonalization*/

  double *spin_real_cor;/**<@brief */
  double *charge_real_cor;/**<@brief */
  double *loc_spin_z;/**<@brief */
  double Target_energy;/**<@brief */
  double Target_CG_energy;/**<@brief */
};/*struct PhysList*/

//For Boost
struct BoostList {
  int flgBoost;/**<@brief */
  long unsigned int R0, W0, num_pivot, ishift_nspin;
  unsigned int NumarrayJ;/**<@brief */
  double complex ***arrayJ;/**<@brief */
  double complex vecB[3];/**<@brief */
  int **list_6spin_star;/**<@brief */
  int ***list_6spin_pair;/**<@brief */
};/*struct BoostList*/

/*=================================================================================================*/
struct BindStruct {
  struct DefineList Def;/**<@brief */
  struct CheckList Check;/**<@brief */
  struct LargeList Large;/**<@brief */
  struct PhysList Phys;/**<@brief */
  struct BoostList Boost;/**<@brief */
};/*struct BindStruct*/
/*=================================================================================================*/
struct EDMainCalStruct {
  struct BindStruct Bind;/**<@brief */
};/*struct EDMainCalStruct*/

struct TimeKeepStruct {
  time_t tstart;/**<@brief */
  time_t tnow;/**<@brief */
  time_t tend;/**<@brief */
};/*struct TimeKeepStruct*/

/*global variables---------------------------------------------*/
struct EDMainCalStruct X;
/*-------------------------------------------------------------*/

#endif /* HPHI_STRUCT_H */
