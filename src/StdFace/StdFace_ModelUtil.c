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
@brief Various utility for constructing models
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "StdFace_vals.h"
#include <string.h>
#ifdef MPI
#include <mpi.h>
#endif

/**
@brief MPI Abortation wrapper
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_exit(int errorcode//!< [in]
)
{
  int ierr = 0;
  fflush(stdout);
  fflush(stderr);
#ifdef MPI
  fprintf(stdout, "\n\n #######  You DO NOT have to WORRY about the following MPI-ERROR MESSAGE.  #######\n\n");
  ierr = MPI_Abort(MPI_COMM_WORLD, errorcode);
  ierr = MPI_Finalize();
#endif
  if (ierr != 0) fprintf(stderr, "\n  MPI_Finalize() = %d\n\n", ierr);
  exit(errorcode);
}
/**
@brief Add transfer to the list
set StdIntList::trans and StdIntList::transindx and
increment StdIntList::ntrans
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_trans(
  struct StdIntList *StdI,//!<[inout]
  double complex trans0,//!<[in] Hopping integral @f$t, mu@f$, etc.
  int isite,//!<[in] @f$i@f$ for @f$c_{i \sigma}^\dagger@f$
  int ispin,//!<[in] @f$\sigma@f$ for @f$c_{i \sigma}^\dagger@f$
  int jsite,//!<[in] @f$j@f$ for @f$c_{j \sigma'}@f$
  int jspin//!<[in] @f$\sigma'@f$ for @f$c_{j \sigma'}@f$
)
{
  StdI->trans[StdI->ntrans] = trans0;
  StdI->transindx[StdI->ntrans][0] = isite;
  StdI->transindx[StdI->ntrans][1] = ispin;
  StdI->transindx[StdI->ntrans][2] = jsite; 
  StdI->transindx[StdI->ntrans][3] = jspin;
  StdI->ntrans = StdI->ntrans + 1;
}/*void StdFace_trans*/
/**
@brief Add Hopping for the both spin
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Hopping(
struct StdIntList *StdI,//!<[inout]
  double complex trans0,//!<[in] Hopping integral @f$t@f$
  int isite,//!<[in] @f$i@f$ for @f$c_{i \sigma}^\dagger@f$
  int jsite//!<[in] @f$j@f$ for @f$c_{j \sigma}@f$
)
{
  int ispin;
  /**@brief
   Both @f$c_{i \sigma}^\daggerc_{j \sigma}@f$ and
  @f$c_{j \sigma}^\daggerc_{i \sigma}@f$ for every spin channel
  (@f$\sigma@f$) is specified
  */
  for (ispin = 0; ispin < 2; ispin++) {
    StdFace_trans(StdI, trans0, jsite, ispin, isite, ispin);
    StdFace_trans(StdI, conj(trans0), isite, ispin, jsite, ispin);
  }/*for (ispin = 0; ispin < 2; ispin++)*/
}/*void StdFace_Hopping*/
/**
@brief Add intra-Coulomb, magnetic field, chemical potential for the
itenerant electron
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_HubbardLocal(
  struct StdIntList *StdI,//!<[inout]
  double mu0,//!<[in] Chemical potential
  double h0,//!<[in] Longitudinal magnetic feild
  double Gamma0,//!<[in] Transvers magnetic feild
  double U0,//!<[in] Intra-site Coulomb potential
  int isite//!<[in] i for @f$c_{i \sigma}^\dagger@f$
)
{
  StdFace_trans(StdI, mu0 + 0.5 * h0, isite, 0, isite, 0);
  StdFace_trans(StdI, mu0 - 0.5 * h0, isite, 1, isite, 1);
  StdFace_trans(StdI, -0.5 * Gamma0, isite, 1, isite, 0);
  StdFace_trans(StdI, -0.5 * Gamma0, isite, 0, isite, 1);
  /**@brief
  Set StdIntList::Cintra and StdIntList::CintraIndx
  with the input argument and increase the number
  of them (StdIntList::NCintra).
  */
  StdI->Cintra[StdI->NCintra] = U0;
  StdI->CintraIndx[StdI->NCintra][0] = isite;
  StdI->NCintra += 1;
}/*void StdFace_HubbardLocal*/
/**
@brief Add longitudinal and transvars magnetic field to the list
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_MagField(
  struct StdIntList *StdI,//!<[inout]
  int S2,//!<[in] Spin moment in @f$i@f$ site
  double h,//!<[in] Longitudinal magnetic field @f$h@f$
  double Gamma,//!<[in] Transvars magnetic field @f$h@f$
  int isite//!<[in] @f$i@f$ for @f$c_{i \sigma}^\dagger@f$
)
{
  int ispin;
  double S, Sz;
 
  S = (double)S2 * 0.5;
  /**@brief
   Use Bogoliubov representation.
  */
  for (ispin = 0; ispin <= S2; ispin++){
    /**@brief
    Londitudinal part
    @f[
     \sum_{\sigma = -S}^{S} -h\sigma c_{i \sigma}^\dagger c_{i \sigma}
    @f]
    */
    Sz = (double)ispin - S;
    StdFace_trans(StdI, -h * Sz, isite, ispin, isite, ispin);
    /**@brief
    Transvars part
    @f[
    -\Gamma \frac{S_i^+ + S_i^-}{2} =
    \sum_{\sigma = -S}^{S-1} -\frac{\Gamma}{2}
    \sqrt{S(S+1) - \sigma(\sigma+1)}
    (\sigma c_{i \sigma+ 1}^\dagger c_{i \sigma} + 
    \sigma c_{i \sigma}^\dagger c_{i \sigma+1})
    @f]
    */
    if (ispin < S2) {
      StdFace_trans(StdI, -0.5 * Gamma * sqrt(S*(S + 1.0) - Sz*(Sz + 1.0)),
        isite, ispin + 1, isite, ispin);
      StdFace_trans(StdI, -0.5 * Gamma * sqrt(S*(S + 1.0) - Sz*(Sz + 1.0)),
        isite, ispin, isite, ispin + 1);
    }/*if (ispin < S2)*/
  }/*for (ispin = 0; ispin <= S2; ispin++)*/
}/*void StdFace_MagField*/
/**
@brief Add interaction (InterAll) to the list
Set StdIntList::intr and StdIntList::intrindx and
increase the number of that (StdIntList::nintr).
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_intr(
  struct StdIntList *StdI,//!<[inout]
  double complex intr0,//!<[in] Interaction @f$U, V, J@f$, etc.
  int site1,//!<[in] @f$i_1@f$ for @f$c_{i_1 \sigma_1}^\dagger@f$
  int spin1,//!<[in] @f$sigma1_1@f$ for @f$c_{i_1 \sigma_1}^\dagger@f$
  int site2,//!<[in] @f$i_2@f$ for @f$c_{i_2 \sigma_2}@f$
  int spin2,//!<[in] @f$sigma1_2@f$ for @f$c_{i_2 \sigma_2}@f$
  int site3,//!<[in] @f$i_3@f$ for @f$c_{i_3 \sigma_3}^\dagger@f$
  int spin3,//!<[in] @f$sigma1_3@f$ for @f$c_{i_3 \sigma_3}^\dagger@f$
  int site4,//!<[in] @f$i_2@f$ for @f$c_{i_2 \sigma_2}@f$
  int spin4//!<[in] @f$sigma1_2@f$ for @f$c_{i_2 \sigma_2}@f$
)
{
  StdI->intr[StdI->nintr] = intr0;
  StdI->intrindx[StdI->nintr][0] = site1; StdI->intrindx[StdI->nintr][1] = spin1;
  StdI->intrindx[StdI->nintr][2] = site2; StdI->intrindx[StdI->nintr][3] = spin2;
  StdI->intrindx[StdI->nintr][4] = site3; StdI->intrindx[StdI->nintr][5] = spin3;
  StdI->intrindx[StdI->nintr][6] = site4; StdI->intrindx[StdI->nintr][7] = spin4;
  StdI->nintr = StdI->nintr + 1;
}/*void StdFace_intr*/
/**
@brief Treat J as a 3*3 matrix [(6S + 1)*(6S' + 1) interactions]
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_GeneralJ(
struct StdIntList *StdI,//!<[inout]
  double J[3][3],//!<[in] The Spin interaction @f$J_x, J_{xy}@f$, ...
  int Si2,//!<[in] Spin moment in @f$i@f$ site
  int Sj2,//!<[in] Spin moment in @f$j@f$ site
  int isite,//!<[in] @f$i@f$ of @f$S_i@f$
  int jsite//!<[in] @f$j@f$ of @f$S_j@f$
)
{
  int ispin, jspin, ZGeneral, ExGeneral;
  double Si, Sj, Siz, Sjz;
  double complex intr0;
  /**@brief
   If both spins are S=1/2 ...
  */
  ZGeneral = 1;
  ExGeneral = 1;
  if (Si2 == 1 || Sj2 == 1) {
    /**@brief
    Set StdIntList::Hund, StdIntList::HundIndx, StdIntList::Cinter, 
    StdIntList::CinterIndx for @f$S_{iz} S_{jz}@f$ term,
    and increase the number of them 
    (StdIntList::NHund, StdIntList::NCinter). And ...
    */
    ZGeneral = 0;

    StdI->Hund[StdI->NHund] = -0.5 * J[2][2];
    StdI->HundIndx[StdI->NHund][0] = isite;
    StdI->HundIndx[StdI->NHund][1] = jsite;
    StdI->NHund += 1;

    StdI->Cinter[StdI->NCinter] = -0.25 * J[2][2];
    StdI->CinterIndx[StdI->NCinter][0] = isite;
    StdI->CinterIndx[StdI->NCinter][1] = jsite;
    StdI->NCinter += 1;

    if (fabs(J[0][1]) < 0.000001 && fabs(J[1][0]) < 0.000001
#if defined(_mVMC)
      && abs(J[0][0] - J[1][1]) < 0.000001
#endif
      ) {
      /**@brief
       If there is no off-diagonal term, set 
       StdIntList::Ex, StdIntList::ExIndx, StdIntList::PairLift, 
       StdIntList::PLIndx and increase the number of them
       (StdIntList::NEx, StdIntList::NPairLift).
      */
      ExGeneral = 0;

#if defined(_mVMC)
      StdI->Ex[StdI->NEx] = - 0.25 * (J[0][0] + J[1][1]);
#else
      if (strcmp(StdI->model, "kondo") == 0)
        StdI->Ex[StdI->NEx] = -0.25 * (J[0][0] + J[1][1]);
      else
        StdI->Ex[StdI->NEx] = 0.25 * (J[0][0] + J[1][1]);
#endif
      StdI->ExIndx[StdI->NEx][0] = isite;
      StdI->ExIndx[StdI->NEx][1] = jsite;
      StdI->NEx += 1;

      StdI->PairLift[StdI->NPairLift] = 0.25 * (J[0][0] - J[1][1]);
      StdI->PLIndx[StdI->NPairLift][0] = isite;
      StdI->PLIndx[StdI->NPairLift][1] = jsite;
      StdI->NPairLift += 1;
    }/*if (fabs(J[0][1]) < 0.000001 && fabs(J[1][0]) < 0.000001)*/
  }
  /**@brief
   If S != 1/2 or there is an off-diagonal interaction, 
   use the InterAll format.
  */
  Si = 0.5 * (double)Si2;
  Sj = 0.5 * (double)Sj2;

  for (ispin = 0; ispin <= Si2; ispin++) {
    Siz = (double)ispin - Si;
    for (jspin = 0; jspin <= Sj2; jspin++) {
      Sjz = (double)jspin - Sj;
      /**@brief (1)
       @f[
       J_z S_{i z} * S_{j z} = J_z \sum_{\sigma, \sigma' = -S}^S
       \sigma \sigma' c_{i\sigma}^\dagger c_{i\sigma}
       c_{j\sigma'}^\dagger c_{j\sigma'}
       @f]
      */
      if (ZGeneral == 1) {
        intr0 = J[2][2] * Siz * Sjz;
        StdFace_intr(StdI, intr0,
          isite, ispin, isite, ispin, jsite, jspin, jsite, jspin);
      }
      /**@brief (2)
      @f[
      I S_i^+ S_j^- + I^* S_j^+ S_i^- = 
      \sum_{\sigma, \sigma' = -S}^{S-1}
      \sqrt{S_i(S_i+1) - \sigma(\sigma+1)}\sqrt{S_j(S_j+1) - \sigma'(\sigma'+1)}
      (I c_{i\sigma+1}^\dagger c_{i\sigma} c_{j\sigma'}^\dagger c_{j\sigma'+1} +
      I^* c_{i\sigma}^\dagger c_{i\sigma+1} c_{j\sigma'+1}^\dagger c_{j\sigma'})
      \\
      I \equiv \frac{J_x + J_y + i(J_{xy} - J_{yx})}{4}
      @f]
      */
      if ((ispin < Si2 && jspin < Sj2) && ExGeneral == 1) {
        intr0 = 0.25 * (J[0][0] + J[1][1] + I*(J[0][1] - J[1][0]))
          * sqrt(Si * (Si + 1.0) - Siz * (Siz + 1.0))
          * sqrt(Sj * (Sj + 1.0) - Sjz * (Sjz + 1.0));
        StdFace_intr(StdI, intr0,
          isite, ispin + 1, isite, ispin, jsite, jspin, jsite, jspin + 1);
        StdFace_intr(StdI, conj(intr0),
          isite, ispin, isite, ispin + 1, jsite, jspin + 1, jsite, jspin);
      }
      /**@brief (3)
      @f[
      I S_i^+ S_j^+ + I^* S_j^- S_i^-
      = \sum_{\sigma, \sigma' = -S}^{S-1}
      \sqrt{S_i(S_i+1) - \sigma(\sigma+1)}\sqrt{S_j(S_j+1) - \sigma'(\sigma'+1)}
      (I c_{i\sigma+1}^\dagger c_{i\sigma} c_{j\sigma'+1}^\dagger c_{j\sigma'} +
      I^* c_{i\sigma}^\dagger c_{i\sigma+1} c_{j\sigma'}^\dagger c_{j\sigma'+1})
      \\
      I \equiv \frac{J_x - J_y - i(J_{xy} + J_{yx})}{4}
      @f]
      */
      if ((ispin < Si2 && jspin < Sj2) && ExGeneral == 1) {
        intr0 = 0.5 * 0.5 * (J[0][0] - J[1][1] - I*(J[0][1] + J[1][0]))
          * sqrt(Si * (Si + 1.0) - Siz * (Siz + 1.0))
          * sqrt(Sj * (Sj + 1.0) - Sjz * (Sjz + 1.0));
        StdFace_intr(StdI, intr0,
          isite, ispin + 1, isite, ispin, jsite, jspin + 1, jsite, jspin);
        StdFace_intr(StdI, conj(intr0),
          isite, ispin, isite, ispin + 1, jsite, jspin, jsite, jspin + 1);
      }
      /**@brief (4)
      @f[
      I S_i^+ S_{j z} + I^* S_{j z} S_i^-= 
      \sum_{\sigma=-S}^{S-1} \sum_{\sigma' = -S}^{S} \sqrt{S_i(S_i+1) - \sigma(\sigma+1)}
      (I c_{i\sigma+1}^\dagger c_{i\sigma} c_{j\sigma'}^\dagger c_{j\sigma'} +
      I^* c_{j\sigma'}^\dagger c_{j\sigma'} c_{i\sigma}^\dagger c_{i\sigma+1})
      \\
      I \equiv \frac{J_{xz} - i J_{yz}}{2}
      @f]
      */
      if (ispin < Si2) {
        intr0 = 0.5 * (J[0][2] - I * J[1][2]) * sqrt(Si * (Si + 1.0) - Siz * (Siz + 1.0)) * Sjz;
        StdFace_intr(StdI, intr0,
          isite, ispin + 1, isite, ispin, jsite, jspin, jsite, jspin);
        StdFace_intr(StdI, conj(intr0),
          jsite, jspin, jsite, jspin, isite, ispin, isite, ispin + 1);
      }/*if (ispin < Si2)*/
      /**@brief (5)
      @f[
       I S_{i z} S_j^+ + I^* S_j^- S_{i z} =
       \sum_{\sigma=-S}^{S} \sum_{\sigma' = -S}^{S-1} \sqrt{S_j(S_j+1) - \sigma'(\sigma'+1)}
       (I c_{i\sigma}^\dagger c_{i\sigma} c_{j\sigma'+1}^\dagger c_{j\sigma'} +
       I^* c_{j\sigma'}^\dagger c_{j\sigma'+1} c_{i\sigma}^\dagger c_{i\sigma})
       \\
       I \equiv \frac{J_{zx} - i J_{zy}}{2}
       @f]
      */
      if (jspin < Sj2) {
        intr0 = 0.5 * (J[2][0] - I * J[2][1]) * Siz * sqrt(Sj * (Sj + 1.0) - Sjz * (Sjz + 1.0));
        StdFace_intr(StdI, intr0,
          isite, ispin, isite, ispin, jsite, jspin + 1, jsite, jspin);
        StdFace_intr(StdI, conj(intr0),
          jsite, jspin, jsite, jspin + 1, isite, ispin, isite, ispin);
      }/*if (jspin < Sj2)*/
    }/*for (jspin = 0; jspin <= Sj2; jspin++)*/
  }/*for (ispin = 0; ispin <= Si2; ispin++)*/
}/*StdFace_GeneralJ*/
/**
@brief Add onsite/offsite Coulomb term to the list
StdIntList::Cinter and StdIntList::CinterIndx,
and increase the number of them (StdIntList::NCinter).
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Coulomb(
struct StdIntList *StdI,//!<[inout]
  double V,//!<[in] Coulomb integral U, V, etc.
  int isite,//!<[in] i of n_i
  int jsite//!<[in] j of n_j
)
{
  StdI->Cinter[StdI->NCinter] = V;
  StdI->CinterIndx[StdI->NCinter][0] = isite;
  StdI->CinterIndx[StdI->NCinter][1] = jsite;
  StdI->NCinter += 1;
}/*void StdFace_Coulomb*/
/**
@brief Print a valiable (real) read from the input file
if it is not specified in the input file (=NaN), 
set the default value.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_PrintVal_d(
  char* valname,//!<[in] Name of the valiable
  double *val,//!<[inout] Valiable to be set 
  double val0//!<[in] The default value
)
{
  if (isnan(*val) == 1) {
    *val = val0;
    fprintf(stdout, "  %15s = %-10.5f  ######  DEFAULT VALUE IS USED  ######\n", valname, *val);
  }
  else fprintf(stdout, "  %15s = %-10.5f\n", valname, *val);
}/*void StdFace_PrintVal_d*/
/**
@brief Print a valiable (real) read from the input file
if it is not specified in the input file (=NaN),
set the default value.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_PrintVal_dd(
  char* valname,//!<[in] Name of the valiable
  double *val,//!<[inout] Valiable to be set
  double val0,//!<[in] The primary default value, possible not to be specified
  double val1//!<[in] The secondary default value
)
{
  if (isnan(*val) == 1) {
    /**@brief
    If the primary default value (val0) is not specified, use secondary default value
    */
    if (isnan(val0) == 1) *val = val1;
    else *val = val0;
    fprintf(stdout, "  %15s = %-10.5f  ######  DEFAULT VALUE IS USED  ######\n", valname, *val);
  }
  else fprintf(stdout, "  %15s = %-10.5f\n", valname, *val);
}/*void StdFace_PrintVal_dd*/
/**
@brief Print a valiable (complex) read from the input file
if it is not specified in the input file (=NaN),
set the default value.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_PrintVal_c(
  char* valname,//!<[in] Name of the valiable
  double complex *val,//!<[inout] Valiable to be set
  double complex val0//!<[in] The default value
)
{
  if (isnan(creal(*val)) == 1) {
    *val = val0;
    fprintf(stdout, "  %15s = %-10.5f %-10.5f  ######  DEFAULT VALUE IS USED  ######\n", valname, creal(*val), cimag(*val));
  }
  else fprintf(stdout, "  %15s = %-10.5f %-10.5f\n", valname, creal(*val), cimag(*val));
}/*void StdFace_PrintVal_c*/
/**
@brief Print a valiable (integer) read from the input file
if it is not specified in the input file (=2147483647, the upper limt of Int)
set the default value.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_PrintVal_i(
  char* valname,//!<[in] Name of the valiable
  int *val,//!<[inout] Valiable to be set
  int val0//!<[in] The default value
)
{
  int NaN_i = 2147483647;/*The upper limt of Int*/

  if (*val == NaN_i) {
    *val = val0;
    fprintf(stdout, "  %15s = %-10d  ######  DEFAULT VALUE IS USED  ######\n", valname, *val);
  }
  else fprintf(stdout, "  %15s = %-10d\n", valname, *val);
}/*void StdFace_PrintVal_i*/
/**
@brief Stop HPhi if a variable (real) not used is specified
in the input file (!=NaN).
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_NotUsed_d(
  char* valname,//!<[in] Name of the valiable
  double val//!<[in]
)
{
  if (isnan(val) == 0) {
    fprintf(stdout, "\n Check !  %s is SPECIFIED but will NOT be USED. \n", valname);
    fprintf(stdout, "            Please COMMENT-OUT this line \n");
    fprintf(stdout, "            or check this input is REALLY APPROPRIATE for your purpose ! \n\n");
    StdFace_exit(-1);
  }
}/*void StdFace_NotUsed_d*/
/**
@brief Stop HPhi if a variable (complex) not used is specified
in the input file (!=NaN).
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_NotUsed_c(
  char* valname,//!<[in] Name of the valiable
  double complex val//!<[in]
)
{
  if (isnan(creal(val)) == 0) {
    fprintf(stdout, "\n Check !  %s is SPECIFIED but will NOT be USED. \n", valname);
    fprintf(stdout, "            Please COMMENT-OUT this line \n");
    fprintf(stdout, "            or check this input is REALLY APPROPRIATE for your purpose ! \n\n");
    StdFace_exit(-1);
  }
}/*void StdFace_NotUsed_c*/
/**
@brief Stop HPhi if variables (real) not used is specified
in the input file (!=NaN).
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_NotUsed_J(
  char* valname,//!<[in] Name of the valiable*/,
  double JAll,//!<[in]*/,
  double J[3][3]//!<[in]
)
{
  int i1, i2;
  char Jname[3][3][10];

  sprintf(Jname[0][0], "%sx", valname);
  sprintf(Jname[0][1], "%sxy", valname);
  sprintf(Jname[0][2], "%sxz", valname);
  sprintf(Jname[1][0], "%syx", valname);
  sprintf(Jname[1][1], "%sy", valname);
  sprintf(Jname[1][2], "%syz", valname);
  sprintf(Jname[2][0], "%szx", valname);
  sprintf(Jname[2][1], "%szy", valname);
  sprintf(Jname[2][2], "%sz", valname);

  StdFace_NotUsed_d(valname, JAll);

  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      StdFace_NotUsed_d(Jname[i1][i2], J[i1][i2]);
    }/*for (j = 0; j < 3; j++)*/
  }/*for (i = 0; i < 3; i++)*/
}/*void StdFace_NotUsed_J*/
/**
@brief Stop HPhi if a variable (integer) not used is specified
in the input file (!=2147483647, the upper limt of Int).
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_NotUsed_i(
  char* valname,//!<[in] Name of the valiable
  int val//!<[in]
)
{
  int NaN_i = 2147483647;

  if (val != NaN_i) {
    fprintf(stdout, "\n Check !  %s is SPECIFIED but will NOT be USED. \n", valname);
    fprintf(stdout, "            Please COMMENT-OUT this line \n");
    fprintf(stdout, "            or check this input is REALLY APPROPRIATE for your purpose ! \n\n");
    StdFace_exit(-1);
  }
}/*void StdFace_NotUsed_i*/
/**
@brief Stop HPhi if a variable (integer) which must be specified
is absent in the input file (=2147483647, the upper limt of Int).
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_RequiredVal_i(
  char* valname,//!<[in] Name of the valiable
  int val//!<[in]
)
{
  int NaN_i = 2147483647;

  if (val == NaN_i){
    fprintf(stdout, "ERROR ! %s is NOT specified !\n", valname);
    StdFace_exit(-1);
  }
  else fprintf(stdout, "  %15s = %-3d\n", valname, val);
}/*void StdFace_RequiredVal_i*/
/**
@brief Move a site into the original supercell if it is outside the 
original supercell.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StdFace_FoldSite(
  struct StdIntList *StdI,//!<[inout]
  int iCellV[3],//!<[in] The fractional coordinate of a site
  int nBox[3], //!<[out] the index of supercell
  int iCellV_fold[3]/**<[out] The fractional coordinate of a site 
                    which is moved into the original cell*/
)
{
  int ii, jj, iCellV_frac[3];
  /**@brief
  (1) Transform to fractional coordinate (times NCell).
  */
  for (ii = 0; ii < 3; ii++) {
    iCellV_frac[ii] = 0;
    for (jj = 0; jj < 3; jj++)iCellV_frac[ii] += StdI->rbox[ii][jj] * iCellV[jj];
  }
  /**@brief
  (2) Search which supercell contains this cell.
  */
  for (ii = 0; ii < 3; ii++)
    nBox[ii] = (iCellV_frac[ii] + StdI->NCell * 1000) / StdI->NCell - 1000;
  /**@brief
  (3) Fractional coordinate (times NCell) in the original supercell
  */
  for (ii = 0; ii < 3; ii++)
    iCellV_frac[ii] -= StdI->NCell*(nBox[ii]);
  /**/
  for (ii = 0; ii < 3; ii++) {
    iCellV_fold[ii] = 0;
    for (jj = 0; jj < 3; jj++) iCellV_fold[ii] += StdI->box[jj][ii] * iCellV_frac[jj];
    iCellV_fold[ii] = (iCellV_fold[ii] + StdI->NCell * 1000) / StdI->NCell - 1000;
  }/*for (ii = 0; ii < 3; ii++)*/
}/*static void StdFace_FoldSite*/
/**
@brief Initialize the super-cell where simulation is performed.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_InitSite(
  struct StdIntList *StdI,//!<[inout]
  FILE *fp,//!<[in] File pointer to lattice.gp
  int dim//!<[in] dimension of system, if = 2, print lattice.gp
)
{
  int bound[3][2], edge, ii, jj;
  int ipos;
  int nBox[3], iCellV_fold[3], iCellV[3];
  double pos[4][2], xmin, xmax/*, offset[2], scale*/;

  fprintf(stdout, "\n  @ Super-Lattice setting\n\n");
  /**@brief
  (1) Check Input parameters about the shape of super-cell
  */
  if (
    (StdI->L != StdI->NaN_i || StdI->W != StdI->NaN_i || StdI->Height != StdI->NaN_i)
    && 
    (StdI->box[0][0] != StdI->NaN_i || StdI->box[0][1] != StdI->NaN_i || StdI->box[0][2] != StdI->NaN_i ||
     StdI->box[1][0] != StdI->NaN_i || StdI->box[1][1] != StdI->NaN_i || StdI->box[1][2] != StdI->NaN_i ||
     StdI->box[2][0] != StdI->NaN_i || StdI->box[2][1] != StdI->NaN_i || StdI->box[2][2] != StdI->NaN_i)
    )
  {
    fprintf(stdout, "\nERROR ! (L, W, Height) and (a0W, ..., a2H) conflict !\n\n");
    StdFace_exit(-1);
  }
  else if (StdI->L != StdI->NaN_i || StdI->W != StdI->NaN_i || StdI->Height != StdI->NaN_i)
  {
    StdFace_PrintVal_i("L", &StdI->L, 1);
    StdFace_PrintVal_i("W", &StdI->W, 1);
    StdFace_PrintVal_i("Height", &StdI->Height, 1);
    for (ii = 0; ii < 3; ii++) for (jj = 0; jj < 3; jj++)
      StdI->box[ii][jj] = 0;
    StdI->box[0][0] = StdI->W;
    StdI->box[1][1] = StdI->L;
    StdI->box[2][2] = StdI->Height;
  }
  else
  {
    StdFace_PrintVal_i("a0W", &StdI->box[0][0], 1);
    StdFace_PrintVal_i("a0L", &StdI->box[0][1], 0);
    StdFace_PrintVal_i("a0H", &StdI->box[0][2], 0);
    StdFace_PrintVal_i("a1W", &StdI->box[1][0], 0);
    StdFace_PrintVal_i("a1L", &StdI->box[1][1], 1);
    StdFace_PrintVal_i("a1H", &StdI->box[1][2], 0);
    StdFace_PrintVal_i("a2W", &StdI->box[2][0], 0);
    StdFace_PrintVal_i("a2L", &StdI->box[2][1], 0);
    StdFace_PrintVal_i("a2H", &StdI->box[2][2], 1);
  }
  /*
   Parameters for the 3D system will not used.
  */
  if (dim == 2) {
    StdI->direct[0][2] = 0.0;
    StdI->direct[1][2] = 0.0;
    StdI->direct[2][0] = 0.0;
    StdI->direct[2][1] = 0.0;
    StdI->direct[2][2] = 1.0;
  }
  /**@brief
  (2) Define the phase factor at each boundary.
  Set anti-period flag (StdIntList::AntiPeriod).
  */
  if (dim == 2) StdI->phase[2] = 0.0;
  for (ii = 0; ii < 3; ii++) {
    StdI->ExpPhase[ii] = cos(StdI->pi180 * StdI->phase[ii]) + I*sin(StdI->pi180 * StdI->phase[ii]);
    if (cabs(StdI->ExpPhase[ii] + 1.0) < 0.000001) StdI->AntiPeriod[ii] = 1;
    else StdI->AntiPeriod[ii] = 0;
  }
  /**@brief
  (3) Malloc StdIntList::tau, intrinsic structure of unit-cell
  */
  StdI->tau = (double **)malloc(sizeof(double*) * StdI->NsiteUC);
  for (ii = 0; ii < StdI->NsiteUC; ii++) {
    StdI->tau[ii] = (double *)malloc(sizeof(double) * 3);
  }
  /**@brief
  (4) Calculate reciprocal lattice vectors (times StdIntList::NCell, the number of unit-cells)
  for folding sites. store in StdIntList::rbox.
  */
  StdI->NCell = 0;
  for (ii = 0; ii < 3; ii++) {
    StdI->NCell += StdI->box[0][ii]
      * StdI->box[1][(ii + 1) % 3]
      * StdI->box[2][(ii + 2) % 3]
      - StdI->box[0][ii]
      * StdI->box[1][(ii + 2) % 3]
      * StdI->box[2][(ii + 1) % 3];
  }
  printf("         Number of Cell = %d\n", abs(StdI->NCell));
  if (StdI->NCell == 0) {
    StdFace_exit(-1);
  }

  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < 3; jj++) {
      StdI->rbox[ii][jj] = StdI->box[(ii + 1) % 3][(jj + 1) % 3] * StdI->box[(ii + 2) % 3][(jj + 2) % 3]
                         - StdI->box[(ii + 1) % 3][(jj + 2) % 3] * StdI->box[(ii + 2) % 3][(jj + 1) % 3];
    }
  }
  if (StdI->NCell < 0) {
    for (ii = 0; ii < 3; ii++)
      for (jj = 0; jj < 3; jj++)
        StdI->rbox[ii][jj] *= -1;
    StdI->NCell *= -1;
  }/*if (StdI->NCell < 0)*/
  /**@brief
  (5) Find Cells in the super-cell
  (5-1) Find the lower- and the upper bound for surrowndinf the super-cell
  */
  for (ii = 0; ii < 3; ii++) {
    bound[ii][0] = 0;
    bound[ii][1] = 0;
    for (nBox[2] = 0; nBox[2] < 2; nBox[2]++) {
      for (nBox[1] = 0; nBox[1] < 2; nBox[1]++) {
        for (nBox[0] = 0; nBox[0] < 2; nBox[0]++) {
          edge = 0;
          for (jj = 0; jj < 3; jj++) edge += nBox[jj] * StdI->box[jj][ii];
          if (edge < bound[ii][0]) bound[ii][0] = edge;
          if (edge > bound[ii][1]) bound[ii][1] = edge;
        }
      }
    }
  }
  /**@brief
  (5-2) Find cells in the super-cell (StdIntList::Cell, the fractional coordinate)
  */
  StdI->Cell = (int **)malloc(sizeof(int*) * StdI->NCell);
  for (ii = 0; ii < StdI->NCell; ii++) {
    StdI->Cell[ii] = (int *)malloc(sizeof(int) * 3);
  }/*for (ii = 0; ii < StdI->NCell; ii++)*/
  jj = 0;
  for (iCellV[2] = bound[2][0]; iCellV[2] <= bound[2][1]; iCellV[2]++) {
    for (iCellV[1] = bound[1][0]; iCellV[1] <= bound[1][1]; iCellV[1]++) {
      for (iCellV[0] = bound[0][0]; iCellV[0] <= bound[0][1]; iCellV[0]++) {
        StdFace_FoldSite(StdI, iCellV, nBox, iCellV_fold);
        if (nBox[0] == 0 && nBox[1] == 0 && nBox[2] == 0) {
          for (ii = 0; ii < 3; ii++)
            StdI->Cell[jj][ii] = iCellV[ii];
          jj += 1;
        }/*if (lUC == 1)*/
      }/*for (iCellV[0] = bound[0][0]; iCellV[0] <= bound[0][1]; iCellV[0]++*/
    }/*for (iCellV[1] = bound[1][0]; iCellV[1] <= bound[1][1]; iCellV[1]++)*/
  }/*for (iCellV[2] = bound[2][0]; iCellV[2] <= bound[2][1]; iCellV[2]++)*/
  /**@brief
  (6) For 2D system, print header of lattice.gp
  */
  if (dim == 2) {
    pos[0][0] = 0.0;
    pos[0][1] = 0.0;
    pos[1][0] = StdI->direct[0][0] * (double)StdI->box[0][0] + StdI->direct[1][0] * (double)StdI->box[0][1];
    pos[1][1] = StdI->direct[0][1] * (double)StdI->box[0][0] + StdI->direct[1][1] * (double)StdI->box[0][1];
    pos[2][0] = StdI->direct[0][0] * (double)StdI->box[1][0] + StdI->direct[1][0] * (double)StdI->box[1][1];
    pos[2][1] = StdI->direct[0][1] * (double)StdI->box[1][0] + StdI->direct[1][1] * (double)StdI->box[1][1];
    pos[3][0] = pos[1][0] + pos[2][0];
    pos[3][1] = pos[1][1] + pos[2][1];
    /**/
    xmin = 0.0;
    xmax = 0.0;
    for (ipos = 0; ipos < 4; ipos++) {
      if (pos[ipos][0] < xmin) xmin = pos[ipos][0];
      if (pos[ipos][0] > xmax) xmax = pos[ipos][0];
      if (pos[ipos][1] < xmin) xmin = pos[ipos][1];
      if (pos[ipos][1] > xmax) xmax = pos[ipos][1];
    }
    xmin -= 2.0;
    xmax += 2.0;

    fprintf(fp, "#set terminal pdf color enhanced \\\n");
    fprintf(fp, "#dashed dl 1.0 size 20.0cm, 20.0cm \n");
    fprintf(fp, "#set output \"lattice.pdf\"\n");
    fprintf(fp, "set xrange [%f: %f]\n", xmin, xmax);
    fprintf(fp, "set yrange [%f: %f]\n", xmin, xmax);
    fprintf(fp, "set size square\n");
    fprintf(fp, "unset key\n");
    fprintf(fp, "unset tics\n");
    fprintf(fp, "unset border\n");

    fprintf(fp, "set style line 1 lc 1 lt 1\n");
    fprintf(fp, "set style line 2 lc 5 lt 1\n");
    fprintf(fp, "set style line 3 lc 0 lt 1\n");

    fprintf(fp, "set arrow from %f, %f to %f, %f nohead front ls 3\n", pos[0][0], pos[0][1], pos[1][0], pos[1][1]);
    fprintf(fp, "set arrow from %f, %f to %f, %f nohead front ls 3\n", pos[1][0], pos[1][1], pos[3][0], pos[3][1]);
    fprintf(fp, "set arrow from %f, %f to %f, %f nohead front ls 3\n", pos[3][0], pos[3][1], pos[2][0], pos[2][1]);
    fprintf(fp, "set arrow from %f, %f to %f, %f nohead front ls 3\n", pos[2][0], pos[2][1], pos[0][0], pos[0][1]);
  }/*if (dim == 2)*/
}/*void StdFace_InitSite2D*/
/**
@brief Find the index of transfer and interaction
*/
void StdFace_FindSite(
  struct StdIntList *StdI,//!<[inout]
  int iW,//!<[in] position of initial site
  int iL,//!<[in] position of initial site
  int iH,//!<[in] position of initial site
  int diW,//!<[in] Translation from the initial site
  int diL,//!<[in] Translation from the initial site
  int diH,//!<[in] Translation from the initial site
  int isiteUC,//!<[in] Intrinsic site index of initial site
  int jsiteUC,//!<[in] Intrinsic site index of final site
  int *isite,//!<[out] initial site
  int *jsite,//!<[out] final site
  double complex *Cphase//!<[out] Boundary phase, if it across boundary
)
{
  int iCell, jCell, kCell, ii;
  int nBox[3], jCellV[3];
  /**/
  jCellV[0] = iW + diW;
  jCellV[1] = iL + diL;
  jCellV[2] = iH + diH;
  StdFace_FoldSite(StdI, jCellV, nBox, jCellV);
  *Cphase = 1.0;
  for (ii = 0; ii < 3; ii++) *Cphase *= cpow(StdI->ExpPhase[ii], (double)nBox[ii]);
  /**/
  for (kCell = 0; kCell < StdI->NCell; kCell++) {
    if (jCellV[0] == StdI->Cell[kCell][0] &&
        jCellV[1] == StdI->Cell[kCell][1] &&
        jCellV[2] == StdI->Cell[kCell][2]) 
    {
      jCell = kCell;
    }
    if (iW == StdI->Cell[kCell][0] &&
        iL == StdI->Cell[kCell][1] &&
        iH == StdI->Cell[kCell][2])
    {
      iCell = kCell;
    }
  }/*for (iCell = 0; iCell < StdI->NCell; iCell++)*/
  *isite = iCell * StdI->NsiteUC + isiteUC;
  *jsite = jCell * StdI->NsiteUC + jsiteUC;
  if (strcmp(StdI->model, "kondo") == 0) {
    *isite += StdI->NCell * StdI->NsiteUC;
    *jsite += StdI->NCell * StdI->NsiteUC;
  }
}/*void StdFace_FindSite*/
/**
@brief Set Label in the gnuplot display (Only used in 2D system)
*/
void StdFace_SetLabel(
  struct StdIntList *StdI,//!<[inout]
  FILE *fp,//!<[in] File pointer to lattice.gp
  int iW,//!<[in] position of initial site
  int iL,//!<[in] position of initial site
  int diW,//!<[in] Translation from the initial site
  int diL,//!<[in] Translation from the initial site
  int isiteUC,//!<[in] Intrinsic site index of initial site
  int jsiteUC,//!<[in] Intrinsic site index of final site 
  int *isite,//!<[out] initial site 
  int *jsite,//!<[out] final site 
  int connect,//!<[in] 1 for nearest neighbor, 2 for 2nd nearest
  double complex *Cphase//!<[out] Boundary phase, if it across boundary
)
{
  double xi, yi, xj, yj;
  /**@brief
  First print the reversed one
  */
  StdFace_FindSite(StdI, iW, iL, 0, -diW, -diL, 0, jsiteUC, isiteUC, isite, jsite, Cphase);

  xi = StdI->direct[0][0] * ((double)iW + StdI->tau[jsiteUC][0])
     + StdI->direct[1][0] * ((double)iL + StdI->tau[jsiteUC][1]);
  yi = StdI->direct[0][1] * ((double)iW + StdI->tau[jsiteUC][0])
     + StdI->direct[1][1] * ((double)iL + StdI->tau[jsiteUC][1]);

  xj = StdI->direct[0][0] * ((double)(iW - diW) + StdI->tau[isiteUC][0])
     + StdI->direct[1][0] * ((double)(iL - diL) + StdI->tau[isiteUC][1]);
  yj = StdI->direct[0][1] * ((double)(iW - diW) + StdI->tau[isiteUC][0])
     + StdI->direct[1][1] * ((double)(iL - diL) + StdI->tau[isiteUC][1]);

  if (*isite < 10)fprintf(fp, "set label \"%1d\" at %f, %f center front\n", *isite, xi, yi);
  else            fprintf(fp, "set label \"%2d\" at %f, %f center front\n", *isite, xi, yi);
  if (*jsite < 10)fprintf(fp, "set label \"%1d\" at %f, %f center front\n", *jsite, xj, yj);
  else            fprintf(fp, "set label \"%2d\" at %f, %f center front\n", *jsite, xj, yj);
  fprintf(fp, "set arrow from %f, %f to %f, %f nohead ls %d\n", xi, yi, xj, yj, connect);
  /**@brief
  Then print the normal one, these are different when they cross boundary.
  */
  StdFace_FindSite(StdI, iW, iL, 0, diW, diL, 0, isiteUC, jsiteUC, isite, jsite, Cphase);

  xi = StdI->direct[1][0] * ((double)iL + StdI->tau[isiteUC][1])
     + StdI->direct[0][0] * ((double)iW + StdI->tau[isiteUC][0]);
  yi = StdI->direct[1][1] * ((double)iL + StdI->tau[isiteUC][1])
     + StdI->direct[0][1] * ((double)iW + StdI->tau[isiteUC][0]);

  xj = StdI->direct[0][0] * ((double)(iW + diW) + StdI->tau[jsiteUC][0])
     + StdI->direct[1][0] * ((double)(iL + diL) + StdI->tau[jsiteUC][1]);
  yj = StdI->direct[0][1] * ((double)(iW + diW) + StdI->tau[jsiteUC][0])
     + StdI->direct[1][1] * ((double)(iL + diL) + StdI->tau[jsiteUC][1]);

  if (*isite < 10)fprintf(fp, "set label \"%1d\" at %f, %f center front\n", *isite, xi, yi);
  else            fprintf(fp, "set label \"%2d\" at %f, %f center front\n", *isite, xi, yi);
  if (*jsite < 10)fprintf(fp, "set label \"%1d\" at %f, %f center front\n", *jsite, xj, yj);
  else            fprintf(fp, "set label \"%2d\" at %f, %f center front\n", *jsite, xj, yj);
  fprintf(fp, "set arrow from %f, %f to %f, %f nohead ls %d\n", xi, yi, xj, yj, connect);
}/*void StdFace_SetLabel*/
/**
@brief Print lattice.xsf (XCrysDen format) 
*/
void StdFace_PrintXSF(struct StdIntList *StdI) {
  FILE *fp;
  int ii, jj, kk, isite, iCell;
  double vec[3];

  fp = fopen("lattice.xsf", "w");
  fprintf(fp, "CRYSTAL\n");
  fprintf(fp, "PRIMVEC\n");
  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < 3; jj++) {
      vec[jj] = 0.0;
      for (kk = 0; kk < 3; kk++)
        vec[jj] += (double)StdI->box[ii][kk] * StdI->direct[kk][jj];
    }
    fprintf(fp, "%15.5f %15.5f %15.5f\n", vec[0], vec[1], vec[2]);
  }
  fprintf(fp, "PRIMCOORD\n");
  fprintf(fp, "%d 1\n", StdI->NCell * StdI->NsiteUC);
  for (iCell = 0; iCell < StdI->NCell; iCell++) {
    for (isite = 0; isite < StdI->NsiteUC; isite++) {
      for (jj = 0; jj < 3; jj++) {
        vec[jj] = 0.0;
        for (kk = 0; kk < 3; kk++)
          vec[jj] += ((double)StdI->Cell[iCell][kk] + StdI->tau[isite][kk])
          * StdI->direct[kk][jj];
      }
      fprintf(fp, "H %15.5f %15.5f %15.5f\n", vec[0], vec[1], vec[2]);
    }
  }
  fflush(fp);
  fclose(fp);
}/*void StdFace_PrintXSF*/
/**
@brief Input nearest-neighbor spin-spin interaction
*/
void StdFace_InputSpinNN(
  struct StdIntList *StdI,//!<[inout]
  double J0[3][3],//!<[in] The anisotropic spin interaction
  double J0All,//!<[in] The isotropic interaction
  char *J0name//!<[in] The name of this spin interaction (e.g. J1)
) 
{
  int i1, i2, i3, i4;
  char Jname[3][3][10]; 
  
  strcpy(Jname[0][0], "x\0");
  strcpy(Jname[0][1], "xy\0");
  strcpy(Jname[0][2], "xz\0");
  strcpy(Jname[1][0], "yx\0");
  strcpy(Jname[1][1], "y\0");
  strcpy(Jname[1][2], "yz\0");
  strcpy(Jname[2][0], "zx\0");
  strcpy(Jname[2][1], "zy\0");
  strcpy(Jname[2][2], "z\0");

  if (isnan(StdI->JAll) == 0 && isnan(J0All)  == 0) {
    fprintf(stdout, "\n ERROR! J and %s conflict !\n\n", J0name);
    StdFace_exit(-1);
  }
  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      if (isnan(StdI->JAll) == 0 && isnan(StdI->J[i1][i2]) == 0) {
        fprintf(stdout, "\n ERROR! J and J%s conflict !\n\n", Jname[i1][i2]);
        StdFace_exit(-1);
      }
      else if (isnan(J0All) == 0 && isnan(StdI->J[i1][i2]) == 0) {
        fprintf(stdout, "\n ERROR! %s and J%s conflict !\n\n",
          J0name, Jname[i1][i2]);
        StdFace_exit(-1);
      }
      else if (isnan(J0All) == 0 && isnan(J0[i1][i2]) == 0) {
        fprintf(stdout, "\n ERROR! %s and %s%s conflict !\n\n", J0name,
          J0name, Jname[i1][i2]);
        StdFace_exit(-1);
      }
      else if (isnan(J0[i1][i2]) == 0 && isnan(StdI->JAll) == 0) {
        fprintf(stdout, "\n ERROR! %s%s and J conflict !\n\n",
          J0name, Jname[i1][i2]);
        StdFace_exit(-1);
      }
    }/*for (j = 0; j < 3; j++)*/
  }/*for (i = 0; i < 3; i++)*/
 
  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      for (i3 = 0; i3 < 3; i3++) {
        for (i4 = 0; i4 < 3; i4++) {
          if (isnan(J0[i1][i2]) == 0 && isnan(StdI->J[i3][i4]) == 0) {
            fprintf(stdout, "\n ERROR! %s%s and J%s conflict !\n\n", 
              J0name, Jname[i1][i2], Jname[i3][i4]);
            StdFace_exit(-1);
          }
        }/*for (i4 = 0; i4 < 3; i4++)*/
      }/*for (i3 = 0; i3 < 3; i3++)*/
    }/*for (j = 0; j < 3; j++)*/
  }/*for (i = 0; i < 3; i++)*/

  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      if (isnan(J0[i1][i2]) == 0)
        fprintf(stdout, "  %14s%s = %-10.5f\n", J0name, Jname[i1][i2], J0[i1][i2]);
      else if (isnan(StdI->J[i1][i2]) == 0) {
        J0[i1][i2] = StdI->J[i1][i2];
        fprintf(stdout, "  %14s%s = %-10.5f\n", J0name, Jname[i1][i2], J0[i1][i2]);
      }
      else if (i1 == i2 && isnan(J0All) == 0) {
        J0[i1][i2] = J0All;
        fprintf(stdout, "  %14s%s = %-10.5f\n", J0name, Jname[i1][i2], J0[i1][i2]);
      }
      else if (i1 == i2 && isnan(StdI->JAll) == 0) {
        J0[i1][i2] = StdI->JAll;
        fprintf(stdout, "  %14s%s = %-10.5f\n", J0name, Jname[i1][i2], J0[i1][i2]);
      }
      else {
        J0[i1][i2] = 0.0;
      }
    }/*for (i2 = 0; i2 < 3; i2++)*/
  }/*for (i = 0; i < 3; i++)*/
}/*void StdFace_InputSpinNN*/
/**
@brief Input spin-spin interaction other than nearest-neighbor
*/
void StdFace_InputSpin(
  struct StdIntList *StdI,//!<[inout] 
  double Jp[3][3],//!<[in] Fully anisotropic spin interaction
  double JpAll,//!<[in] The isotropic interaction
  char *Jpname//!<The name of this spin interaction(e.g.J')
)
{
  int i1, i2;
  char Jname[3][3][10];

  strcpy(Jname[0][0], "x\0");
  strcpy(Jname[0][1], "xy\0");
  strcpy(Jname[0][2], "xz\0");
  strcpy(Jname[1][0], "yx\0");
  strcpy(Jname[1][1], "y\0");
  strcpy(Jname[1][2], "yz\0");
  strcpy(Jname[2][0], "zx\0");
  strcpy(Jname[2][1], "zy\0");
  strcpy(Jname[2][2], "z\0");

  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      if (isnan(JpAll) == 0 && isnan(Jp[i1][i2]) == 0) {
        fprintf(stdout, "\n ERROR! %s and %s%s conflict !\n\n", Jpname,
          Jpname, Jname[i1][i2]);
        StdFace_exit(-1);
      }
    }/*for (j = 0; j < 3; j++)*/
  }/*for (i = 0; i < 3; i++)*/

  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      if (isnan(Jp[i1][i2]) == 0)
        fprintf(stdout, "  %14s%s = %-10.5f\n", Jpname, Jname[i1][i2], Jp[i1][i2]);
      else if (i1 == i2 && isnan(JpAll) == 0) {
        Jp[i1][i2] = JpAll;
        fprintf(stdout, "  %14s%s = %-10.5f\n", Jpname, Jname[i1][i2], Jp[i1][i2]);
      }
      else {
        Jp[i1][i2] = 0.0;
      }
    }/*for (i2 = 0; i2 < 3; i2++)*/
  }/*for (i = 0; i < 3; i++)*/
}/*void StdFace_InputSpin*/
/**
@brief Input off-site Coulomb interaction from the 
input file, if it is not specified, use the default value (0
or the isotropic Coulomb interaction StdIntList::V).
*/
void StdFace_InputCoulombV(
  struct StdIntList *StdI,//!<[inout]
  double *V0,//<!<[in]
  char *V0name//!<[in] E.g. V1
)
{
  if (isnan(StdI->V) == 0 && isnan(*V0) == 0) {
    fprintf(stdout, "\n ERROR! V and %s conflict !\n\n", V0name);
    StdFace_exit(-1);
  }
  else if (isnan(*V0) == 0)
    fprintf(stdout, "  %15s = %-10.5f\n", V0name, *V0);
  else if (isnan(StdI->V) == 0) {
    *V0 = StdI->V;
    fprintf(stdout, "  %15s = %-10.5f\n", V0name, *V0);
  }
  else {
    *V0 = 0.0;
  }
}/*void StdFace_InputCoulombV*/
/**
@brief Input hopping integral from the
input file, if it is not specified, use the default value(0
or the isotropic hopping StdIntList::V).
*/
void StdFace_InputHopp(
  struct StdIntList *StdI,//!<[inout]
  double complex *t0,//!<[in]
  char *t0name//!<[in] E.g. t1
)
{
  if (isnan(creal(StdI->t)) == 0 && isnan(creal(*t0)) == 0) {
    fprintf(stdout, "\n ERROR! t and %s conflict !\n\n", t0name);
    StdFace_exit(-1);
  }
  else if (isnan(creal(*t0)) == 0)
    fprintf(stdout, "  %15s = %-10.5f %-10.5f\n", t0name, creal(*t0), cimag(*t0));
  else if (isnan(creal(StdI->t)) == 0) {
    *t0 = StdI->t;
    fprintf(stdout, "  %15s = %-10.5f %-10.5f\n", t0name, creal(*t0), cimag(*t0));
  }
  else {
    *t0 = 0.0;
  }
}/*void StdFace_InputHopp*/
/**
@brief Print geometry of sites for the pos-process of correlation function
*/
void StdFace_PrintGeometry(struct StdIntList *StdI) {
  FILE *fp;
  int isite, iCell;

  fp = fopen("geometry.dat", "w");

  fprintf(fp, "%25.15e %25.15e %25.15e\n", StdI->direct[0][0], StdI->direct[0][1], StdI->direct[0][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n", StdI->direct[1][0], StdI->direct[1][1], StdI->direct[1][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n", StdI->direct[2][0], StdI->direct[2][1], StdI->direct[2][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n", StdI->phase[0], StdI->phase[1], StdI->phase[2]);
  fprintf(fp, "%d %d %d\n", StdI->box[0][0], StdI->box[0][1], StdI->box[0][2]);
  fprintf(fp, "%d %d %d\n", StdI->box[1][0], StdI->box[1][1], StdI->box[1][2]);
  fprintf(fp, "%d %d %d\n", StdI->box[2][0], StdI->box[2][1], StdI->box[2][2]);

  for (iCell = 0; iCell < StdI->NCell; iCell++) {
    for (isite = 0; isite < StdI->NsiteUC; isite++) {
      fprintf(fp, "%25.15e %25.15e %25.15e\n",
        StdI->tau[isite][0] + (double)StdI->Cell[iCell][0],
        StdI->tau[isite][1] + (double)StdI->Cell[iCell][1],
        StdI->tau[isite][2] + (double)StdI->Cell[iCell][2]);
    }/*for (isite = 0; isite < StdI->NsiteUC; isite++)*/
  }/* for (iCell = 0; iCell < StdI->NCell; iCell++)*/
  if (strcmp(StdI->model, "kondo") == 0) {
    for (iCell = 0; iCell < StdI->NCell; iCell++) {
      for (isite = 0; isite < StdI->NsiteUC; isite++) {
        fprintf(fp, "%25.15e %25.15e %25.15e\n",
          StdI->tau[isite][0] + (double)StdI->Cell[iCell][0],
          StdI->tau[isite][1] + (double)StdI->Cell[iCell][1],
          StdI->tau[isite][2] + (double)StdI->Cell[iCell][2]);
      }/*for (isite = 0; isite < StdI->NsiteUC; isite++)*/
    }/* for (iCell = 0; iCell < StdI->NCell; iCell++)*/
  }
  fflush(fp);
  fclose(fp);
}/*void StdFace_PrintGeometry()*/
/**
@brief Malloc Arrays for interactions
*/
void StdFace_MallocInteractions(
  struct StdIntList *StdI,//!<[inout]
  int ntransMax,//!<[in] upper limit of the number of transfer
  int nintrMax//!<[in] upper limit of the number of interaction
) {
  int ii;
  /**@brief
  (1) Transfer StdIntList::trans, StdIntList::transindx
  */
  StdI->transindx = (int **)malloc(sizeof(int*) * ntransMax);
  StdI->trans = (double complex *)malloc(sizeof(double complex) * ntransMax);
  for (ii = 0; ii < ntransMax; ii++) {
    StdI->transindx[ii] = (int *)malloc(sizeof(int) * 4);
  }
  StdI->ntrans = 0;
  /**@brief
  (2) InterAll StdIntList::intr, StdIntList::intrindx
  */
  StdI->intrindx = (int **)malloc(sizeof(int*) * nintrMax);
  StdI->intr = (double complex *)malloc(sizeof(double complex) * nintrMax);
  for (ii = 0; ii < nintrMax; ii++) {
    StdI->intrindx[ii] = (int *)malloc(sizeof(int) * 8);
  }
  StdI->nintr = 0;
  /**@brief
  (3) Coulomb intra StdIntList::Cintra, StdIntList::CintraIndx
  */
  StdI->CintraIndx = (int **)malloc(sizeof(int*) * nintrMax);
  StdI->Cintra = (double *)malloc(sizeof(double) * nintrMax);
  for (ii = 0; ii < nintrMax; ii++) {
    StdI->CintraIndx[ii] = (int *)malloc(sizeof(int) * 1);
  }
  StdI->NCintra = 0;
  /**@brief
  (4) Coulomb inter StdIntList::Cinter, StdIntList::CinterIndx
  */
  StdI->CinterIndx = (int **)malloc(sizeof(int*) * nintrMax);
  StdI->Cinter = (double *)malloc(sizeof(double) * nintrMax);
  for (ii = 0; ii < nintrMax; ii++) {
    StdI->CinterIndx[ii] = (int *)malloc(sizeof(int) * 2);
  }
  StdI->NCinter = 0;
  /**@brief
  (5) Hund StdIntList::Hund, StdIntList::HundIndx
  */
  StdI->HundIndx = (int **)malloc(sizeof(int*) * nintrMax);
  StdI->Hund = (double *)malloc(sizeof(double) * nintrMax);
  for (ii = 0; ii < nintrMax; ii++) {
    StdI->HundIndx[ii] = (int *)malloc(sizeof(int) * 2);
  }
  StdI->NHund = 0;
  /**@brief
  (6) Excahnge StdIntList::Ex, StdIntList::ExIndx
  */
  StdI->ExIndx = (int **)malloc(sizeof(int*) * nintrMax);
  StdI->Ex = (double *)malloc(sizeof(double) * nintrMax);
  for (ii = 0; ii < nintrMax; ii++) {
    StdI->ExIndx[ii] = (int *)malloc(sizeof(int) * 2);
  }
  StdI->NEx = 0;
  /**@brief
  (7) PairLift StdIntList::PairLift, StdIntList::PLIndx
  */
  StdI->PLIndx = (int **)malloc(sizeof(int*) * nintrMax);
  StdI->PairLift = (double *)malloc(sizeof(double) * nintrMax);
  for (ii = 0; ii < nintrMax; ii++) {
    StdI->PLIndx[ii] = (int *)malloc(sizeof(int) * 2);
  }
  StdI->NPairLift = 0;
}/*void StdFace_MallocInteractions*/
#if defined(_mVMC)
/**
@brief Define whether the specified site is in the unit cell or not.
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StdFace_FoldSiteSub(
  struct StdIntList *StdI,//!<[inout]
  int iCellV[3],//!<[in]
  int nBox[3], //!<[out]
  int iCellV_fold[3]//!<[out]
)
{
  int ii, jj, iCellV_frac[3];
  /**@brief
  (1) Transform to fractional coordinate (times NCell)
  */
  for (ii = 0; ii < 3; ii++) {
    iCellV_frac[ii] = 0;
    for (jj = 0; jj < 3; jj++)iCellV_frac[ii] += StdI->rboxsub[ii][jj] * iCellV[jj];
  }
  /**@brief
  (2) Which supercell contains this cell
  */
  for (ii = 0; ii < 3; ii++)
    nBox[ii] = (iCellV_frac[ii] + StdI->NCellsub * 1000) / StdI->NCellsub - 1000;
  /**@brief
  (3) Fractional coordinate (times NCell) in the original supercell
  */
  for (ii = 0; ii < 3; ii++)
    iCellV_frac[ii] -= StdI->NCellsub*(nBox[ii]);
  /**/
  for (ii = 0; ii < 3; ii++) {
    iCellV_fold[ii] = 0;
    for (jj = 0; jj < 3; jj++) iCellV_fold[ii] += StdI->boxsub[jj][ii] * iCellV_frac[jj];
    iCellV_fold[ii] = (iCellV_fold[ii] + StdI->NCellsub * 1000) / StdI->NCellsub - 1000;
  }
}/*static void StdFace_FoldSiteSub*/
/**
@brief Print Quantum number projection
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Proj(struct StdIntList *StdI)
{
  FILE *fp;
  int jsite, iCell, jCell, kCell;
  int nBox[3], iCellV[3], jCellV[3], ii;
  int iSym;
  int **Sym, **Anti;

  Sym = (int **)malloc(sizeof(int*) * StdI->nsite);
  Anti = (int **)malloc(sizeof(int*) * StdI->nsite);
  for (jsite = 0; jsite < StdI->nsite; jsite++) {
    Sym[jsite] = (int *)malloc(sizeof(int) * StdI->nsite);
    Anti[jsite] = (int *)malloc(sizeof(int) * StdI->nsite);
  }
  /**pbrief
  (1) Define translation operator in sub lattice
  */
  StdI->NSym = 0;
  for (iCell = 0; iCell < StdI->NCell; iCell++) {

    StdFace_FoldSiteSub(StdI, StdI->Cell[iCell], nBox, iCellV);

    StdFace_FoldSite(StdI, iCellV, nBox, iCellV);

    if (iCellV[0] == StdI->Cell[iCell][0] && 
        iCellV[1] == StdI->Cell[iCell][1] && 
        iCellV[2] == StdI->Cell[iCell][2]) {
      /**@brief
      (2) Translation operator in sub lattice
      */
      for (jCell = 0; jCell < StdI->NCell; jCell++) {

        for (ii = 0; ii < 3; ii++)jCellV[ii] = StdI->Cell[jCell][ii] + iCellV[ii];
        StdFace_FoldSite(StdI, jCellV, nBox, jCellV);

        for (kCell = 0; kCell < StdI->NCell; kCell++) {
          if (jCellV[0] == StdI->Cell[kCell][0] && 
              jCellV[1] == StdI->Cell[kCell][1] &&
              jCellV[2] == StdI->Cell[kCell][2]) 
          {

            for (jsite = 0; jsite < StdI->NsiteUC; jsite++) {

              Sym[StdI->NSym][jCell*StdI->NsiteUC + jsite] = kCell*StdI->NsiteUC + jsite;
              Anti[StdI->NSym][jCell*StdI->NsiteUC + jsite]
                = StdI->AntiPeriod[0] * nBox[0]
                + StdI->AntiPeriod[1] * nBox[1]
                + StdI->AntiPeriod[2] * nBox[2];

              if (strcmp(StdI->model, "kondo") == 0) {
                Sym[StdI->NSym][StdI->nsite / 2 + jCell*StdI->NsiteUC + jsite] = StdI->nsite / 2 + kCell*StdI->NsiteUC + jsite;
                Anti[StdI->NSym][StdI->nsite / 2 + jCell*StdI->NsiteUC + jsite]
                  = StdI->AntiPeriod[0] * nBox[0]
                  + StdI->AntiPeriod[1] * nBox[1]
                  + StdI->AntiPeriod[2] * nBox[2];
              }/*if (strcmp(StdI->model, "kondo") == 0)*/

            }/*for (jsite = 0; jsite < StdI->NsiteUC; jsite++)*/

          }/*if (jWfold == StdI->Cell[kCell][0] && jLfold == StdI->Cell[kCell][1])*/
        }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/
      }/*for (jCell = 0; jCell < StdI->NCell; jCell++)*/
      StdI->NSym += 1;
    }/*if (iWfold == iW && iLfold == iL)*/
  }/*for (iCell = 0; iCell < StdI->NCell; iCell++)*/

  fp = fopen("qptransidx.def", "w");
  fprintf(fp, "=============================================\n");
  fprintf(fp, "NQPTrans %10d\n", StdI->NSym);
  fprintf(fp, "=============================================\n");
  fprintf(fp, "======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n");
  fprintf(fp, "=============================================\n");
  for (iSym = 0; iSym < StdI->NSym; iSym++) {
    fprintf(fp, "%d %10.5f\n", iSym, 1.0);
  }
  for (iSym = 0; iSym < StdI->NSym; iSym++) {
    for (jsite = 0; jsite < StdI->nsite; jsite++) {
      if (Anti[iSym][jsite] % 2 == 0) Anti[iSym][jsite] = 1;
      else Anti[iSym][jsite] = -1;
      if (StdI->AntiPeriod[0] == 1 || StdI->AntiPeriod[1] == 1 || StdI->AntiPeriod[2] == 1) {
        fprintf(fp, "%5d  %5d  %5d  %5d\n", iSym, jsite, Sym[iSym][jsite], Anti[iSym][jsite]);
      }
      else {
        fprintf(fp, "%5d  %5d  %5d\n", iSym, jsite, Sym[iSym][jsite]);
      }
    }
  }
  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    qptransidx.def is written.\n");

  for (jsite = 0; jsite < StdI->nsite; jsite++) {
    free(Anti[jsite]);
    free(Sym[jsite]);
  }
  free(Sym);
  free(Anti);
}/*void StdFace_Proj(struct StdIntList *StdI)*/
/**
@brief Initialize sub Cell
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StdFace_InitSiteSub(struct StdIntList *StdI)
{
  int ii, jj, kk, prod;
  /**@brief
  (1) check Input parameters
  */
  if ((StdI->Lsub != StdI->NaN_i || StdI->Wsub != StdI->NaN_i || StdI->Hsub != StdI->NaN_i)
    && (StdI->boxsub[0][0] != StdI->NaN_i || StdI->boxsub[0][1] != StdI->NaN_i || StdI->boxsub[0][2] != StdI->NaN_i ||
        StdI->boxsub[1][0] != StdI->NaN_i || StdI->boxsub[1][1] != StdI->NaN_i || StdI->boxsub[1][2] != StdI->NaN_i ||
        StdI->boxsub[2][0] != StdI->NaN_i || StdI->boxsub[2][1] != StdI->NaN_i || StdI->boxsub[2][2] != StdI->NaN_i))
  {
    fprintf(stdout, "\nERROR ! (Lsub, Wsub, Hsub) and (a0Wsub, ..., a2Hsub) conflict !\n\n");
    StdFace_exit(-1);
  }
  else if (StdI->Wsub != StdI->NaN_i || StdI->Lsub != StdI->NaN_i || StdI->Hsub != StdI->NaN_i) {
    StdFace_PrintVal_i("Lsub", &StdI->Lsub, 1);
    StdFace_PrintVal_i("Wsub", &StdI->Wsub, 1);
    StdFace_PrintVal_i("Hsub", &StdI->Hsub, 1);
    for (ii = 0; ii < 3; ii++) for (jj = 0; jj < 3; jj++)
      StdI->boxsub[ii][jj] = 0;
    StdI->boxsub[0][0] = StdI->Wsub;
    StdI->boxsub[1][1] = StdI->Lsub;
    StdI->boxsub[2][2] = StdI->Hsub;
  }
  else {
    StdFace_PrintVal_i("a0Wsub", &StdI->boxsub[0][0], StdI->box[0][0]);
    StdFace_PrintVal_i("a0Lsub", &StdI->boxsub[0][1], StdI->box[0][1]);
    StdFace_PrintVal_i("a0Hsub", &StdI->boxsub[0][2], StdI->box[0][2]);
    StdFace_PrintVal_i("a1Wsub", &StdI->boxsub[1][0], StdI->box[1][0]);
    StdFace_PrintVal_i("a1Lsub", &StdI->boxsub[1][1], StdI->box[1][1]);
    StdFace_PrintVal_i("a1Hsub", &StdI->boxsub[1][2], StdI->box[1][2]);
    StdFace_PrintVal_i("a2Wsub", &StdI->boxsub[2][0], StdI->box[2][0]);
    StdFace_PrintVal_i("a2Lsub", &StdI->boxsub[2][1], StdI->box[2][1]);
    StdFace_PrintVal_i("a2Hsub", &StdI->boxsub[2][2], StdI->box[2][2]);
  }
  /**@brief
  (2) Calculate reciprocal lattice vectors (times NCellsub)
  */
  StdI->NCellsub = 0;
  for (ii = 0; ii < 3; ii++) {
    StdI->NCellsub += StdI->boxsub[0][ii]
      * StdI->boxsub[1][(ii + 1) % 3]
      * StdI->boxsub[2][(ii + 2) % 3]
      - StdI->boxsub[0][ii]
      * StdI->boxsub[1][(ii + 2) % 3]
      * StdI->boxsub[2][(ii + 1) % 3];
  }
  printf("         Number of Cell in the sublattice: %d\n", abs(StdI->NCellsub));
  if (StdI->NCellsub == 0) {
    StdFace_exit(-1);
  }

  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < 3; jj++) {
      StdI->rboxsub[ii][jj] = StdI->boxsub[(ii + 1) % 3][(jj + 1) % 3] * StdI->boxsub[(ii + 2) % 3][(jj + 2) % 3]
        - StdI->boxsub[(ii + 1) % 3][(jj + 2) % 3] * StdI->boxsub[(ii + 2) % 3][(jj + 1) % 3];
    }
  }
  if (StdI->NCellsub < 0) {
    for (ii = 0; ii < 3; ii++)
      for (jj = 0; jj < 3; jj++)
        StdI->rboxsub[ii][jj] *= -1;
    StdI->NCellsub *= -1;
  }/*if (StdI->NCell < 0)*/
  /**@brief
  (4) Check the sublattice is commensurate or not
  */
  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < 3; jj++) {
      prod = 0.0;
      for (kk = 0; kk < 3; kk++) prod += StdI->rboxsub[ii][kk] * (double)StdI->box[jj][kk];
      if (prod % StdI->NCellsub != 0) {
        printf("\n ERROR ! Sublattice is INCOMMENSURATE !\n\n");
        StdFace_exit(-1);
      }/*if (prod % StdI->NCellsub != 0)*/
    }
  }
}/*void StdFace_InitSiteSub*/
/**
@brief Generate orbitalindex
*/
void StdFace_generate_orb(struct StdIntList *StdI) {
  int iCell, jCell, kCell, iCell2, jCell2, iOrb, isite, jsite, Anti;
  int nBox[3], iCellV[3], jCellV[3], dCellV[3], ii;
  int **CellDone;

  StdFace_InitSiteSub(StdI);

  StdI->Orb = (int **)malloc(sizeof(int*) * StdI->nsite);
  StdI->AntiOrb = (int **)malloc(sizeof(int*) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++) {
    StdI->Orb[isite] = (int *)malloc(sizeof(int) * StdI->nsite);
    StdI->AntiOrb[isite] = (int *)malloc(sizeof(int) * StdI->nsite);
  }
  CellDone = (int **)malloc(sizeof(int*) * StdI->NCell);
  for (iCell = 0; iCell < StdI->NCell; iCell++) {
    CellDone[iCell] = (int *)malloc(sizeof(int) * StdI->NCell);
    for (jCell = 0; jCell < StdI->NCell; jCell++) {
      CellDone[iCell][jCell] = 0;
    }
  }

  iOrb = 0;
  for (iCell = 0; iCell < StdI->NCell; iCell++) {

    StdFace_FoldSiteSub(StdI, StdI->Cell[iCell], nBox, iCellV);

    StdFace_FoldSite(StdI, iCellV, nBox, iCellV);

    for (kCell = 0; kCell < StdI->NCell; kCell++) {
      if (iCellV[0] == StdI->Cell[kCell][0] && 
          iCellV[1] == StdI->Cell[kCell][1] &&
          iCellV[2] == StdI->Cell[kCell][2]) 
      {
        iCell2 = kCell;
      }
    }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/

    for (jCell = 0; jCell < StdI->NCell; jCell++) {

      for (ii = 0; ii < 3; ii++)
        jCellV[ii] = StdI->Cell[jCell][ii] + iCellV[ii] - StdI->Cell[iCell][ii];

      StdFace_FoldSite(StdI, jCellV, nBox, jCellV);

      for (kCell = 0; kCell < StdI->NCell; kCell++) {
        if (jCellV[0] == StdI->Cell[kCell][0] &&
            jCellV[1] == StdI->Cell[kCell][1] &&
            jCellV[2] == StdI->Cell[kCell][2]) 
        {
          jCell2 = kCell;
        }
      }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/
      /*
       AntiPeriodic factor
      */
      for (ii = 0; ii < 3; ii++)
        dCellV[ii] = StdI->Cell[jCell][ii] - StdI->Cell[iCell][ii];
      StdFace_FoldSite(StdI, dCellV, nBox, dCellV);
      Anti = 0;
      for (ii = 0; ii < 3; ii++)Anti += StdI->AntiPeriod[ii] * nBox[ii];
      if (Anti % 2 == 0) Anti = 1;
      else Anti = -1;

      for (isite = 0; isite < StdI->NsiteUC; isite++) {
        for (jsite = 0; jsite < StdI->NsiteUC; jsite++) {
 
          if (CellDone[iCell2][jCell2] == 0) {
            StdI->Orb[iCell2*StdI->NsiteUC + isite][jCell2*StdI->NsiteUC + jsite] = iOrb;
            StdI->AntiOrb[iCell2*StdI->NsiteUC + isite][jCell2*StdI->NsiteUC + jsite] = Anti;
            iOrb += 1;
          }
          StdI->Orb[iCell*StdI->NsiteUC + isite][jCell*StdI->NsiteUC + jsite]
            = StdI->Orb[iCell2*StdI->NsiteUC + isite][jCell2*StdI->NsiteUC + jsite];
          StdI->AntiOrb[iCell*StdI->NsiteUC + isite][jCell*StdI->NsiteUC + jsite] = Anti;

          if (strcmp(StdI->model, "kondo") == 0) {
            if (CellDone[iCell2][jCell2] == 0) {
              StdI->Orb[StdI->nsite / 2 + iCell2*StdI->NsiteUC + isite]
                       [                  jCell2*StdI->NsiteUC + jsite] = iOrb;
              StdI->AntiOrb[StdI->nsite / 2 + iCell2*StdI->NsiteUC + isite]
                           [                  jCell2*StdI->NsiteUC + jsite] = Anti;
              iOrb += 1;
              StdI->Orb[                  iCell2*StdI->NsiteUC + isite]
                       [StdI->nsite / 2 + jCell2*StdI->NsiteUC + jsite] = iOrb;
              StdI->AntiOrb[                  iCell2*StdI->NsiteUC + isite]
                           [StdI->nsite / 2 + jCell2*StdI->NsiteUC + jsite] = Anti;
              iOrb += 1;
              StdI->Orb[StdI->nsite / 2 + iCell2*StdI->NsiteUC + isite]
                       [StdI->nsite / 2 + jCell2*StdI->NsiteUC + jsite] = iOrb;
              StdI->AntiOrb[StdI->nsite / 2 + iCell2*StdI->NsiteUC + isite]
                           [StdI->nsite / 2 + jCell2*StdI->NsiteUC + jsite] = Anti;
              iOrb += 1;
            }
            StdI->Orb[StdI->nsite / 2 + iCell*StdI->NsiteUC + isite]
                     [                  jCell*StdI->NsiteUC + jsite]
            = StdI->Orb[StdI->nsite / 2 + iCell2*StdI->NsiteUC + isite]
                       [                  jCell2*StdI->NsiteUC + jsite];
            StdI->AntiOrb[StdI->nsite / 2 + iCell*StdI->NsiteUC + isite]
                         [                  jCell*StdI->NsiteUC + jsite] = Anti;
            StdI->Orb[                  iCell*StdI->NsiteUC + isite]
                     [StdI->nsite / 2 + jCell*StdI->NsiteUC + jsite]
            = StdI->Orb[                  iCell2*StdI->NsiteUC + isite]
                       [StdI->nsite / 2 + jCell2*StdI->NsiteUC + jsite];
            StdI->AntiOrb[iCell*StdI->NsiteUC + isite]
                         [StdI->nsite / 2 + jCell*StdI->NsiteUC + jsite] = Anti;
            StdI->Orb[StdI->nsite / 2 + iCell*StdI->NsiteUC + isite]
                     [StdI->nsite / 2 + jCell*StdI->NsiteUC + jsite]
              = StdI->Orb[StdI->nsite / 2 + iCell2*StdI->NsiteUC + isite]
                         [StdI->nsite / 2 + jCell2*StdI->NsiteUC + jsite];
              StdI->AntiOrb[StdI->nsite / 2 + iCell*StdI->NsiteUC + isite]
                           [StdI->nsite / 2 + jCell*StdI->NsiteUC + jsite] = Anti;
          }/*if (strcmp(StdI->model, "kondo") == 0)*/

        }/*for (jsite = 0; jsite < StdI->NsiteUC; jsite++)*/
      }/*for (isite = 0; isite < StdI->NsiteUC; isite++)*/
      CellDone[iCell2][jCell2] = 1;
    }/*for (jCell = 0; jCell < StdI->NCell; jCell++)*/

  }/*for (iCell = 0; iCell < StdI->NCell; iCell++)*/
  StdI->NOrb = iOrb;

  for (iCell = 0; iCell < StdI->NCell; iCell++) free(CellDone[iCell]);
  free(CellDone);
}/*void StdFace_generate_orb*/
/**
@brief Output Jastrow
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void PrintJastrow(struct StdIntList *StdI) {
  FILE *fp;
  int isite, jsite, isiteUC, jsiteUC, revarsal, isite1, jsite1, iorb;
  int NJastrow, iJastrow;
  int dCell, iCell;//, jCell, dCellv[3];
  int **Jastrow;
  double complex Cphase;

  Jastrow = (int **)malloc(sizeof(int*) * StdI->nsite);
  for (isite = 0; isite < StdI->nsite; isite++) 
    Jastrow[isite] = (int *)malloc(sizeof(int) * StdI->nsite);

  if (abs(StdI->NMPTrans) == 1 || StdI->NMPTrans == StdI->NaN_i) {
    /**@brief
    (1) Copy Orbital index
    */
    for (isite = 0; isite < StdI->nsite; isite++) {
      for (jsite = 0; jsite < StdI->nsite; jsite++) {
        Jastrow[isite][jsite] = StdI->Orb[isite][jsite];
      }/*for (jsite = 0; jsite < isite; jsite++)*/
    }/*for (isite = 0; isite < StdI->nsite; isite++)*/
    /**@brief
    (2) Symmetrize
    */
    for (iorb = 0; iorb < StdI->NOrb; iorb++) {
      for (isite = 0; isite < StdI->nsite; isite++) {
        for (jsite = 0; jsite < StdI->nsite; jsite++) {
          if (Jastrow[isite][jsite] == iorb) {
            Jastrow[jsite][isite] = Jastrow[isite][jsite];
          }
        }/*for (jsite = 0; jsite < isite; jsite++)*/
      }/*for (isite = 0; isite < StdI->nsite; isite++)*/
    }/*for (iorb = 0; iorb < StdI->NOrb; iorb++)*/
     /**/
    if (strcmp(StdI->model, "hubbard") == 0) NJastrow = 0;
    else NJastrow = -1;
    for (isite = 0; isite < StdI->nsite; isite++) {
      /*
      For Local spin
      */
      if (StdI->locspinflag[isite] != 0) {
        for (jsite = 0; jsite < StdI->nsite; jsite++) {
          Jastrow[isite][jsite] = -1;
          Jastrow[jsite][isite] = -1;
        }
        continue;
      }
      /**/
      for (jsite = 0; jsite < isite; jsite++) {
        if (Jastrow[isite][jsite] >= 0) {
          iJastrow = Jastrow[isite][jsite];
          NJastrow -= 1;
          for (isite1 = 0; isite1 < StdI->nsite; isite1++) {
            for (jsite1 = 0; jsite1 < StdI->nsite; jsite1++) {
              if (Jastrow[isite1][jsite1] == iJastrow)
                Jastrow[isite1][jsite1] = NJastrow;
            }/*for (jsite1 = 0; jsite1 < StdI->nsite; jsite1++)*/
          }/*for (isite1 = 0; isite1 < StdI->nsite; isite1++)*/
        }/*if (Jastrow[isite][jsite] >= 0)*/
      }/*for (jsite = 0; jsite < isite; jsite++)*/
    }/*for (isite = 0; isite < StdI->nsite; isite++)*/
     /**/
    NJastrow = -NJastrow;
    for (isite = 0; isite < StdI->nsite; isite++) {
      for (jsite = 0; jsite < StdI->nsite; jsite++) {
        Jastrow[isite][jsite] = -1 - Jastrow[isite][jsite];
      }/*for (jsite = 0; jsite < isite; jsite++)*/
    }/*for (isite = 0; isite < StdI->nsite; isite++)*/
  }/*if (abs(StdI->NMPTrans) == 1)*/
  else {

    if (strcmp(StdI->model, "spin") == 0) {
      NJastrow = 1;

      for (isite = 0; isite < StdI->nsite; isite++) {
        for (jsite = 0; jsite < StdI->nsite; jsite++) {
          Jastrow[isite][jsite] = 0;
        }/*for (jsite = 0; jsite < StdI->nsite; jsite++)*/
      }/*for (isite = 0; isite < StdI->nsite; isite++)*/
    }/*if (strcmp(StdI->model, "spin") == 0)*/
    else {

      NJastrow = 0;

      if (strcmp(StdI->model, "kondo") == 0) {
        /*
        Local spin - itererant electron part
        */
        for (isite = 0; isite < StdI->nsite; isite++) {
          for (jsite = 0; jsite < StdI->nsite / 2; jsite++) {
            Jastrow[isite][jsite] = 0;
            Jastrow[jsite][isite] = 0;
          }/*for (jsite = 0; jsite < StdI->nsite; jsite++)*/
        }/*for (isite = 0; isite < StdI->nsite; isite++)*/

        NJastrow += 1;
      }/*if (strcmp(StdI->model, "kondo") == 0)*/

      for (dCell = 0; dCell < StdI->NCell; dCell++) {
        StdFace_FindSite(StdI,
          0, 0, 0,
          -StdI->Cell[dCell][0], -StdI->Cell[dCell][1], -StdI->Cell[dCell][2],
          0, 0, &isite, &jsite, &Cphase);
        if (strcmp(StdI->model, "kondo") == 0) jsite += -StdI->NCell * StdI->NsiteUC;
        iCell = jsite / StdI->NsiteUC;
        if (iCell < dCell) {
          /*
          If -R has been already done, skip.
          */
          continue;
        }
        else if (iCell == dCell) {
          /*
          If revarsal symmetry [Fold(-R) = R], J(R,i,j) = J(R,j,i)
          */
          revarsal = 1;
        }
        else revarsal = 0;

        for (isiteUC = 0; isiteUC < StdI->NsiteUC; isiteUC++) {
          for (jsiteUC = 0; jsiteUC < StdI->NsiteUC; jsiteUC++) {
            if (revarsal == 1 && jsiteUC > isiteUC) continue;/*If [Fold(-R) = R]*/
            if (isiteUC == jsiteUC &&
              StdI->Cell[dCell][0] == 0 &&
              StdI->Cell[dCell][1] == 0 &&
              StdI->Cell[dCell][2] == 0) continue;/*Diagonal*/

            for (iCell = 0; iCell < StdI->NCell; iCell++) {
              StdFace_FindSite(StdI,
                StdI->Cell[iCell][0], StdI->Cell[iCell][1], StdI->Cell[iCell][2],
                StdI->Cell[dCell][0], StdI->Cell[dCell][1], StdI->Cell[dCell][2],
                isiteUC, jsiteUC, &isite, &jsite, &Cphase);

              Jastrow[isite][jsite] = NJastrow;
              Jastrow[jsite][isite] = NJastrow;

            }/*for (iCell = 0; iCell < StdI->NCell; iCell++)*/

            NJastrow += 1;

          }/*for (jsiteUC = 0; jsiteUC < StdI->NsiteUC; jsiteUC++)*/
        }/*for (isiteUC = 0; isiteUC < StdI->NsiteUC; isiteUC++)*/
      }/*for (dCell = 0; dCell < StdI->NCell; dCell++)*/
    }/*if (strcmp(StdI->model, "spin") != 0)*/
  }/*if (abs(StdI->NMPTrans) != 1)*/
    
  fp = fopen("jastrowidx.def", "w");
  fprintf(fp, "=============================================\n");
  fprintf(fp, "NJastrowIdx %10d\n", NJastrow);
  fprintf(fp, "ComplexType %10d\n", 0);
  fprintf(fp, "=============================================\n");
  fprintf(fp, "=============================================\n");

  for (isite = 0; isite < StdI->nsite; isite++) {
    for (jsite = 0; jsite < StdI->nsite; jsite++) {
      if (isite == jsite) continue;
      fprintf(fp, "%5d  %5d  %5d\n", isite, jsite, Jastrow[isite][jsite]);
    }/*for (jsite = 0; jsite < isite; jsite++)*/
  }/*for (isite = 0; isite < StdI->nsite; isite++)*/

  for (iJastrow = 0; iJastrow < NJastrow; iJastrow++){
    if (strcmp(StdI->model, "hubbard") == 0 || iJastrow > 0)
      fprintf(fp, "%5d  %5d\n", iJastrow, 1);
    else
      fprintf(fp, "%5d  %5d\n", iJastrow, 0);
  }/*for (iJastrow = 0; iJastrow < NJastrow; iJastrow++)*/
  fflush(fp);
  fclose(fp);
  fprintf(stdout, "    jastrowidx.def is written.\n");

  for (isite = 0; isite < StdI->nsite; isite++) free(Jastrow[isite]);
  free(Jastrow);
}/*void PrintJastrow*/
#endif
