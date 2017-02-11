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
*
* MPI Abortation wrapper
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_exit(int errorcode /**< [in]*/)
{
  int ierr;
  fflush(stdout);
  fflush(stderr);
#ifdef MPI
  fprintf(stdout, "\n\n #######  You DO NOT have to WORRY about the following MPI-ERROR MESSAGE.  #######\n\n");
  ierr = MPI_Abort(MPI_COMM_WORLD, errorcode);
  ierr = MPI_Finalize();
  if (ierr != 0) fprintf(stderr, "\n  MPI_Finalize() = %d\n\n", ierr);
#endif
  exit(errorcode);
}

/**
 *
 * Add transfer to the list
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_trans(
struct StdIntList *StdI,
  double complex trans0 /**< [in] Hopping integral t, mu, etc. */,
  int isite /**< [in] i for c_{i sigma}^dagger*/, 
  int ispin /**< [in] sigma for c_{i sigma}^dagger*/,
  int jsite /**< [in] j for c_{j sigma'}*/,
  int jspin /**< [in] sigma' for c_{j sigma'}*/)
{
  StdI->trans[StdI->ntrans] = trans0;
  StdI->transindx[StdI->ntrans][0] = isite;
  StdI->transindx[StdI->ntrans][1] = ispin;
  StdI->transindx[StdI->ntrans][2] = jsite; 
  StdI->transindx[StdI->ntrans][3] = jspin;
  StdI->ntrans = StdI->ntrans + 1;
}

/**
*
* Add Hopping and Local potential for the both spin
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Hopping(
struct StdIntList *StdI,
  double complex trans0 /**< [in] Hopping integral t, mu, etc. */,
  int isite /**< [in] i for c_{i sigma}^dagger*/,
  int jsite /**< [in] j for c_{j sigma'}*/,
  int loff
  )
{
  int ispin;

  for (ispin = 0; ispin < 2; ispin++) {
    StdFace_trans(StdI, trans0, jsite, ispin, isite, ispin);
    if(loff == 1)
      StdFace_trans(StdI, conj(trans0), isite, ispin, jsite, ispin);
  }/*for (ispin = 0; ispin < 2; ispin++)*/

}

/**
*
* Add Longitudinal magnetic field to the list
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_MagField(
struct StdIntList *StdI,
  int S2 /**< [in] Spin moment in i site*/,
  double h /**< [in] Longitudinal magnetic field h. */,
  double Gamma /**< [in] Transvars magnetic field h. */,
  int isite /**< [in] i for c_{i sigma}^dagger*/)
{
  int ispin;
  double S, Sz;
 
  S = (double)S2 * 0.5;

  for (ispin = 0; ispin <= S2; ispin++){
    Sz = (double)ispin - S;
    StdFace_trans(StdI, -h * Sz, isite, ispin, isite, ispin);

    if (ispin < S2) {
      StdFace_trans(StdI, -0.5 * Gamma * sqrt(S*(S + 1.0) - Sz*(Sz + 1.0)),
        isite, ispin + 1, isite, ispin);
      StdFace_trans(StdI, -0.5 * Gamma * sqrt(S*(S + 1.0) - Sz*(Sz + 1.0)),
        isite, ispin, isite, ispin + 1);
    }
  }
}

/**
 *
 * Add interaction to the list
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_intr(
struct StdIntList *StdI,
  double complex intr0 /**< [in] Interaction U, V, J, etc.*/,
  int site1 /**< [in] i1 for c_{i1 sigma1}^dagger*/,
  int spin1 /**< [in] sigma11 for c_{i1 sigma1}^dagger*/,
  int site2 /**< [in] i2 for c_{i2 sigma2}*/,
  int spin2 /**< [in] sigma12 for c_{i2 sigma2}*/,
  int site3 /**< [in] i3 for c_{i3 sigma3}^dagger*/,
  int spin3 /**< [in] sigma13 for c_{i3 sigma3}^dagger*/,
  int site4 /**< [in] i2 for c_{i2 sigma2}*/,
  int spin4 /**< [in] sigma12 for c_{i2 sigma2}*/)
{
  StdI->intr[StdI->nintr] = intr0;
  StdI->intrindx[StdI->nintr][0] = site1; StdI->intrindx[StdI->nintr][1] = spin1;
  StdI->intrindx[StdI->nintr][2] = site2; StdI->intrindx[StdI->nintr][3] = spin2;
  StdI->intrindx[StdI->nintr][4] = site3; StdI->intrindx[StdI->nintr][5] = spin3;
  StdI->intrindx[StdI->nintr][6] = site4; StdI->intrindx[StdI->nintr][7] = spin4;
  StdI->nintr = StdI->nintr + 1;
}

/**
*
* Treat J as a 3*3 matrix [(6S + 1)*(6S' + 1) interactions]
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_GeneralJ(
struct StdIntList *StdI,
  double J[3][3],
  int Si2 /**< [in] Spin moment in i site*/,
  int Sj2 /**< [in] Spin moment in j site*/,
  int isite /**< [in] i of S_i */,
  int jsite /**< [in] j of S_j */)
{
  int ispin, jspin, ZGeneral, ExGeneral;
  double Si, Sj, Siz, Sjz;
  double complex intr0;
  /*
   Only For S=1/2 system WO off-diagonal term
  */
  ZGeneral = 1;
  ExGeneral = 1;
  if (Si2 == 1 || Sj2 == 1) {

    ZGeneral = 0;

    StdI->Hund[StdI->NHund] = -0.5 * J[2][2];
    StdI->HundIndx[StdI->NHund][0] = isite;
    StdI->HundIndx[StdI->NHund][1] = jsite;
    StdI->NHund += 1;

    StdI->Cinter[StdI->NCinter] = -0.25 * J[2][2];
    StdI->CinterIndx[StdI->NCinter][0] = isite;
    StdI->CinterIndx[StdI->NCinter][1] = jsite;
    StdI->NCinter += 1;

    if (J[0][1] < 0.000001 && J[1][0] < 0.000001) {

      ExGeneral = 0;

      StdI->Ex[StdI->NEx] = - 0.25 * (J[0][0] + J[1][1]);
      StdI->ExIndx[StdI->NEx][0] = isite;
      StdI->ExIndx[StdI->NEx][1] = jsite;
      StdI->NEx += 1;

      StdI->PairLift[StdI->NPairLift] = 0.25 * (J[0][0] - J[1][1]);
      StdI->PLIndx[StdI->NPairLift][0] = isite;
      StdI->PLIndx[StdI->NPairLift][1] = jsite;
      StdI->NPairLift += 1;
    }
  }
  /*
   For S != 1/2 spin or off-diagonal interaction
  */
  Si = 0.5 * (double)Si2;
  Sj = 0.5 * (double)Sj2;

  for (ispin = 0; ispin <= Si2; ispin++) {
    Siz = (double)ispin - Si;
    for (jspin = 0; jspin <= Sj2; jspin++) {
      Sjz = (double)jspin - Sj;
      /*
       S_{i z} * S_{j z}
      */
      if (ZGeneral == 1) {
        intr0 = J[2][2] * Siz * Sjz;
        StdFace_intr(StdI, intr0,
          isite, ispin, isite, ispin, jsite, jspin, jsite, jspin);
      }
      /*
       S_i^+ S_j^- + S_j^+ S_i^-
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
      /*
       S_i^+ S_j^+ + S_j^- S_i^-
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
      /*
       S_i^+ S_{j z} + S_{j z} S_i^-
      */
      if (ispin < Si2) {
        intr0 = 0.5 * (J[0][2] - I * J[1][2]) * sqrt(Si * (Si + 1.0) - Siz * (Siz + 1.0)) * Sjz;
        StdFace_intr(StdI, intr0,
          isite, ispin + 1, isite, ispin, jsite, jspin, jsite, jspin);
        StdFace_intr(StdI, conj(intr0),
          jsite, jspin, jsite, jspin, isite, ispin, isite, ispin + 1);
      }/*if (ispin < Si2)*/
      /*
       S_{i z} S_j^+ + S_j^- S_{i z}
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
 *
 * Add onsite/offsite Coulomb term to the list
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_Coulomb(
struct StdIntList *StdI,
  double V /**< [in] Coulomb integral U, V, etc.*/,
  int isite /**< [in] i of n_i */,
  int jsite /**< [in] j of n_j */)
{
  StdI->Cinter[StdI->NCinter] = V;
  StdI->CinterIndx[StdI->NCinter][0] = isite;
  StdI->CinterIndx[StdI->NCinter][1] = jsite;
  StdI->NCinter += 1;
}

/**
 *
 * Print a valiable (real) read from the input file
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_PrintVal_d(
  char* valname /**< [in] Name of the valiable*/, 
  double *val /**< [inout] Valiable to be set*/, 
  double val0 /**< [in] The default value*/)
{
  if (isnan(*val) == 1) {
    *val = val0;
    fprintf(stdout, "  %15s = %-10.5f  ######  DEFAULT VALUE IS USED  ######\n", valname, *val);
  }
  else fprintf(stdout, "  %15s = %-10.5f\n", valname, *val);
}

/**
*
* Print a valiable (real) read from the input file
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_PrintVal_dd(
  char* valname /**< [in] Name of the valiable*/,
  double *val /**< [inout] Valiable to be set*/,
  double val0 /**< [in] The default value*/,
  double val1 /**< [in] The default value*/)
{
  if (isnan(*val) == 1) {
    if (isnan(val0) == 1) *val = val1;
    else *val = val0;
    fprintf(stdout, "  %15s = %-10.5f  ######  DEFAULT VALUE IS USED  ######\n", valname, *val);
  }
  else fprintf(stdout, "  %15s = %-10.5f\n", valname, *val);
}

/**
*
* Print a valiable (complex) read from the input file
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_PrintVal_c(
  char* valname /**< [in] Name of the valiable*/,
  double complex *val /**< [inout] Valiable to be set*/,
  double complex val0 /**< [in] The default value*/)
{
  if (isnan(creal(*val)) == 1) {
    *val = val0;
    fprintf(stdout, "  %15s = %-10.5f %-10.5f  ######  DEFAULT VALUE IS USED  ######\n", valname, creal(*val), cimag(*val));
  }
  else fprintf(stdout, "  %15s = %-10.5f %-10.5f\n", valname, creal(*val), cimag(*val));
}

/**
 *
 * Print a valiable (integer) read from the input file
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_PrintVal_i(
  char* valname /**< [in] Name of the valiable*/,
  int *val /**< [inout] Valiable to be set*/,
  int val0 /**< [in] The default value*/)
{
  int NaN_i = 2147483647;

  if (*val == NaN_i) {
    *val = val0;
    fprintf(stdout, "  %15s = %-10d  ######  DEFAULT VALUE IS USED  ######\n", valname, *val);
  }
  else fprintf(stdout, "  %15s = %-10d\n", valname, *val);
}

/**
 *
 * Stop HPhi if a variable (real) not used is specified
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_NotUsed_d(
  char* valname /**< [in] Name of the valiable*/,
  double val /**< [in]*/)
{
  if (isnan(val) == 0) {
    fprintf(stdout, "\n Check !  %s is SPECIFIED but will NOT be USED. \n", valname);
    fprintf(stdout, "            Please COMMENT-OUT this line \n");
    fprintf(stdout, "            or check this input is REALLY APPROPRIATE for your purpose ! \n\n");
    StdFace_exit(-1);
  }
}

/**
*
* Stop HPhi if a variable (real) not used is specified
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_NotUsed_c(
  char* valname /**< [in] Name of the valiable*/,
  double complex val /**< [in]*/)
{
  if (isnan(creal(val)) == 0) {
    fprintf(stdout, "\n Check !  %s is SPECIFIED but will NOT be USED. \n", valname);
    fprintf(stdout, "            Please COMMENT-OUT this line \n");
    fprintf(stdout, "            or check this input is REALLY APPROPRIATE for your purpose ! \n\n");
    StdFace_exit(-1);
  }
}

/**
*
* Stop HPhi if variables (real) not used is specified
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_NotUsed_J(
  char* valname /**< [in] Name of the valiable*/,
  double JAll /**< [in]*/,
  double J[3][3] /**< [in]*/)
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

 }

/**
 *
 * Stop HPhi if a variable (integer) not used is specified
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_NotUsed_i(
  char* valname /**< [in] Name of the valiable*/,
  int val /**< [in]*/)
{
  int NaN_i = 2147483647;

  if (val != NaN_i) {
    fprintf(stdout, "\n Check !  %s is SPECIFIED but will NOT be USED. \n", valname);
    fprintf(stdout, "            Please COMMENT-OUT this line \n");
    fprintf(stdout, "            or check this input is REALLY APPROPRIATE for your purpose ! \n\n");
    StdFace_exit(-1);
  }
}

/**
 *
 * Stop HPhi if a variable (integer) which must be specified
 * is absent in the input file.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_RequiredVal_i(
  char* valname /**< [in] Name of the valiable*/,
  int val /**< [in]*/)
{
  int NaN_i = 2147483647;

  if (val == NaN_i){
    fprintf(stdout, "ERROR ! %s is NOT specified !\n", valname);
    StdFace_exit(-1);
  }
  else fprintf(stdout, "  %15s = %-3d\n", valname, val);
}
/**
*
* Define whether the specified site is in the unit cell or not.
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StdFace_FoldSite(struct StdIntList *StdI,
  int iCellV[3], int nBox[3], int iCellV_fold[3])
{
  int ii, jj, iCellV_frac[3];
  /*
  Transform to fractional coordinate (times NCell)
  */
  for (ii = 0; ii < 3; ii++) {
    iCellV_frac[ii] = 0;
    for (jj = 0; jj < 3; jj++)iCellV_frac[ii] += StdI->rbox[ii][jj] * iCellV[jj];
  }
  /*
  Which supercell contains this cell
  */
  for (ii = 0; ii < 3; ii++)
    nBox[ii] = (iCellV_frac[ii] + StdI->NCell * 1000) / StdI->NCell - 1000;
  /*
  Fractional coordinate (times NCell) in the original supercell
  */
  for (ii = 0; ii < 3; ii++)
    iCellV_frac[ii] -= StdI->NCell*(nBox[ii]);
  /**/
  for (ii = 0; ii < 3; ii++) {
    iCellV_fold[ii] = 0;
    for (jj = 0; jj < 3; jj++) iCellV_fold[ii] += StdI->box[jj][ii] * iCellV_frac[jj];
    iCellV_fold[ii] = (iCellV_fold[ii] + StdI->NCell * 1000) / StdI->NCell - 1000;
  }
}
/**
*
* Initialize Cell
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_InitSite(struct StdIntList *StdI, FILE *fp, int dim)
{
  int bound[3][2], edge, ii, jj;
  int ipos;
  int nBox[3], iCellV_fold[3], iCellV[3];
  double pos[4][2], xmin, xmax/*, offset[2], scale*/;
  /*
   check Input parameters
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
  /*
   Define the phase factor at the boundary
  */
  if (dim == 2) StdI->phase[2] = 0.0;
  for (ii = 0; ii < 3; ii++) {
    StdI->ExpPhase[ii] = cos(StdI->pi180 * StdI->phase[ii]) + I*sin(StdI->pi180 * StdI->phase[ii]);
    if (cabs(StdI->ExpPhase[ii] + 1.0) < 0.000001) StdI->AntiPeriod[ii] = 1;
    else StdI->AntiPeriod[ii] = 0;
  }
  /*
   Structure in a cell
  */
  StdI->tau = (double **)malloc(sizeof(double*) * StdI->NsiteUC);
  for (ii = 0; ii < StdI->NsiteUC; ii++) {
    StdI->tau[ii] = (double *)malloc(sizeof(double) * 3);
  }
  /*
   Calculate reciprocal lattice vectors (times NCell)
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
  printf("         Number of Cell : %d\n", abs(StdI->NCell));
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
  /*
   Find the lower- and the upper bound for surrowndinf the simulation box
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
  /*
   Find Cells in the Simulation Box
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
  /*
   Initialize gnuplot
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
  }
}/*void StdFace_InitSite2D*/
/**
 * Find the index of transfer and interaction
 */
void StdFace_FindSite(struct StdIntList *StdI,
  int iW, int iL, int iH, int diW, int diL, int diH,
  int isiteUC, int jsiteUC,
  int *isite, int *jsite, double complex *Cphase)
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
 /*
 * Set Label in the gnuplot display
 */
void StdFace_SetLabel(struct StdIntList *StdI, FILE *fp, 
  int iW, int iL, int diW, int diL, int isiteUC, int jsiteUC, 
  int *isite, int *jsite, int connect, double complex *Cphase)
{
  int iCell, jCell, kCell;
  int jCellV[3], nBox[2], jCellV_fold[3];
  double xi, yi, xj, yj;
  /*
   Reversed
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
  /*
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
 * Print lattice.xsf (XCrysDen format) 
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
  fclose(fp);
}/*void StdFace_PrintXSF*/
/*
 * Input nearest-neighbor spin-spin interaction
 */
void StdFace_InputSpinNN(struct StdIntList *StdI, double J0[3][3], 
  double J0All, char *J0name) 
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

}

void StdFace_InputSpin(struct StdIntList *StdI, double Jp[3][3],
  double JpAll, char *Jpname)
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

}

void StdFace_InputCoulombV(struct StdIntList *StdI, double *V0, char *V0name)
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

}

void StdFace_InputHopp(struct StdIntList *StdI, double complex *t0, char *t0name)
{

  if (isnan(creal(StdI->t)) == 0 && isnan(creal(*t0)) == 0) {
    fprintf(stdout, "\n ERROR! t and %s conflict !\n\n", t0name);
    StdFace_exit(-1);
  }
  else if (isnan(creal(*t0)) == 0)
    fprintf(stdout, "  %15s = %-10.5f\n", t0name, creal(*t0));
  else if (isnan(creal(StdI->t)) == 0) {
    *t0 = StdI->t;
    fprintf(stdout, "  %15s = %-10.5f\n", t0name, creal(*t0));
  }
  else {
    *t0 = 0.0;
  }

}/*void StdFace_InputHopp*/
/*
 Print geometry of sites for the pos-process of correlation function
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
  fclose(fp);

}/*void StdFace_PrintGeometry()*/
/*
 * Malloc Arrays for interactions
 */
void StdFace_MallocInteractions(struct StdIntList *StdI) {
  int ii;
  /*
   Transfer
  */
  StdI->transindx = (int **)malloc(sizeof(int*) * StdI->ntrans);
  StdI->trans = (double complex *)malloc(sizeof(double complex) * StdI->ntrans);
  for (ii = 0; ii < StdI->ntrans; ii++) {
    StdI->transindx[ii] = (int *)malloc(sizeof(int) * 4);
  }
  /*
   InterAll
  */
  StdI->intrindx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->intr = (double complex *)malloc(sizeof(double complex) * StdI->nintr);
  for (ii = 0; ii < StdI->nintr; ii++) {
    StdI->intrindx[ii] = (int *)malloc(sizeof(int) * 8);
  }
  /*
  Coulomb intra
  */
  StdI->CintraIndx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->Cintra = (double *)malloc(sizeof(double) * StdI->nintr);
  for (ii = 0; ii < StdI->nintr; ii++) {
    StdI->CintraIndx[ii] = (int *)malloc(sizeof(int) * 1);
  }
  /*
  Coulomb inter
  */
  StdI->CinterIndx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->Cinter = (double *)malloc(sizeof(double) * StdI->nintr);
  for (ii = 0; ii < StdI->nintr; ii++) {
    StdI->CinterIndx[ii] = (int *)malloc(sizeof(int) * 2);
  }
  /*
  Hund
  */
  StdI->HundIndx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->Hund = (double *)malloc(sizeof(double) * StdI->nintr);
  for (ii = 0; ii < StdI->nintr; ii++) {
    StdI->HundIndx[ii] = (int *)malloc(sizeof(int) * 2);
  }
  /*
  Excahnge
  */
  StdI->ExIndx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->Ex = (double *)malloc(sizeof(double) * StdI->nintr);
  for (ii = 0; ii < StdI->nintr; ii++) {
    StdI->ExIndx[ii] = (int *)malloc(sizeof(int) * 2);
  }
  /*
  PairLift
  */
  StdI->PLIndx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->PairLift = (double *)malloc(sizeof(double) * StdI->nintr);
  for (ii = 0; ii < StdI->nintr; ii++) {
    StdI->PLIndx[ii] = (int *)malloc(sizeof(int) * 2);
  }

  StdI->NCintra = 0;
  StdI->NCinter = 0;
  StdI->NHund = 0;
  StdI->NEx = 0;
  StdI->NPairLift = 0;
}/*void StdFace_MallocInteractions*/

#if defined(_mVMC)
 /**
 *
 * Define whether the specified site is in the unit cell or not.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void StdFace_FoldSiteSub(struct StdIntList *StdI,
  int iCellV[3], int nBox[3], int iCellV_fold[3])
{
  int ii, jj, iCellV_frac[3];
  /*
  Transform to fractional coordinate (times NCell)
  */
  for (ii = 0; ii < 3; ii++) {
    iCellV_frac[ii] = 0;
    for (jj = 0; jj < 3; jj++)iCellV_frac[ii] += StdI->rboxsub[ii][jj] * iCellV[jj];
  }
  /*
  Which supercell contains this cell
  */
  for (ii = 0; ii < 3; ii++)
    nBox[ii] = (iCellV_frac[ii] + StdI->NCellsub * 1000) / StdI->NCellsub - 1000;
  /*
  Fractional coordinate (times NCell) in the original supercell
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
*
* Print Quantum number projection
*
* @author Mitsuaki Kawamura (The University of Tokyo)
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
  /*
  Define translation operator in sub lattice
  */
  StdI->NSym = 0;
  for (iCell = 0; iCell < StdI->NCell; iCell++) {

    StdFace_FoldSiteSub(StdI, StdI->Cell[iCell], nBox, iCellV);

    StdFace_FoldSite(StdI, iCellV, nBox, iCellV);

    if (iCellV[0] == StdI->Cell[iCell][0] && 
        iCellV[1] == StdI->Cell[iCell][1] && 
        iCellV[2] == StdI->Cell[iCell][2]) {
      /*
      Translation operator in sub lattice
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
*
* Initialize sub Cell
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
static void StdFace_InitSiteSub(struct StdIntList *StdI)
{
  int ii, jj, kk, prod;
  int bound[3][2], iCellV[3], nBox[3], iCellV_fold[3];
  /*
  check Input parameters
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
  /*
  Calculate reciprocal lattice vectors (times NCellsub)
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
  /*
   Check : Is the sublattice commensurate ?
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
/*
 *Generate orbitalindex
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
}
#endif
