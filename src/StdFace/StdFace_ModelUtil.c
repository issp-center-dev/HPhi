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
  if (*val > 9999.0) {
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
  if (*val > 9999.0) {
    if (val0 > 9999.0) *val = val1;
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
  if (creal(*val) > 9999.0) {
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
  if (*val == 9999) {
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
  if (val < 9999.0) {
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
  if (creal(val) < 9999.0) {
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
  if (val != 9999) {
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
  if (val == 9999){
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
static void StdFace_FoldSite2D(struct StdIntList *StdI,
  int iCellV[2], int nBox[2], int iCellV_fold[2])
{
  int iCellV_frac[2], ii;
  /*
   Transform to fractional coordinate (times NCell)
  */
  for (ii = 0; ii < 2; ii++)
    iCellV_frac[ii] = StdI->rbox[ii][0] * iCellV[0] + StdI->rbox[ii][1] * iCellV[1];
  /*
   Which supercell contains this cell
  */
  for (ii = 0; ii < 2; ii++)
    nBox[ii] = (iCellV_frac[ii] + StdI->NCell * 1000) / StdI->NCell - 1000;
  /*
   Fractional coordinate (times NCell) in the original supercell
  */
  for (ii = 0; ii < 2; ii++)
    iCellV_frac[ii] -= StdI->NCell*(nBox[ii]);
  /**/
  for (ii = 0; ii < 2; ii++) {
    iCellV_fold[ii] = StdI->box[0][ii] * iCellV_frac[0] + StdI->box[1][ii] * iCellV_frac[1];
    iCellV_fold[ii] = (iCellV_fold[ii] + StdI->NCell * 1000) / StdI->NCell - 1000;
  }
}

/**
*
* Initialize Cell
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_InitSite2D(struct StdIntList *StdI, FILE *fp)
{
  int bound[2][2], edge, ii, jj;
  int ipos;
  int nBox[2], iCellV_fold[2], iCellV[2];
  double pos[4][2], xmin, xmax/*, offset[2], scale*/;
  /*
   check Input parameters
  */
  if ((StdI->L != 9999 || StdI->W != 9999)
    && (StdI->box[0][1] != 9999 || StdI->box[0][0] != 9999 || StdI->box[1][1] != 9999 || StdI->box[1][0] != 9999)) {
    fprintf(stdout, "\nERROR ! (L, W) and (a0W, a0L, a1W, a1L) conflict !\n\n");
    StdFace_exit(-1);
  }
  else if (StdI->L != 9999 || StdI->W != 9999) {
    if (StdI->L == 9999 || StdI->W == 9999) {
      fprintf(stdout, "\nERROR ! Both W and L MUST be specified !\n\n");
      StdFace_exit(-1);
    }
    StdFace_PrintVal_i("L", &StdI->L, 1);
    StdFace_PrintVal_i("W", &StdI->W, 1);
    StdI->box[0][0] = StdI->W;
    StdI->box[0][1] = 0;
    StdI->box[1][0] = 0;
    StdI->box[1][1] = StdI->L;
  }
  else if (StdI->box[0][0] != 9999 || StdI->box[0][1] != 9999 || 
           StdI->box[1][0] != 9999 || StdI->box[1][1] != 9999) 
  {
    if (StdI->box[0][0] == 9999 || StdI->box[0][1] == 9999 ||
        StdI->box[1][0] == 9999 || StdI->box[1][1] == 9999) {
      fprintf(stdout, "\nERROR ! Both a0W, a0L, a1W, a1L MUST be specified !\n\n");
      StdFace_exit(-1);
    }
    StdFace_PrintVal_i("a0W", &StdI->box[0][0], 1);
    StdFace_PrintVal_i("a0L", &StdI->box[0][1], 0);
    StdFace_PrintVal_i("a1W", &StdI->box[1][0], 0);
    StdFace_PrintVal_i("a1L", &StdI->box[1][1], 1);
  }
  /*
   Structure in a cell
  */
  StdI->tau = (double **)malloc(sizeof(double*) * StdI->NsiteUC);
  for (ii = 0; ii < StdI->NsiteUC; ii++) {
    StdI->tau[ii] = (double *)malloc(sizeof(double) * 2);
  }
  /*
   Calculate reciprocal lattice vectors (times NCell)
  */
  StdI->NCell = StdI->box[0][0] * StdI->box[1][1] - StdI->box[0][1] * StdI->box[1][0];
  printf("         Number of Cell : %d\n", abs(StdI->NCell));
  if (StdI->NCell == 0) {
    StdFace_exit(-1);
  }

  StdI->rbox[0][0] = StdI->box[1][1];
  StdI->rbox[0][1] = - StdI->box[1][0];
  StdI->rbox[1][0] = - StdI->box[0][1];
  StdI->rbox[1][1] = StdI->box[0][0];
  if (StdI->NCell < 0) {
    for (ii = 0; ii < 2; ii++) 
      for (jj = 0; jj < 2; jj++)
        StdI->rbox[ii][jj] *= -1;
    StdI->NCell *= -1;
  }/*if (StdI->NCell < 0)*/
  /*
   Find the lower- and the upper bound for surrowndinf the simulation box
  */
  for (ii = 0; ii < 2; ii++) {
    bound[ii][0] = 0;
    bound[ii][1] = 0;
    for (nBox[1] = 0; nBox[1] < 2; nBox[1]++) {
      for (nBox[0] = 0; nBox[0] < 2; nBox[0]++) {
        edge = nBox[0] * StdI->box[0][ii] + nBox[1] * StdI->box[1][ii];
        if (edge < bound[ii][0]) bound[ii][0] = edge;
        if (edge > bound[ii][1]) bound[ii][1] = edge;
      }
    }
  }
  /*
   Find Cells in the Simulation Box
  */
  StdI->Cell = (int **)malloc(sizeof(int*) * StdI->NCell);
  for (ii = 0; ii < StdI->NCell; ii++) {
    StdI->Cell[ii] = (int *)malloc(sizeof(int) * 2);
  }/*for (iCell = 0; iCell < (Wmax - Wmin + 1) * (Lmax - Lmin + 1); iCell++)*/

  jj = 0;
  for (iCellV[1] = bound[1][0]; iCellV[1] <= bound[1][1]; iCellV[1]++) {
    for (iCellV[0] = bound[0][0]; iCellV[0] <= bound[0][1]; iCellV[0]++) {
      StdFace_FoldSite2D(StdI, iCellV, nBox, iCellV_fold);
      if (nBox[0] == 0 && nBox[1] == 0) {
        StdI->Cell[jj][0] = iCellV[0];
        StdI->Cell[jj][1] = iCellV[1];
        jj += 1;
      }/*if (lUC == 1)*/
    }/*for (iW = Wmin; iW <= Wmax; iW++*/
  }/*for (iL = Lmin; iL <= Lmax; iL++)*/
  /*
   Initialize gnuplot
  */
  pos[0][0] = 0.0;
  pos[0][1] = 0.0;
  pos[1][0] = StdI->direct[0][0] * (double)StdI->box[0][0] + StdI->direct[1][0] * (double)StdI->box[0][1];
  pos[1][1] = StdI->direct[0][1] * (double)StdI->box[0][0] + StdI->direct[1][1] * (double)StdI->box[0][1];
  pos[2][0] = StdI->direct[0][0] * (double)StdI->box[1][0] + StdI->direct[1][0] * (double)StdI->box[1][1];
  pos[2][1] = StdI->direct[0][1] * (double)StdI->box[1][0] + StdI->direct[1][1] * (double)StdI->box[1][1];
  pos[3][0] = pos[1][0] + pos[2][0];
  pos[3][1] = pos[1][1] + pos[2][1];
  /**/
  /*
  scale = sqrt((double)((StdI->box[0][0] + StdI->box[1][0])*(StdI->box[0][0] + StdI->box[1][0])
                      + (StdI->box[0][1] + StdI->box[1][1])*(StdI->box[0][1] + StdI->box[1][1])));
  scale = 0.5 / scale;
  offset[0] = pos[3][0] * scale;
  offset[1] = pos[3][1] * scale;
  
  for (ipos = 0; ipos < 4; ipos++) {
    pos[ipos][0] -= offset[0];
    pos[ipos][1] -= offset[1];
  }
  */
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

}/*void StdFace_InitSite2D*/
/*
 * Set Label in the gnuplot display
 */
void StdFace_SetLabel(struct StdIntList *StdI, FILE *fp, 
  int iW, int iL, int diW, int diL, int isiteUC, int jsiteUC, 
  int *isite, int *jsite, int connect, double complex *phase)
{
  int iCell, jCell, kCell;
  int nBox[2];
  int jCellV_fold[2], jCellV[2];
  double xi, yi, xj, yj;
  /*
   Reversed
  */
  jCellV[0] = iW - diW;
  jCellV[1] = iL - diL;
  StdFace_FoldSite2D(StdI, jCellV, nBox, jCellV_fold);
  /**/
  for (kCell = 0; kCell < StdI->NCell; kCell++) {
    if (jCellV_fold[0] == StdI->Cell[kCell][0] && jCellV_fold[1] == StdI->Cell[kCell][1]) {
      jCell = kCell;
    }
    if (iW == StdI->Cell[kCell][0] && iL == StdI->Cell[kCell][1]) {
      iCell = kCell;
    }
  }/*for (iCell = 0; iCell < StdI->NCell; iCell++)*/
  *isite = iCell * StdI->NsiteUC + jsiteUC;
  *jsite = jCell * StdI->NsiteUC + isiteUC;
  if (strcmp(StdI->model, "kondo") == 0) {
    *isite += StdI->NCell * StdI->NsiteUC;
    *jsite += StdI->NCell * StdI->NsiteUC;
  }

  xi = StdI->direct[0][0] * ((double)iW + StdI->tau[jsiteUC][0])
     + StdI->direct[1][0] * ((double)iL + StdI->tau[jsiteUC][1]);
  yi = StdI->direct[0][1] * ((double)iW + StdI->tau[jsiteUC][0])
     + StdI->direct[1][1] * ((double)iL + StdI->tau[jsiteUC][1]);

  xj = StdI->direct[0][0] * ((double)jCellV[0] + StdI->tau[isiteUC][0])
     + StdI->direct[1][0] * ((double)jCellV[1] + StdI->tau[isiteUC][1]);
  yj = StdI->direct[0][1] * ((double)jCellV[0] + StdI->tau[isiteUC][0])
     + StdI->direct[1][1] * ((double)jCellV[1] + StdI->tau[isiteUC][1]);

  if (*isite < 10)fprintf(fp, "set label \"%1d\" at %f, %f center front\n", *isite, xi, yi);
  else fprintf(fp, "set label \"%2d\" at %f, %f center front\n", *isite, xi, yi);
  if (*jsite < 10)fprintf(fp, "set label \"%1d\" at %f, %f center front\n", *jsite, xj, yj);
  else fprintf(fp, "set label \"%2d\" at %f, %f center front\n", *jsite, xj, yj);
  fprintf(fp, "set arrow from %f, %f to %f, %f nohead ls %d\n", xi, yi, xj, yj, connect);
  /*
  */
  jCellV[0] = iW + diW;
  jCellV[1] = iL + diL;
  StdFace_FoldSite2D(StdI, jCellV, nBox, jCellV_fold);
  *phase = cpow(StdI->ExpPhase0, (double)nBox[0]) * cpow(StdI->ExpPhase1, (double)nBox[1]);
  /**/
  for (kCell = 0; kCell < StdI->NCell; kCell++) {
    if (jCellV_fold[0] == StdI->Cell[kCell][0] && jCellV_fold[1] == StdI->Cell[kCell][1]) {
      jCell = kCell;
    }
    if (iW == StdI->Cell[kCell][0] && iL == StdI->Cell[kCell][1]) {
      iCell = kCell;
    }
  }/*for (iCell = 0; iCell < StdI->NCell; iCell++)*/
  *isite = iCell * StdI->NsiteUC + isiteUC;
  *jsite = jCell * StdI->NsiteUC + jsiteUC;
  if (strcmp(StdI->model, "kondo") == 0) {
    *isite += StdI->NCell * StdI->NsiteUC;
    *jsite += StdI->NCell * StdI->NsiteUC;
  }

  xi = StdI->direct[1][0] * ((double)iL + StdI->tau[isiteUC][1])
     + StdI->direct[0][0] * ((double)iW + StdI->tau[isiteUC][0]);
  yi = StdI->direct[1][1] * ((double)iL + StdI->tau[isiteUC][1])
     + StdI->direct[0][1] * ((double)iW + StdI->tau[isiteUC][0]);

  xj = StdI->direct[0][0] * ((double)jCellV[0] + StdI->tau[jsiteUC][0])
     + StdI->direct[1][0] * ((double)jCellV[1] + StdI->tau[jsiteUC][1]);
  yj = StdI->direct[0][1] * ((double)jCellV[0] + StdI->tau[jsiteUC][0])
     + StdI->direct[1][1] * ((double)jCellV[1] + StdI->tau[jsiteUC][1]);

  if(*isite < 10)fprintf(fp, "set label \"%1d\" at %f, %f center front\n", *isite, xi, yi);
  else fprintf(fp, "set label \"%2d\" at %f, %f center front\n", *isite, xi, yi);
  if (*jsite < 10)fprintf(fp, "set label \"%1d\" at %f, %f center front\n", *jsite, xj, yj);
  else fprintf(fp, "set label \"%2d\" at %f, %f center front\n", *jsite, xj, yj);
  fprintf(fp, "set arrow from %f, %f to %f, %f nohead ls %d\n", xi, yi, xj, yj, connect);
}/*void StdFace_SetLabel*/
 /**
 *
 * Define whether the specified site is in the unit cell or not.
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void StdFace_FoldSite3D(struct StdIntList *StdI,
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
void StdFace_InitSite3D(struct StdIntList *StdI, FILE *fp)
{
  int ii, jj, edge;
  int bound[3][2], iCellV[3], nBox[3], iCellV_fold[3];
  /*
  check Input parameters
  */
  if ((StdI->L != 9999 || StdI->W != 9999 || StdI->Height != 9999)
    && (StdI->box[0][1] != 9999 || StdI->box[0][0] != 9999 || StdI->box[0][2] != 9999 || 
        StdI->box[1][1] != 9999 || StdI->box[1][0] != 9999 || StdI->box[1][2] != 9999 ||
        StdI->box[2][1] != 9999 || StdI->box[2][0] != 9999 || StdI->box[2][2] != 9999)) 
  {
    fprintf(stdout, "\nERROR ! (L, W, Height) and (a0W, ..., a2H) conflict !\n\n");
    StdFace_exit(-1);
  }
  else if (StdI->W != 9999 || StdI->L != 9999 || StdI->Height != 9999) {
    if (StdI->L == 9999 || StdI->W == 9999 || StdI->Height == 9999) {
      fprintf(stdout, "\nERROR ! All of W, L, Height MUST be specified !\n\n");
      StdFace_exit(-1);
    }
    StdFace_PrintVal_i("L", &StdI->L, 1);
    StdFace_PrintVal_i("W", &StdI->W, 1);
    StdFace_PrintVal_i("Height", &StdI->Height, 1);
    for (ii = 0; ii < 3; ii++) for (jj = 0; jj < 3; jj++)
      StdI->box[ii][jj] = 0;
    StdI->box[0][0] = StdI->W;
    StdI->box[1][1] = StdI->L;
    StdI->box[2][2] = StdI->Height;
  }
  else if (
    StdI->box[0][0] != 9999 || StdI->box[0][1] != 9999 || StdI->box[0][2] != 9999 || 
    StdI->box[1][0] != 9999 || StdI->box[1][1] != 9999 || StdI->box[1][2] != 9999 ||
    StdI->box[2][0] != 9999 || StdI->box[2][1] != 9999 || StdI->box[2][2] != 9999) {
    if (StdI->box[0][0] == 9999 || StdI->box[0][1] == 9999 || StdI->box[0][2] == 9999 ||
        StdI->box[1][0] == 9999 || StdI->box[1][1] == 9999 || StdI->box[1][2] == 9999 ||
        StdI->box[2][0] == 9999 || StdI->box[2][1] == 9999 || StdI->box[2][2] == 9999) {
      fprintf(stdout, "\nERROR ! All of a0W, a0L, a0H, a1W, ... MUST be specified !\n\n");
      StdFace_exit(-1);
    }
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
   Initialize gnuplot
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
  Calculate the number of Unit Cell
  */
  StdI->Cell = (int **)malloc(sizeof(int*) * StdI->NCell);
  for (ii = 0; ii < StdI->NCell; ii++) {
    StdI->Cell[ii] = (int *)malloc(sizeof(int) * 3);
  }/*for (ii = 0; ii < (Wmax - Wmin + 1) * (Lmax - Lmin + 1); ii++)*/

  jj = 0;
  for (iCellV[2] = bound[2][0]; iCellV[2] <= bound[2][1]; iCellV[2]++) {
    for (iCellV[1] = bound[1][0]; iCellV[1] <= bound[1][1]; iCellV[1]++) {
      for (iCellV[0] = bound[0][0]; iCellV[0] <= bound[0][1]; iCellV[0]++) {
        StdFace_FoldSite3D(StdI, iCellV, nBox, iCellV_fold);
        if (nBox[0] == 0 && nBox[1] == 0 && nBox[2] == 0) {
          for (ii = 0; ii < 3; ii++)
            StdI->Cell[jj][ii] = iCellV[ii];
          jj += 1;
        }/*if (lUC == 1)*/
      }/*for (iCellV[0] = bound[0][0]; iCellV[0] <= bound[0][1]; iCellV[0]++*/
    }/*for (iCellV[1] = bound[1][0]; iCellV[1] <= bound[1][1]; iCellV[1]++)*/
  }/*for (iCellV[2] = bound[2][0]; iCellV[2] <= bound[2][1]; iCellV[2]++)*/
}/*void StdFace_InitSite3D*/
 /**
  * Set Label in the gnuplot display
  */
void StdFace_FindSite3d(struct StdIntList *StdI,
  int iW, int iL, int iH, int diW, int diL, int diH, 
  int isiteUC, int jsiteUC,
  int *isite, int *jsite, double complex *phase)
{
  int iCell, jCell, kCell;
  int nBox[3];
  int jCellV_fold[3], jCellV[3];
  /**/
  jCellV[0] = iW + diW;
  jCellV[1] = iL + diL;
  jCellV[2] = iH + diH;
  StdFace_FoldSite2D(StdI, jCellV, nBox, jCellV_fold);
  *phase = cpow(StdI->ExpPhase0, (double)nBox[0]) 
         * cpow(StdI->ExpPhase1, (double)nBox[1])
         * cpow(StdI->ExpPhase2, (double)nBox[2]);
  /**/
  for (kCell = 0; kCell < StdI->NCell; kCell++) {
    if (jCellV_fold[0] == StdI->Cell[kCell][0] && 
        jCellV_fold[1] == StdI->Cell[kCell][1] &&
        jCellV_fold[2] == StdI->Cell[kCell][2]) {
      jCell = kCell;
    }
    if (iW == StdI->Cell[kCell][0] && 
        iL == StdI->Cell[kCell][1] &&
        iH == StdI->Cell[kCell][2]) {
      iCell = kCell;
    }
  }/*for (iCell = 0; iCell < StdI->NCell; iCell++)*/
  *isite = iCell * StdI->NsiteUC + isiteUC;
  *jsite = jCell * StdI->NsiteUC + jsiteUC;
  if (strcmp(StdI->model, "kondo") == 0) {
    *isite += StdI->NCell * StdI->NsiteUC;
    *jsite += StdI->NCell * StdI->NsiteUC;
  }
}/*void StdFace_FindSite3d*/
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
        vec[jj] += StdI->direct[kk][jj] * (double)StdI->box[ii][jj];
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
          vec[jj] += StdI->direct[kk][jj] * 
          ((double)StdI->Cell[iCell][kk] + StdI->tau[isite][kk]);
      }
      fprintf(fp, "X %15.5f %15.5f %15.5f\n", vec[0], vec[1], vec[2]);
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

  if (StdI->JAll < 9999.0 && J0All < 9999.0) {
    fprintf(stdout, "\n ERROR! J and %s conflict !\n\n", J0name);
    StdFace_exit(-1);
  }
  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      if (StdI->JAll < 9999.0 && StdI->J[i1][i2] < 9999.0) {
        fprintf(stdout, "\n ERROR! J and J%s conflict !\n\n", Jname[i1][i2]);
        StdFace_exit(-1);
      }
      else if (J0All < 9999.0 && StdI->J[i1][i2] < 9999.0) {
        fprintf(stdout, "\n ERROR! %s and J%s conflict !\n\n",
          J0name, Jname[i1][i2]);
        StdFace_exit(-1);
      }
      else if (J0All < 9999.0 && J0[i1][i2] < 9999.0) {
        fprintf(stdout, "\n ERROR! %s and %s%s conflict !\n\n", J0name,
          J0name, Jname[i1][i2]);
        StdFace_exit(-1);
      }
      else if (J0[i1][i2] < 9999.0 && StdI->JAll < 9999.0) {
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
          if (J0[i1][i2] < 9999.0 && StdI->J[i3][i4] < 9999.0) {
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
      if (J0[i1][i2] < 9999.0)
        fprintf(stdout, "  %14s%s = %-10.5f\n", J0name, Jname[i1][i2], J0[i1][i2]);
      else if (StdI->J[i1][i2] < 9999.0) {
        J0[i1][i2] = StdI->J[i1][i2];
        fprintf(stdout, "  %14s%s = %-10.5f\n", J0name, Jname[i1][i2], J0[i1][i2]);
      }
      else if (i1 == i2 && J0All < 9999.0) {
        J0[i1][i2] = J0All;
        fprintf(stdout, "  %14s%s = %-10.5f\n", J0name, Jname[i1][i2], J0[i1][i2]);
      }
      else if (i1 == i2 && StdI->JAll < 9999.0) {
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
      if (JpAll < 9999.0 && Jp[i1][i2] < 9999.0) {
        fprintf(stdout, "\n ERROR! %s and %s%s conflict !\n\n", Jpname,
          Jpname, Jname[i1][i2]);
        StdFace_exit(-1);
      }
    }/*for (j = 0; j < 3; j++)*/
  }/*for (i = 0; i < 3; i++)*/

  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      if (Jp[i1][i2] < 9999.0)
        fprintf(stdout, "  %14s%s = %-10.5f\n", Jpname, Jname[i1][i2], Jp[i1][i2]);
      else if (i1 == i2 && JpAll < 9999.0) {
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
  
  if (StdI->V < 9999.0 && *V0 < 9999.0) {
    fprintf(stdout, "\n ERROR! V and %s conflict !\n\n", V0name);
    StdFace_exit(-1);
  }
  else if (*V0 < 9999.0)
    fprintf(stdout, "  %15s = %-10.5f\n", V0name, *V0);
  else if (StdI->V < 9999.0) {
    *V0 = StdI->V;
    fprintf(stdout, "  %15s = %-10.5f\n", V0name, *V0);
  }
  else {
    *V0 = 0.0;
  }

}

void StdFace_InputHopp(struct StdIntList *StdI, double complex *t0, char *t0name)
{

  if (creal(StdI->t) < 9999.0 && creal(*t0) < 9999.0) {
    fprintf(stdout, "\n ERROR! t and %s conflict !\n\n", t0name);
    StdFace_exit(-1);
  }
  else if (creal(*t0) < 9999.0)
    fprintf(stdout, "  %15s = %-10.5f\n", t0name, creal(*t0));
  else if (creal(StdI->t) < 9999.0) {
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

  fprintf(fp, "%25.15e %25.15e\n", StdI->direct[0][0], StdI->direct[0][1]);
  fprintf(fp, "%25.15e %25.15e\n", StdI->direct[1][0], StdI->direct[1][1]);
  fprintf(fp, "%d %d\n", StdI->box[0][0], StdI->box[0][1]);
  fprintf(fp, "%d %d\n", StdI->box[1][0], StdI->box[1][1]);

  for (iCell = 0; iCell < StdI->NCell; iCell++) {
    for (isite = 0; isite < StdI->NsiteUC; isite++) {
      fprintf(fp, "%25.15e %25.15e\n",
        StdI->tau[isite][0] + (double)StdI->Cell[iCell][0],
        StdI->tau[isite][1] + (double)StdI->Cell[iCell][1]);
    }/*for (isite = 0; isite < StdI->NsiteUC; isite++)*/
  }/* for (iCell = 0; iCell < StdI->NCell; iCell++)*/
  if (strcmp(StdI->model, "kondo") == 0) {
    for (iCell = 0; iCell < StdI->NCell; iCell++) {
      for (isite = 0; isite < StdI->NsiteUC; isite++) {
        fprintf(fp, "%25.15e %25.15e\n",
          StdI->tau[isite][0] + (double)StdI->Cell[iCell][0],
          StdI->tau[isite][1] + (double)StdI->Cell[iCell][1]);
      }/*for (isite = 0; isite < StdI->NsiteUC; isite++)*/
    }/* for (iCell = 0; iCell < StdI->NCell; iCell++)*/
  }
  fclose(fp);

}/*void StdFace_PrintGeometry()*/
/*
 * Malloc Arrays for interactions
 */
void StdFace_MallocInteractions(struct StdIntList *StdI) {
  int kintr;
  /*
  Coulomb intra
  */
  StdI->CintraIndx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->Cintra = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++) {
    StdI->CintraIndx[kintr] = (int *)malloc(sizeof(int) * 1);
  }
  /*
  Coulomb inter
  */
  StdI->CinterIndx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->Cinter = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++) {
    StdI->CinterIndx[kintr] = (int *)malloc(sizeof(int) * 2);
  }
  /*
  Hund
  */
  StdI->HundIndx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->Hund = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++) {
    StdI->HundIndx[kintr] = (int *)malloc(sizeof(int) * 2);
  }
  /*
  Excahnge
  */
  StdI->ExIndx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->Ex = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++) {
    StdI->ExIndx[kintr] = (int *)malloc(sizeof(int) * 2);
  }
  /*
  PairLift
  */
  StdI->PLIndx = (int **)malloc(sizeof(int*) * StdI->nintr);
  StdI->PairLift = (double *)malloc(sizeof(double) * StdI->nintr);
  for (kintr = 0; kintr < StdI->nintr; kintr++) {
    StdI->PLIndx[kintr] = (int *)malloc(sizeof(int) * 2);
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
void StdFace_FoldSite2Dsub(struct StdIntList *StdI,
  int iW, int iL, int *iCell0, int *iCell1, int *iWfold, int *iLfold)
{
  int x0, x1, xW, xL;
  /*
  Transform to fractional coordinate (times NCell)
  */
  x0 = StdI->bW0sub * iW + StdI->bL0sub * iL;
  x1 = StdI->bW1sub * iW + StdI->bL1sub * iL;
  /*
  Which supercell contains this cell
  */
  *iCell0 = (x0 + StdI->NCellsub * 1000) / StdI->NCellsub - 1000;
  *iCell1 = (x1 + StdI->NCellsub * 1000) / StdI->NCellsub - 1000;
  /*
  Fractional coordinate (times NCell) in the original supercell
  */
  x0 -= StdI->NCellsub*(*iCell0);
  x1 -= StdI->NCellsub*(*iCell1);
  /**/
  xW = StdI->a0Wsub * x0 + StdI->a1Wsub * x1;
  xL = StdI->a0Lsub * x0 + StdI->a1Lsub * x1;
  *iWfold = (xW + StdI->NCellsub * 1000) / StdI->NCellsub - 1000;
  *iLfold = (xL + StdI->NCellsub * 1000) / StdI->NCellsub - 1000;
}
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
  int iWfold, iLfold, jWfold, jLfold, iW, iL;
  int NotUse1, NotUse2;
  int UnitNum0, UnitNum1;
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
    iW = StdI->Cell[iCell][0];
    iL = StdI->Cell[iCell][1];

    StdFace_FoldSite2Dsub(StdI, iW, iL, &NotUse1, &NotUse2, &iWfold, &iLfold);

    StdFace_FoldSite2D(StdI, iWfold, iLfold, &NotUse1, &NotUse2, &iWfold, &iLfold);

    if (iWfold == iW && iLfold == iL) {
      /*
      Translation operator in sub lattice
      */
      for (jCell = 0; jCell < StdI->NCell; jCell++) {

        jWfold = StdI->Cell[jCell][0] + iW;
        jLfold = StdI->Cell[jCell][1] + iL;
        StdFace_FoldSite2D(StdI, jWfold, jLfold, &UnitNum0, &UnitNum1, &jWfold, &jLfold);

        for (kCell = 0; kCell < StdI->NCell; kCell++) {
          if (jWfold == StdI->Cell[kCell][0] && jLfold == StdI->Cell[kCell][1]) {

            for (jsite = 0; jsite < StdI->NsiteUC; jsite++) {

              Sym[StdI->NSym][jCell*StdI->NsiteUC + jsite] = kCell*StdI->NsiteUC + jsite;
              Anti[StdI->NSym][jCell*StdI->NsiteUC + jsite] =
                StdI->AntiPeriod0 * UnitNum0 + StdI->AntiPeriod1 * UnitNum1;

              if (strcmp(StdI->model, "kondo") == 0) {
                Sym[StdI->NSym][StdI->nsite / 2 + jCell*StdI->NsiteUC + jsite] = StdI->nsite / 2 + kCell*StdI->NsiteUC + jsite;
                Anti[StdI->NSym][StdI->nsite / 2 + jCell*StdI->NsiteUC + jsite]
                  = StdI->AntiPeriod0 * UnitNum0 + StdI->AntiPeriod1 * UnitNum1;
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
      if (StdI->AntiPeriod0 == 1 || StdI->AntiPeriod1 == 1) {
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
}
/**
*
* Initialize sub Cell
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_InitSite2DSub(struct StdIntList *StdI)
{
  int ii, prod[4];
  /*
  check Input parameters
  */
  if ((StdI->Lsub != 9999 || StdI->Wsub != 9999)
    && (StdI->a0Lsub != 9999 || StdI->a0Wsub != 9999 || StdI->a1Lsub != 9999 || StdI->a1Wsub != 9999)) {
    fprintf(stdout, "\nERROR ! (Lsub, Wsub) and (a0Wsub, a0Lsub, a1Wsub, a1Lsub) conflict !\n\n");
    StdFace_exit(-1);
  }
  else if (StdI->Lsub != 9999 || StdI->Wsub != 9999) {
    if (StdI->Lsub == 9999) {
      fprintf(stdout, "\nERROR ! Wsub is specified, but Lsub is NOT specified !\n\n");
      StdFace_exit(-1);
    }
    else if (StdI->Wsub == 9999) {
      fprintf(stdout, "\nERROR ! Lsub is specified, but Wsub is NOT specified !\n\n");
      StdFace_exit(-1);
    }
    StdFace_PrintVal_i("Lsub", &StdI->Lsub, 1);
    StdFace_PrintVal_i("Wsub", &StdI->Wsub, 1);
    StdI->a0Wsub = StdI->Wsub;
    StdI->a0Lsub = 0;
    StdI->a1Wsub = 0;
    StdI->a1Lsub = StdI->Lsub;
  }
  else if (StdI->a0Lsub != 9999 || StdI->a0Wsub != 9999 || StdI->a1Lsub != 9999 || StdI->a1Wsub != 9999) {
    if (StdI->a0Wsub == 9999) {
      fprintf(stdout, "\nERROR ! a0Wsub is NOT specified !\n\n");
      StdFace_exit(-1);
    }
    else if (StdI->box[0][1] == 9999) {
      fprintf(stdout, "\nERROR ! a0Lsub is NOT specified !\n\n");
      StdFace_exit(-1);
    }
    else if (StdI->box[1][0] == 9999) {
      fprintf(stdout, "\nERROR ! a1Wsub is NOT specified !\n\n");
      StdFace_exit(-1);
    }
    else if (StdI->box[1][1] == 9999) {
      fprintf(stdout, "\nERROR ! a1Lsub is NOT specified !\n\n");
      StdFace_exit(-1);
    }
    StdFace_PrintVal_i("a0Wsub", &StdI->a0Wsub, 1);
    StdFace_PrintVal_i("a0Lsub", &StdI->a0Lsub, 0);
    StdFace_PrintVal_i("a1Wsub", &StdI->a1Wsub, 0);
    StdFace_PrintVal_i("a1Lsub", &StdI->a1Lsub, 1);
  }
  else {
    StdFace_PrintVal_i("a0Wsub", &StdI->a0Wsub, StdI->box[0][0]);
    StdFace_PrintVal_i("a0Lsub", &StdI->a0Lsub, StdI->box[0][1]);
    StdFace_PrintVal_i("a1Wsub", &StdI->a1Wsub, StdI->box[1][0]);
    StdFace_PrintVal_i("a1Lsub", &StdI->a1Lsub, StdI->box[1][1]);
  }
  /*
  Calculate reciplocal sublattice vectors (times NCellsub)
  */
  StdI->NCellsub = StdI->a0Wsub * StdI->a1Lsub - StdI->a0Lsub * StdI->a1Wsub;
  printf("         Number of Cell in sublattice: %d\n", StdI->NCellsub);
  if (StdI->NCell == 0) {
    StdFace_exit(-1);
  }

  StdI->bW0sub = StdI->a1Lsub;
  StdI->bL0sub = -StdI->a1Wsub;
  StdI->bW1sub = -StdI->a0Lsub;
  StdI->bL1sub = StdI->a0Wsub;
  if (StdI->NCellsub < 0) {
    StdI->bW0sub *= -1;
    StdI->bL0sub *= -1;
    StdI->bW1sub *= -1;
    StdI->bL1sub *= -1;
    StdI->NCellsub *= -1;
  }/*if (StdI->NCell < 0)*/
  /*
   Check : Is the sublattice commensurate ?
  */
  prod[0] = StdI->bW0sub * (double)StdI->box[0][0] + StdI->bL0sub * (double)StdI->box[0][1];
  prod[1] = StdI->bW1sub * (double)StdI->box[0][0] + StdI->bL1sub * (double)StdI->box[0][1];
  prod[2] = StdI->bW0sub * (double)StdI->box[1][0] + StdI->bL0sub * (double)StdI->box[1][1];
  prod[3] = StdI->bW1sub * (double)StdI->box[1][0] + StdI->bL1sub * (double)StdI->box[1][1];

  for (ii = 0; ii < 4; ii++){
    if (prod[ii] % StdI->NCellsub != 0) {
      printf("\n ERROR ! Sublattice is INCOMMENSURATE !\n\n");
      StdFace_exit(-1);
    }/*if (fabs(prod[ii]) > 0.00001)*/
  }/*for (ii = 0; ii < 4; ii++)*/

}/*void StdFace_InitSite2DSub*/
/*
 *Generate orbitalindex
*/
void StdFace_generate_orb(struct StdIntList *StdI) {
  int iCell, jCell, kCell, iW, iL, jW, jL, iCell2, jCell2;
  int NotUse1, NotUse2, iWfold, iLfold, jWfold, jLfold, iOrb;
  int isite, jsite;
  int dW, dL, dWfold, dLfold, UnitNum0, UnitNum1, Anti;
  int **CellDone;

  StdFace_InitSite2DSub(StdI);

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
    /**/
    iW = StdI->Cell[iCell][0];
    iL = StdI->Cell[iCell][1];

    StdFace_FoldSite2Dsub(StdI, iW, iL, &NotUse1, &NotUse2, &iWfold, &iLfold);

    StdFace_FoldSite2D(StdI, iWfold, iLfold, &NotUse1, &NotUse2, &iWfold, &iLfold);

    for (kCell = 0; kCell < StdI->NCell; kCell++) {
      if (iWfold == StdI->Cell[kCell][0] && iLfold == StdI->Cell[kCell][1]) {
        iCell2 = kCell;
      }
    }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/

    for (jCell = 0; jCell < StdI->NCell; jCell++) {
      /**/
      jW = StdI->Cell[jCell][0];
      jL = StdI->Cell[jCell][1];

      jW = jW + iWfold - iW;
      jL = jL + iLfold - iL;

      StdFace_FoldSite2D(StdI, jW, jL, &NotUse1, &NotUse2, &jWfold, &jLfold);

      for (kCell = 0; kCell < StdI->NCell; kCell++) {
        if (jWfold == StdI->Cell[kCell][0] && jLfold == StdI->Cell[kCell][1]) {
          jCell2 = kCell;
        }
      }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/
      /*
       AntiPeriodic factor
      */
      dW = StdI->Cell[jCell][0] - StdI->Cell[iCell][0];
      dL = StdI->Cell[jCell][1] - StdI->Cell[iCell][1];
      StdFace_FoldSite2D(StdI, dW, dL, &UnitNum0, &UnitNum1, &dWfold, &dLfold);
      Anti = StdI->AntiPeriod0 * UnitNum0 + StdI->AntiPeriod1 * UnitNum1;
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
