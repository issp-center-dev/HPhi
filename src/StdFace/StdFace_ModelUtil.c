/*
HPhi  -  Quantum Lattice Model Simulator
Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima

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

/**
 *
 * Add transfer to the list
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_trans(
struct StdIntList *StdI,
  _Dcomplex trans0 /**< [in] Hopping integral t, mu, etc. */,
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
  _Dcomplex trans0 /**< [in] Hopping integral t, mu, etc. */,
  int isite /**< [in] i for c_{i sigma}^dagger*/,
  int jsite /**< [in] j for c_{j sigma'}*/,
  char *mode
  )
{
  int ispin;

  for (ispin = 0; ispin < 2; ispin++) {
    StdFace_trans(StdI, trans0, isite, ispin, jsite, ispin);
    if(mode = "hopp") 
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
  _Dcomplex intr0 /**< [in] Interaction U, V, J, etc.*/,
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
  double *J,
  int Si2 /**< [in] Spin moment in i site*/,
  int Sj2 /**< [in] Spin moment in j site*/,
  int isite /**< [in] i of S_i */,
  int jsite /**< [in] j of S_j */)
{
  int ispin, jspin;
  double Si, Sj, Siz, Sjz;
  _Dcomplex intr0;

  Si = 0.5 * (double)Si2;
  Sj = 0.5 * (double)Sj2;

  for (ispin = 0; ispin <= Si2; ispin++) {
    Siz = (double)ispin - Si;
    for (jspin = 0; jspin <= Sj2; jspin++) {
      Sjz = (double)jspin - Sj;
      /*
       S_{i z} * S_{j z}
      */
      intr0 = J[2][2] * Siz * Sjz;
      StdFace_intr(StdI, intr0,
        isite, ispin, isite, ispin, jsite, jspin, jsite, jspin);
      /*
       S_i^+ S_j^- + S_j^+ S_i^-
      */
      if (ispin < Si2 && jspin < Sj2) {
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
      if (ispin < Si2 && jspin < Sj2) {
        intr0 = 0.5 * 0.5 * (J[0][0] - J[1][1] + I*(J[0][1] + J[0][1]))
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
        intr0 = (J[0][2] - I * J[1][2]) * sqrt(Si * (Si + 1.0) - Siz * (Siz + 1.0)) * Sjz;
        StdFace_intr(StdI, intr0,
          isite, ispin + 1, isite, ispin, jsite, jspin, jsite, jspin);
        StdFace_intr(StdI, conj(intr0),
          jsite, jspin, jsite, jspin, isite, ispin, isite, ispin + 1);
      }/*if (ispin < Si2)*/
      /*
       S_{i z} S_j^+ + S_j^- S_{i z}
      */
      if (jspin < Sj2) {
        intr0 = (J[0][2] - I * J[1][2]) * Siz * sqrt(Sj * (Sj + 1.0) - Sjz * (Sjz + 1.0));
        StdFace_intr(StdI, intr0,
          isite, ispin, isite, ispin, jsite, jspin + 1, jsite, jspin);
        StdFace_intr(StdI, conj(intr0),
          jsite, jspin, jsite, jspin + 1, isite, ispin, isite, ispin);
      }/*if (jspin < Sj2)*/

    }/*for (jspin = 0; jspin <= Sj2; jspin++)*/
  }/*for (ispin = 0; ispin <= Si2; ispin++)*/
}

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
  int ispin, jspin;
  for (ispin = 0; ispin < 2; ispin++){
    for (jspin = 0; jspin < 2; jspin++){
      StdFace_intr(StdI, V, isite, ispin, isite, ispin, jsite, jspin, jsite, jspin);
    }
  }
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
* Print a valiable (complex) read from the input file
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_PrintVal_c(
  char* valname /**< [in] Name of the valiable*/,
  _Dcomplex *val /**< [inout] Valiable to be set*/,
  _Dcomplex val0 /**< [in] The default value*/)
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
    fprintf(stderr, "\n Check !  %s is SPECIFIED but will NOT be USED. \n", valname);
    fprintf(stderr, "            Please COMMENT-OUT this line \n");
    fprintf(stderr, "            or check this input is REALLY APPROPRIATE for your purpose ! \n\n");
    exit(-1);
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
  _Dcomplex val /**< [in]*/)
{
  if (creal(val) < 9999.0) {
    fprintf(stderr, "\n Check !  %s is SPECIFIED but will NOT be USED. \n", valname);
    fprintf(stderr, "            Please COMMENT-OUT this line \n");
    fprintf(stderr, "            or check this input is REALLY APPROPRIATE for your purpose ! \n\n");
    exit(-1);
  }
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
    fprintf(stderr, "\n Check !  %s is SPECIFIED but will NOT be USED. \n", valname);
    fprintf(stderr, "            Please COMMENT-OUT this line \n");
    fprintf(stderr, "            or check this input is REALLY APPROPRIATE for your purpose ! \n\n");
    exit(-1);
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
    fprintf(stderr, "ERROR ! %s is NOT specified !\n", valname);
    exit(-1);
  }
  else fprintf(stdout, "  %15s = %-3d\n", valname, val);
}

/**
*
* Define whether the specified site is in the unit cell or not.
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_FoldSite2D(struct StdIntList *StdI, 
  int iW, int iL, int *iCell0, int *iCell1, int *iWfold, int *iLfold)
{
  double x0, x1;

  /*
   Transform to fractional coordinate
  */
  x0 = StdI->bW0 * (double)iW + StdI->bL0 * (double)iL;
  x1 = StdI->bW1 * (double)iW + StdI->bL1 * (double)iL;
  /**/
  *iCell0 = (int)floor(x0 + 1.0e-8);
  *iCell1 = (int)floor(x1 + 1.0e-8);
  /**/
  x0 -= (double)*iCell0;
  x1 -= (double)*iCell1;
  /**/
  *iWfold = (int)((double)StdI->a0W * x0 + (double)StdI->a1W * x1 + 1.0e-8);
  *iLfold = (int)((double)StdI->a0L * x0 + (double)StdI->a1L * x1 + 1.0e-8);
}

/**
*
* Initialize Cell
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_InitSite2D(struct StdIntList *StdI, FILE *fp)
{
  int Wmin, Wmax, Lmin, Lmax;
  int iW, iL, ipos;
  int iCell, iCell0, iCell1, iWfold, iLfold;
  double pos[4][2], xmin, xmax, det;

  if (StdI->a0W * StdI->a1L - StdI->a0L * StdI->a1W == 0) {
    exit(-1);
  }

  det = (double)(StdI->a0W * StdI->a1L - StdI->a0L * StdI->a1W);
  StdI->bW0 = (double)StdI->a1L / det;
  StdI->bL0 = - (double)StdI->a1W / det;
  StdI->bW1 = - (double)StdI->a0L / det;
  StdI->bL1 = (double)StdI->a0W / det;

  Wmax = 0;
  if (StdI->a0W > Wmax) Wmax = StdI->a0W;
  if (StdI->a1W > Wmax) Wmax = StdI->a1W;
  if (StdI->a0W + StdI->a1W > Wmax) Wmax = StdI->a0W + StdI->a1W;
  /**/
  Wmin = 0;
  if (StdI->a0W < Wmin) Wmin = StdI->a0W;
  if (StdI->a1W < Wmin) Wmin = StdI->a1W;
  if (StdI->a0W + StdI->a1W < Wmin) Wmin = StdI->a0W + StdI->a1W;
  /**/
  Lmax = 0;
  if (StdI->a0L > Lmax) Lmax = StdI->a0L;
  if (StdI->a1L > Lmax) Lmax = StdI->a1L;
  if (StdI->a0L + StdI->a1L > Lmax) Lmax = StdI->a0L + StdI->a1L;
  /**/
  Lmin = 0;
  if (StdI->a0L < Lmin) Lmin = StdI->a0L;
  if (StdI->a1L < Lmin) Lmin = StdI->a1L;
  if (StdI->a0L + StdI->a1L < Lmin) Lmin = StdI->a0L + StdI->a1L;

  StdI->Cell = (int **)malloc(sizeof(int*) * (Wmax - Wmin + 1) * (Lmax - Lmin + 1));
  for (iCell = 0; iCell < (Wmax - Wmin + 1) * (Lmax - Lmin + 1); iCell++) {
    StdI->Cell[iCell] = (int *)malloc(sizeof(int) * 2);
  }/*for (iCell = 0; iCell < (Wmax - Wmin + 1) * (Lmax - Lmin + 1); iCell++)*/

  StdI->NCell = 0;
  for (iL = Lmin; iL <= Lmax; iL++) {
    for (iW = Wmin; iW <= Wmax; iW++) {
      StdFace_FoldSite2D(StdI, iW, iL, &iCell0, &iCell1, &iWfold, &iLfold);
      if (iCell0 == 0 && iCell1 == 0) {
        StdI->Cell[StdI->NCell][0] = iW;
        StdI->Cell[StdI->NCell][1] = iL;
        StdI->NCell += 1;
      }/*if (lUC == 1)*/
    }/*for (iW = Wmin; iW <= Wmax; iW++*/
  }/*for (iL = Lmin; iL <= Lmax; iL++)*/
  /*
   Initialize gnuplot
  */
  pos[0][0] = 0.0;
  pos[0][1] = 0.0;
  pos[1][0] = StdI->Lx;
  pos[1][1] = StdI->Ly;
  pos[2][0] = StdI->Wx;
  pos[2][1] = StdI->Wy;
  pos[3][0] = StdI->Wx + StdI->Lx;
  pos[3][1] = StdI->Wy + StdI->Ly;
  /**/
  xmin = 0.0;
  xmax = 0.0;
  for (ipos = 0; ipos < 4; ipos++) {
    if (pos[ipos][0] < xmin) xmin = pos[ipos][0];
    if (pos[ipos][0] > xmax) xmax = pos[ipos][0];
    if (pos[ipos][1] < xmin) xmin = pos[ipos][1];
    if (pos[ipos][1] > xmax) xmax = pos[ipos][1];
  }
  xmin -= 5.0;
  xmax += 5.0;

  fprintf(fp, "set xrange [%f: %f]\n", xmin, xmax);
  fprintf(fp, "set yrange [%f: %f]\n", xmin, xmax);
  fprintf(fp, "set size square\n", xmin, xmax);

  fprintf(fp, "set style line 1 lc 1 lt 1\n");
  fprintf(fp, "set style line 2 lc 3 lt 0\n");

  fprintf(fp, "set arrow from %f, %f to %f, %f nohead\n", pos[0][0], pos[0][1], pos[1][0], pos[1][1]);
  fprintf(fp, "set arrow from %f, %f to %f, %f nohead\n", pos[1][0], pos[1][1], pos[2][0], pos[2][1]);
  fprintf(fp, "set arrow from %f, %f to %f, %f nohead\n", pos[2][0], pos[2][1], pos[3][0], pos[3][1]);
  fprintf(fp, "set arrow from %f, %f to %f, %f nohead\n", pos[3][0], pos[3][1], pos[0][0], pos[0][1]);

}

void StdFace_SetLabel(struct StdIntList *StdI, FILE *fp, 
  int iW, int iL, int jW, int jL, 
  double xiW, double xiL, double xjW, double xjL, 
  int isite, int *jsite, int connect, int NsiteUC, int isiteUC, char *model)
{
  int iCell, jCell, jCell0, jCell1, jWfold, jLfold;
  double xi, yi, xj, yj;

  StdFace_FoldSite2D(StdI, iW, iL, &jCell0, &jCell1, &jWfold, &jLfold);
  /**/
  for (iCell = 0; iCell < StdI->NCell; iCell++) {
    if (jWfold == StdI->Cell[iCell][0] && jLfold == StdI->Cell[iCell][1]) {
      jCell = iCell;
      break;
    }
  }/*for (iCell = 0; iCell < StdI->NCell; iCell++)*/
  jsite = jCell * NsiteUC + isiteUC;
  if (model == "kondo") jsite += StdI->NCell * NsiteUC;

  xi = StdI->Lx * xiL + StdI->Wx * xiW;
  yi = StdI->Ly * xiL + StdI->Wy * xiW;

  xj = StdI->Lx * xjL + StdI->Wx * xjW;
  yj = StdI->Ly * xjL + StdI->Wy * xjW;

  fprintf(fp, "set label \"%2d\" at %f, %f center front\n", isite, xi, yi);
  fprintf(fp, "set label \"%2d\" at %f, %f center front\n", jsite, xj, yj);
  fprintf(fp, "set arrow from %f, %f to %f, %f nohead ls %d\n", xi, yi, xj, yj, connect);

}
