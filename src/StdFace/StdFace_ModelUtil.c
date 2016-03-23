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
#include "StdFace_vals.h"

/**
 *
 * Add transfer to the list
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_trans(
struct StdIntList *StdI,
  double trans0 /**< [in] Hopping integral t, mu, etc. */,
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
* Add Longitudinal magnetic field to the list
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_MagLong(
struct StdIntList *StdI,
  int S2 /**< [in] Spin moment in i site*/,
  double Mag /**< [in] Longitudinal magnetic field h. */,
  int isite /**< [in] i for c_{i sigma}^dagger*/)
{
  int ispin;
  double S, Sz;
 
  S = (double)S2 * 0.5;

  for (ispin = 0; ispin <= S2; ispin++){
    Sz = (double)ispin - S;
    StdFace_trans(StdI, Mag * Sz, isite, ispin, isite, ispin);
  }
}

/**
*
* Add Transvars magnetic field to the list
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_MagTrans(
struct StdIntList *StdI,
  int S2 /**< [in] Spin moment in i site*/,
  double Mag /**< [in] Transvars magnetic field h. */,
  int isite /**< [in] i for c_{i sigma}^dagger*/)
{
  int ispin;
  double S, Sz;

  S = (double)S2 * 0.5;

  for (ispin = 0; ispin < S2; ispin++){
    Sz = (double)ispin - S;
    StdFace_trans(StdI, 0.5 * Mag * sqrt(S*(S + 1.0) - Sz*(Sz + 1.0)), isite, ispin + 1, isite, ispin);
    StdFace_trans(StdI, 0.5 * Mag * sqrt(S*(S + 1.0) - Sz*(Sz + 1.0)), isite, ispin, isite, ispin + 1);
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
  double intr0 /**< [in] Interaction U, V, J, etc.*/,
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
 * Add Exchange term (S_i^+ S_j^- + S_i^- S_j^+) to the list
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_exchange(
struct StdIntList *StdI,
  int Si2 /**< [in] Spin moment in i site*/,
  int Sj2 /**< [in] Spin moment in j site*/,
  double Jexc /**< [in] Interaction J_{x y}, etc.*/,
  int isite /**< [in] i of S_i */, 
  int jsite /**< [in] j of S_j */)
{
  int ispin, jspin;
  double intr0, Si, Sj, Siz, Sjz;

  Si = 0.5 * (double)Si2;
  Sj = 0.5 * (double)Sj2;

  for (ispin = 0; ispin < Si2; ispin++){
    Siz = (double)ispin - Si;
    for (jspin = 0; jspin < Sj2; jspin++){
      Sjz = (double)jspin - Sj;

      intr0 = 0.5 * Jexc * sqrt(Si * (Si + 1.0) - Siz * (Siz + 1.0)) 
                         * sqrt(Sj * (Sj + 1.0) - Sjz * (Sjz + 1.0));

      StdFace_intr(StdI, intr0,
        isite, ispin + 1, isite, ispin, jsite, jspin, jsite, jspin + 1);
      StdFace_intr(StdI, intr0,
        isite, ispin, isite, ispin + 1, jsite, jspin + 1, jsite, jspin);
    }
  }
}

/**
 *
 * Add S^+ S^+ + S^-S^- term to the list
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_Kitaev(
struct StdIntList *StdI,
  int Si2 /**< [in] Spin moment in i site*/,
  int Sj2 /**< [in] Spin moment in j site*/,
  double Jflip /**< [in] Interaction J_x - J_y, etc.*/,
  int isite /**< [in] i of S_i */,
  int jsite /**< [in] j of S_j */)
{
  int ispin, jspin;
  double intr0, Si, Sj, Siz, Sjz;

  Si = 0.5 * (double)Si2;
  Sj = 0.5 * (double)Sj2;

  for (ispin = 0; ispin < Si2; ispin++){
    Siz = (double)ispin - Si;
    for (jspin = 0; jspin < Sj2; jspin++){
      Sjz = (double)jspin - Sj;

      intr0 = 0.5 * Jflip * sqrt(Si * (Si + 1.0) - Siz * (Siz + 1.0))
                          * sqrt(Sj * (Sj + 1.0) - Sjz * (Sjz + 1.0));

      StdFace_intr(StdI, intr0,
        isite, ispin + 1, isite, ispin, jsite, jspin + 1, jsite, jspin);
      StdFace_intr(StdI, intr0,
        isite, ispin, isite, ispin + 1, jsite, jspin, jsite, jspin + 1);
    }
  }
}

/**
 *
 * Add S_z S_z term to the list
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_SzSz(
struct StdIntList *StdI,
  int Si2 /**< [in] Spin moment in i site*/,
  int Sj2 /**< [in] Spin moment in j site*/,
  double Jexc /**< [in] Interaction J_z, etc.*/,
  int isite /**< [in] i of S_i */,
  int jsite /**< [in] j of S_j */)
{
  int ispin, jspin;
  double intr0;

  for (ispin = 0; ispin <= Si2; ispin++){
    for (jspin = 0; jspin <= Sj2; jspin++){
      intr0 = Jexc * ((double)ispin - (double)Si2 * 0.5)
                   * ((double)jspin - (double)Sj2 * 0.5);
      StdFace_intr(StdI, intr0,
        isite, ispin, isite, ispin, jsite, jspin, jsite, jspin);
    }
  }
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
