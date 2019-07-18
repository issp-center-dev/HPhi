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

#include "Common.h"
#include "HPhiTrans.h"
#include "FileIO.h"
#include "wrapperMPI.h"

/**
 * @file   HPhiTrans.c
 *
 * @brief  Check the inputted transfer integrals.
 *
 * @version 0.1
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */

/** 
 *
 * @brief  Function of checking transfers not to count the same type of operators.\n
 * @note   The same type transfer integrals such as @f$ t^{(1)}_{ij}c_i a_j +  t^{(2)}_{ij}c_i a_j @f$ should be summarized as @f$ t_{ij} c_i a_j @f$ before calling this function.
 *
 * @param X [in] Struct to get the information of the operators of transfer integrals.
 * @retval 0  normally finished
 * @retval -1 unnormally finished
 * 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 *
 */
int HPhiTrans(struct BindStruct *X) {
  FILE *fp_err;
  char sdt_err[D_FileNameMax];

  int i, k;
  int cnt_trans;

  strcpy(sdt_err, cFileNameWarningOnTransfer);
  if (childfopenMPI(sdt_err, "w", &fp_err) != 0) {
    return -1;
  }
  fclose(fp_err);

  //Transefer
  cnt_trans = 0;

  for (i = 0; i < X->Def.EDNTransfer; i++) {
    // eliminate double counting
    for (k = 0; k < cnt_trans; k++) {
      if (X->Def.EDGeneralTransfer[i][1] == X->Def.EDGeneralTransfer[k][1]
          && X->Def.EDGeneralTransfer[i][3] == X->Def.EDGeneralTransfer[k][3]) {
        if (X->Def.EDGeneralTransfer[i][0] == X->Def.EDGeneralTransfer[k][0]
            && X->Def.EDGeneralTransfer[i][2] == X->Def.EDGeneralTransfer[k][2]) {
          sprintf(sdt_err, "%s", cErrTransfer);
          childfopenMPI(sdt_err, "a", &fp_err);
          fprintf(fp_err, cErrDoubleCounting, X->Def.EDGeneralTransfer[k][0], X->Def.EDGeneralTransfer[k][2],
                  X->Def.EDGeneralTransfer[k][1], X->Def.EDGeneralTransfer[k][3]);
          fclose(fp_err);
        }
      }
    }
    cnt_trans += 1;
  }

  return 0;
}

/**
 * @brief  Function of getting transfer with peierls
 * @note Not used now and should be delete in ver.2.1.
 * @param X data list for calculation
 * @param time time
 * @retval 0  normally finished
 * @retval -1 unnormally finished
 *
 * @author Kota Ido (The University of Tokyo)
 */
int TransferWithPeierls(struct BindStruct *X, const double time) {
  int i;
  int ri_x, rj_x;
  int ri_y, rj_y;
  double complex dir;
  const int Mode = (int) (X->Def.ParaLaser[0]);
  const double Avp = X->Def.ParaLaser[1];
  const double omega = X->Def.ParaLaser[2];
  const double time_d = X->Def.ParaLaser[3];
  const double time_c = X->Def.ParaLaser[4];
  const int Lx = (int) (X->Def.ParaLaser[5]);
  const int Ly = (int) (X->Def.ParaLaser[6]);
  const double dirX = X->Def.ParaLaser[7];
  const double dirY = X->Def.ParaLaser[8];
  const double dt = time - time_c;
  const double dt2 = time - (time_c + time_d);
  const double td = time_c / 3.0;
  double VecPot = 0.0;

  if (Mode == 0) {//Gaussian Wave
    VecPot = Avp * cos(omega * dt) * exp(-dt * dt / (2.0 * time_d * time_d));
  } else if (Mode == 1) {//Cosine Wave
    VecPot = Avp * sin(omega * dt);
  } else if (Mode == 2) {//DC Limit
    VecPot = Avp * dt;
  } else if (Mode == 3) {//Pulse
    VecPot = Avp * exp(-time_d * dt);
  } else if (Mode == 4) {//Linear
    if (dt <= 0.0) {
      VecPot = 0.0;
    } else if (dt < time_d) {
      VecPot = Avp * cos(omega * dt) * (dt / time_d);
    } else {
      VecPot = Avp * cos(omega * dt);
    }
  } else if (Mode == 5) {//Linear
    if (time <= 0.0) {
      VecPot = 0.0;
    } else if (time < time_c) {
      VecPot = Avp * cos(omega * dt) * exp(-dt * dt / (2.0 * td * td));
    } else if (time < time_c + time_d) {
      VecPot = Avp * cos(omega * dt);
    } else {
      VecPot = Avp * cos(omega * dt) * exp(-dt2 * dt2 / (2.0 * td * td));
    }
  }

  for (i = 0; i < X->Def.EDNTransfer; i++) {
    ri_x = X->Def.EDGeneralTransfer[i][0] % Lx;
    rj_x = X->Def.EDGeneralTransfer[i][2] % Lx;
    ri_y = X->Def.EDGeneralTransfer[i][0] / Lx;
    rj_y = X->Def.EDGeneralTransfer[i][2] / Lx;
    if (ri_x - rj_x > 1) {
      rj_x += Lx;
    } else if (ri_x - rj_x < -1) {
      rj_x -= Lx;
    }
    if (ri_y - rj_y > 1) {
      rj_y += Ly;
    } else if (ri_y - rj_y < -1) {
      rj_y -= Ly;
    }
    dir = dirX * (ri_x - rj_x) + dirY * (ri_y - rj_y);

    X->Def.EDParaGeneralTransfer[i] = X->Def.ParaGeneralTransfer[i] * cexp(-I * VecPot * dir);
  }

  return 0;
}

/**
 * @brief  Function of getting transfer for quench
 * @note Not used now and should be delete in ver.2.1.
 * @param X data list for calculation
 * @param time time
 * @retval 0  normally finished
 * @retval -1 unnormally finished
 *
 * @author Kota Ido (The University of Tokyo)
 */
int TransferForQuench(struct BindStruct *X, const double time) {
  int i;
  const int Mode = (int) (X->Def.ParaLaser[0]);
  const double Avp = X->Def.ParaLaser[1];
  const double time_d = X->Def.ParaLaser[3];
  const double time_c = X->Def.ParaLaser[4];
  const double dt = time - time_c;
  double Bessel = 0.0;

  if (Mode == 0) {//Gaussian Wave
    if (dt <= 0.0) {
      Bessel = 0.0;
    } else if (dt < time_d) {
      Bessel = j0(Avp * dt / time_d);
    } else {
      Bessel = j0(Avp);
    }
  }

  for (i = 0; i < X->Def.EDNTransfer; i++) {
    X->Def.EDParaGeneralTransfer[i] = X->Def.ParaGeneralTransfer[i] * Bessel;
  }

  return 0;
}

