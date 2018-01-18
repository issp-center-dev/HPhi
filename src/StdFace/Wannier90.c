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
@brief Standard mode for wannier90
*/
#include "StdFace_vals.h"
#include "StdFace_ModelUtil.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
/**
@brief Read Geometry file for wannier90
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static void geometry_W90(
  struct StdIntList *StdI//!<[inout]
)
{
  int isite, ii, ierr;
  char filename[256];
  FILE *fp;

  sprintf(filename, "%s_geom.dat", StdI->CDataFileHead);
  fprintf(stdout, "    Wannier90 Geometry file = %s\n", filename);

  fp = fopen(filename, "r");
  /**@brief
   Direct lattice vector StdIntList::direct
  */
  for (ii = 0; ii < 3; ii++)
    ierr = fscanf(fp, "%lf%lf%lf", &StdI->direct[ii][0], &StdI->direct[ii][1], &StdI->direct[ii][2]);
  if(ierr == EOF) printf("%d\n", ierr);
  /**@brief
   Intrinsic site position StdIntList::tau and its number StdIntList::NsiteUC
  */
  for (isite = 0; isite < StdI->NsiteUC; isite++) free(StdI->tau[isite]);
  free(StdI->tau);
  ierr = fscanf(fp, "%d", &StdI->NsiteUC);
  fprintf(stdout, "    Number of Correlated Sites = %d\n", StdI->NsiteUC);

  StdI->tau = (double **)malloc(sizeof(double*) * StdI->NsiteUC);
  for (ii = 0; ii < StdI->NsiteUC; ii++) StdI->tau[ii] = (double *)malloc(sizeof(double) * 3);

  for (isite = 0; isite < StdI->NsiteUC; isite++)
    ierr = fscanf(fp, "%lf%lf%lf", &StdI->tau[isite][0], &StdI->tau[isite][1], &StdI->tau[isite][2]);
  fclose(fp);

  printf("    Direct lattice vectors:\n");
  for (ii = 0; ii < 3; ii++) printf("      %10.5f %10.5f %10.5f\n",
    StdI->direct[ii][0], StdI->direct[ii][1], StdI->direct[ii][2]);
  printf("    Wannier centres:\n");
  for (isite = 0; isite < StdI->NsiteUC; isite++) printf("      %10.5f %10.5f %10.5f\n",
    StdI->tau[isite][0], StdI->tau[isite][1], StdI->tau[isite][2]);

}/*static void geometry_W90(struct StdIntList *StdI) */
/**
@brief Read Wannier90 hamiltonian file (*_hr) and compute the number of effective term
@author Mitsuaki Kawamura (The University of Tokyo)
*/
static int read_W90_query(
  struct StdIntList *StdI,//!<[inout]
  char *filename,//!<[in] Input file name
  double cutoff//!<[in] Threshold for the Hamiltonian 
) 
{
  FILE *fp;
  int nMat;
  int ierr, nWan, nWSC, iWSC, jWSC, iWan, jWan, iWan0, jWan0, ii;
  double dtmp[2];
  char ctmp[256], *ctmp2;
  double complex ***Mat_tot;
  int **indx_tot;

  fprintf(stdout, "   Wannier90 file = %s\n", filename);
  /*
  Header part
  */
  fp = fopen(filename, "r");
  ctmp2 = fgets(ctmp, 256, fp);
  ierr = fscanf(fp, "%d", &nWan);
  if(ierr == EOF) printf("%d %s\n", ierr, ctmp2);
  ierr = fscanf(fp, "%d", &nWSC);
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    ierr = fscanf(fp, "%d", &ii);
  }
  fprintf(stdout, "           Number of Wannier = %d\n", nWan);
  fprintf(stdout, " Number of Wigner-Seitz Cell = %d\n", nWSC);
  /*
   Allocation of matgrix element and its index
  */
  Mat_tot = (double complex ***)malloc(sizeof(double complex **) * nWSC);
  indx_tot = (int **)malloc(sizeof(int*) * nWSC);
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    Mat_tot[iWSC] = (double complex **)malloc(sizeof(double complex *) * nWan);
    indx_tot[iWSC] = (int *)malloc(sizeof(int) * 3);
    for (iWan = 0; iWan < nWan; iWan++) {
      Mat_tot[iWSC][iWan] = (double complex *)malloc(sizeof(double complex) * nWan);
    }
  }
  /*
  Read body
  */
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < nWan; iWan++) {
      for (jWan = 0; jWan < nWan; jWan++) {
        ierr = fscanf(fp, "%d%d%d%d%d%lf%lf",
          &indx_tot[iWSC][0], &indx_tot[iWSC][1], &indx_tot[iWSC][2], &iWan0, &jWan0,
          &dtmp[0], &dtmp[1]);
        if(iWan0 <= StdI->NsiteUC && jWan0 <= StdI->NsiteUC)
          Mat_tot[iWSC][iWan0 - 1][jWan0 - 1] = dtmp[0] + I * dtmp[1];
      }
    }
    /**@brief
    (1) Apply inversion symmetry and delete duplication
    */
    for (jWSC = 0; jWSC < iWSC; jWSC++) {
      if (
        indx_tot[iWSC][0] == -indx_tot[jWSC][0] &&
        indx_tot[iWSC][1] == -indx_tot[jWSC][1] &&
        indx_tot[iWSC][2] == -indx_tot[jWSC][2]
        )
        for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
          for (jWan = 0; jWan < StdI->NsiteUC; jWan++) {
            Mat_tot[iWSC][iWan][jWan] = 0.0;
          }
        }
    }/*for (jWSC = 0; jWSC < iWSC; jWSC++)*/
    if (indx_tot[iWSC][0] == 0 &&
        indx_tot[iWSC][1] == 0 &&
        indx_tot[iWSC][2] == 0) 
      for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
        for (jWan = 0; jWan < iWan; jWan++) {
          Mat_tot[iWSC][iWan][jWan] = 0.0;
        }
      }
  }/*for (iWSC = 0; iWSC < nWSC; iWSC++)*/
  fclose(fp);
  /**@brief
  (3-1)  Compute the number of terms lerger than cut-off.
  */
  fprintf(stdout, "\n      EFFECTIVE terms:\n");
  nMat = 0;
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
      for (jWan = 0; jWan < StdI->NsiteUC; jWan++) {
        if (cutoff < cabs(Mat_tot[iWSC][iWan][jWan])) {
          fprintf(stdout, "        %5d%5d%5d%5d%5d%12.6f%12.6f\n",
            indx_tot[iWSC][0], indx_tot[iWSC][1], indx_tot[iWSC][2], iWan, jWan,
            creal(Mat_tot[iWSC][iWan][jWan]), cimag(Mat_tot[iWSC][iWan][jWan]));
          nMat += 1;
        }
      }
    }
  }
  fprintf(stdout, "      Total number of EFFECTIVE term = %d\n", nMat);

  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < nWan; iWan++) {
      free(Mat_tot[iWSC][iWan]);
    }
    free(Mat_tot[iWSC]);
    free(indx_tot[iWSC]);
  }
  free(Mat_tot);
  free(indx_tot);

  return nMat;
}/*static int read_W90_query(struct StdIntList *StdI, char *model)*/
 /**
 @brief Read Wannier90 hamiltonian file (*_hr)
 @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void read_W90(
  struct StdIntList *StdI,//!<[inout]
  char *filename,//!<[in] Input file name
  double cutoff,//!<[in] Threshold for the Hamiltonian 
  double complex *Mat,//!<[out] Matrix element
  int **Matindx//!<[out] R, band index of matrix element
)
{
  FILE *fp;
  int nMat;
  int ierr, nWan, nWSC, iWSC, jWSC, iWan, jWan, iWan0, jWan0, ii;
  double dtmp[2];
  char ctmp[256], *ctmp2;
  double complex ***Mat_tot;
  int **indx_tot;
  /*
  Header part
  */
  fp = fopen(filename, "r");
  ctmp2 = fgets(ctmp, 256, fp);
  ierr = fscanf(fp, "%d", &nWan);
  if (ierr == EOF) printf("%d %s\n", ierr, ctmp2);
  ierr = fscanf(fp, "%d", &nWSC);
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    ierr = fscanf(fp, "%d", &ii);
  }
  /*
  Malloc Matrix elements and their indices
  */
  Mat_tot = (double complex ***)malloc(sizeof(double complex **) * nWSC);
  indx_tot = (int **)malloc(sizeof(int*) * nWSC);
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    Mat_tot[iWSC] = (double complex **)malloc(sizeof(double complex *) * nWan);
    indx_tot[iWSC] = (int *)malloc(sizeof(int) * 3);
    for (iWan = 0; iWan < nWan; iWan++) {
      Mat_tot[iWSC][iWan] = (double complex *)malloc(sizeof(double complex) * nWan);
    }
  }
  /*
  Read body
  */
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < nWan; iWan++) {
      for (jWan = 0; jWan < nWan; jWan++) {
        ierr = fscanf(fp, "%d%d%d%d%d%lf%lf",
          &indx_tot[iWSC][0], &indx_tot[iWSC][1], &indx_tot[iWSC][2],
          &iWan0, &jWan0,
          &dtmp[0], &dtmp[1]);
        if (iWan0 <= StdI->NsiteUC && jWan0 <= StdI->NsiteUC)
          Mat_tot[iWSC][iWan0 - 1][jWan0 - 1] = dtmp[0] + I * dtmp[1];
      }
    }
    /**@brief
    (1) Apply inversion symmetry and delete duplication
    */
    for (jWSC = 0; jWSC < iWSC; jWSC++) {
      if (
        indx_tot[iWSC][0] == -indx_tot[jWSC][0] &&
        indx_tot[iWSC][1] == -indx_tot[jWSC][1] &&
        indx_tot[iWSC][2] == -indx_tot[jWSC][2]
        )
        for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
          for (jWan = 0; jWan < StdI->NsiteUC; jWan++) {
            Mat_tot[iWSC][iWan][jWan] = 0.0;
          }
        }
    }/*for (jWSC = 0; jWSC < iWSC; jWSC++)*/
    if (indx_tot[iWSC][0] == 0 &&
      indx_tot[iWSC][1] == 0 &&
      indx_tot[iWSC][2] == 0)
      for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
        for (jWan = 0; jWan < iWan; jWan++) {
          Mat_tot[iWSC][iWan][jWan] = 0.0;
        }
      }
  }/*for (iWSC = 0; iWSC < nWSC; iWSC++)*/
  fclose(fp);
  /**@brief
  Then Store to the hopping Integeral and
  its site index.
  */
  nMat = 0;
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
      for (jWan = 0; jWan < StdI->NsiteUC; jWan++) {
        if (cutoff < cabs(Mat_tot[iWSC][iWan][jWan])) {
          for (ii = 0; ii < 3; ii++) Matindx[nMat][ii] = indx_tot[iWSC][ii];
          Matindx[nMat][3] = iWan;
          Matindx[nMat][4] = jWan;
          Mat[nMat] = Mat_tot[iWSC][iWan][jWan];
          nMat += 1;
        }
      }/*for (jWan = 0; jWan < StdI->NsiteUC; jWan++)*/
    }/*for (iWan = 0; iWan < StdI->NsiteUC; iWan++)*/
  }/*for (iWSC = 0; iWSC < nWSC; iWSC++)*/

  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < nWan; iWan++) {
      free(Mat_tot[iWSC][iWan]);
    }
    free(Mat_tot[iWSC]);
    free(indx_tot[iWSC]);
  }
  free(Mat_tot);
  free(indx_tot);
}/*static int read_W90(struct StdIntList *StdI, char *model)*//**
@brief Setup a Hamiltonian for the Wannier90 *_hr.dat
@author Mitsuaki Kawamura (The University of Tokyo)
*/
void StdFace_Wannier90(
  struct StdIntList *StdI//!<[inout]
)
{
  int isite, jsite, ispin, ntransMax, nintrMax;
  int iL, iW, iH, kCell, it, ii;
  double Jtmp[3][3] = { {0.0} };
  FILE *fp;
  double complex Cphase;
  double dR[3], *Uspin;
  int n_t, n_u, n_j;
  double complex *W90_t, *W90_j, *W90_u;
  int **t_indx, **u_indx, **j_indx;
  char filename[256];

  fprintf(stdout, "\n  @ Wannier90 Geometry \n\n");
  geometry_W90(StdI);

  /**@brief
  (1) Compute the shape of the super-cell and sites in the super-cell
  */
  fp = fopen("lattice.xsf", "w");
  /**/
  StdI->NsiteUC = 1;
  StdFace_InitSite(StdI, fp, 3);
  StdFace_PrintVal_d("phase0", &StdI->phase[0], 0.0);
  StdFace_PrintVal_d("phase1", &StdI->phase[1], 0.0);
  StdFace_PrintVal_d("phase2", &StdI->phase[2], 0.0);
  /*
  Read Hopping
  */
  fprintf(stdout, "\n  @ Wannier90 hopping \n\n");
  StdFace_PrintVal_d("cutoff_t", &StdI->cutoff_t, 1.0e-8);
  sprintf(filename, "%s_hr.dat", StdI->CDataFileHead);
  n_t = read_W90_query(StdI, filename, StdI->cutoff_t);
  W90_t = (double complex *)malloc(sizeof(double complex) * n_t);
  t_indx = (int **)malloc(sizeof(int*) * n_t);
  for (ii = 0; ii < n_t; ii++) t_indx[ii] = (int *)malloc(sizeof(int) * 5);
  read_W90(StdI, filename, StdI->cutoff_t, W90_t, t_indx);
  /*
  Read Coulomb
  */
  fprintf(stdout, "\n  @ Wannier90 Coulomb \n\n");
  StdFace_PrintVal_d("cutoff_u", &StdI->cutoff_u, 1.0e-8);
  sprintf(filename, "%s_ur.dat", StdI->CDataFileHead);
  n_u = read_W90_query(StdI, filename, StdI->cutoff_u);
  W90_u = (double complex *)malloc(sizeof(double complex) * n_u);
  u_indx = (int **)malloc(sizeof(int*) * n_u);
  for (ii = 0; ii < n_u; ii++) u_indx[ii] = (int *)malloc(sizeof(int) * 5);
  read_W90(StdI, filename, StdI->cutoff_t, W90_u, u_indx);
  /*
  Read Hund
  */
  fprintf(stdout, "\n  @ Wannier90 Hund \n\n");
  StdFace_PrintVal_d("cutoff_j", &StdI->cutoff_j, 1.0e-8);
  sprintf(filename, "%s_jr.dat", StdI->CDataFileHead);
  n_j = read_W90_query(StdI, filename, StdI->cutoff_j);
  W90_j = (double complex *)malloc(sizeof(double complex) * n_j);
  j_indx = (int **)malloc(sizeof(int*) * n_j);
  for (ii = 0; ii < n_j; ii++) j_indx[ii] = (int *)malloc(sizeof(int) * 5);
  read_W90(StdI, filename, StdI->cutoff_j, W90_j, j_indx);
  /**@brief
  (2) check & store parameters of Hamiltonian
  */
  fprintf(stdout, "\n  @ Hamiltonian \n\n");
  StdFace_NotUsed_d("K", StdI->K);
  StdFace_PrintVal_d("h", &StdI->h, 0.0);
  StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
  StdFace_NotUsed_d("U", StdI->U);
  /**/
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdFace_NotUsed_i("2S", StdI->S2);
  }/*if (strcmp(StdI->model, "spin") == 0 )*/
  else if (strcmp(StdI->model, "hubbard") == 0) {
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
  }
  else{
    printf("wannier + Kondo is not available !\n");
    StdFace_exit(-1);
  }/*if (model != "spin")*/
  fprintf(stdout, "\n  @ Numerical conditions\n\n");
  /**@brief
  (3) Set local spin flag (StdIntList::locspinflag) and 
  the number of sites (StdIntList::nsite)
  */
  StdI->nsite = StdI->NsiteUC * StdI->NCell;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  /**/
  if(strcmp(StdI->model, "spin") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++) StdI->locspinflag[isite] = StdI->S2;
  else if(strcmp(StdI->model, "hubbard") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++) StdI->locspinflag[isite] = 0;
  /**@brief
  (4) Compute the upper limit of the number of Transfer & Interaction and malloc them.
  */
  if (strcmp(StdI->model, "spin") == 0 ) {
    ntransMax = StdI->nsite * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
    nintrMax = StdI->NCell * (StdI->NsiteUC/*D*/ + n_t/*J*/ + n_u + n_j)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else if (strcmp(StdI->model, "hubbard") == 0) {
    ntransMax = StdI->NCell * 2/*spin*/ * (2 * StdI->NsiteUC/*mu+h+Gamma*/ + n_t * 2/*t*/);
    nintrMax = StdI->NCell * (n_u + n_j + 1);
  }
  /**/
  StdFace_MallocInteractions(StdI, ntransMax, nintrMax);
  /**@brief
  (4.5) For spin system, compute super exchange interaction.
  */
  if (strcmp(StdI->model, "spin") == 0) {
    Uspin = (double *)malloc(sizeof(double) * StdI->NsiteUC);
    for (it = 0; it < n_u; it++)
      if (u_indx[it][0] == 0 && u_indx[it][1] == 0 && u_indx[it][2] == 0
        && u_indx[it][3] == u_indx[it][4])     
        Uspin[u_indx[it][3]] = creal(W90_u[it]);
  }/*if (strcmp(StdI->model, "spin") == 0)*/
  /**@brief
  (5) Set Transfer & Interaction
  */
  for (kCell = 0; kCell < StdI->NCell; kCell++){
    /**/
    iW = StdI->Cell[kCell][0];
    iL = StdI->Cell[kCell][1];
    iH = StdI->Cell[kCell][2];
    /*
     Local term 1
    */
    if (strcmp(StdI->model, "spin") == 0) {
      for (isite = StdI->NsiteUC*kCell; isite < StdI->NsiteUC*(kCell + 1); isite++) {
        StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, isite);
      }
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      for (isite = StdI->NsiteUC*kCell; isite < StdI->NsiteUC*(kCell + 1); isite++) {
        StdFace_HubbardLocal(StdI, StdI->mu, -StdI->h, -StdI->Gamma, 0.0, isite);
      }
    }/*if (strcmp(StdI->model, "spin") != 0 )*/
    /*
     Hopping
    */
    for (it = 0; it < n_t; it++) {
      /*
       Local term
      */
      if (t_indx[it][0] == 0 && t_indx[it][1] == 0 && t_indx[it][2] == 0
        && t_indx[it][3] == t_indx[it][4])
      {
        if (strcmp(StdI->model, "hubbard") == 0) {
          isite = StdI->NsiteUC*kCell + t_indx[it][3];
          for (ispin = 0; ispin < 2; ispin++) {
            StdI->trans[StdI->ntrans] = -W90_t[it];
            StdI->transindx[StdI->ntrans][0] = isite;
            StdI->transindx[StdI->ntrans][1] = ispin;
            StdI->transindx[StdI->ntrans][2] = isite;
            StdI->transindx[StdI->ntrans][3] = ispin;
            StdI->ntrans = StdI->ntrans + 1;
          }/*for (ispin = 0; ispin < 2; ispin++)*/
        }/*if (strcmp(StdI->model, "hubbrad") == 0 )*/
      }/*Local term*/
      else {
        /*
         Non-local term
        */
        StdFace_FindSite(StdI, iW, iL, iH,
          t_indx[it][0], t_indx[it][1], t_indx[it][2],
          t_indx[it][3], t_indx[it][4], &isite, &jsite, &Cphase, dR);
        if (strcmp(StdI->model, "spin") == 0) {
          for (ii = 0; ii < 3; ii++) 
            Jtmp[ii][ii] = 2.0 * W90_t[it] * conj(W90_t[it])
            * (1.0 / Uspin[t_indx[it][3]] + 1.0 / Uspin[t_indx[it][4]]);
          StdFace_GeneralJ(StdI, Jtmp, StdI->S2, StdI->S2, isite, jsite);
        }/*if (strcmp(StdI->model, "spin") == 0 )*/
        else {
          StdFace_Hopping(StdI, - Cphase * W90_t[it], isite, jsite, dR);
        }
      }/*Non-local term*/
    }/*for (it = 0; it < n_t; it++)*/
    /*
    Coulomb integral (U)
    */
    for (it = 0; it < n_u; it++) {
      /*
      Local term
      */
      if (u_indx[it][0] == 0 && u_indx[it][1] == 0 && u_indx[it][2] == 0
        && u_indx[it][3] == u_indx[it][4])
      {
        StdI->Cintra[StdI->NCintra] = creal(W90_u[it]);
        StdI->CintraIndx[StdI->NCintra][0] = StdI->NsiteUC*kCell + u_indx[it][3];
        StdI->NCintra += 1;
      }/*Local term*/
      else {
        /*
        Non-local term
        */
        StdFace_FindSite(StdI, iW, iL, iH,
          u_indx[it][0], u_indx[it][1], u_indx[it][2],
          u_indx[it][3], u_indx[it][4], &isite, &jsite, &Cphase, dR);
        StdFace_Coulomb(StdI, creal(W90_u[it]), isite, jsite);
      }/*Non-local term*/
    }/*for (it = 0; it < n_t; it++)*/
    /*
     Hund coupling (J)
    */
    for (it = 0; it < n_j; it++) {
      /*
      Local term should not be computed
      */
      if (j_indx[it][0] != 0 || j_indx[it][1] != 0 || j_indx[it][2] != 0
        || j_indx[it][3] != j_indx[it][4])
      {
        StdFace_FindSite(StdI, iW, iL, iH,
          j_indx[it][0], j_indx[it][1], j_indx[it][2],
          j_indx[it][3], j_indx[it][4], &isite, &jsite, &Cphase, dR);

        StdI->Hund[StdI->NHund] = creal(W90_j[it]);
        StdI->HundIndx[StdI->NHund][0] = isite;
        StdI->HundIndx[StdI->NHund][1] = jsite;
        StdI->NHund += 1;

        StdI->Ex[StdI->NEx] = creal(W90_j[it]);
        StdI->ExIndx[StdI->NEx][0] = isite;
        StdI->ExIndx[StdI->NEx][1] = jsite;
        StdI->NEx += 1;

        if (strcmp(StdI->model, "hubbard") == 0) {
          StdI->PairHopp[StdI->NPairHopp] = creal(W90_j[it]);
          StdI->PHIndx[StdI->NPairHopp][0] = isite;
          StdI->PHIndx[StdI->NPairHopp][1] = jsite;
          StdI->NPairHopp += 1;
        }
      }/*Non-local term*/
    }/*for (it = 0; it < n_t; it++)*/
  }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/

  fclose(fp);
  StdFace_PrintXSF(StdI);
  StdFace_PrintGeometry(StdI);

  for (it = 0; it < n_t; it++) free(t_indx[it]);
  free(t_indx);
  free(W90_t);
  for (it = 0; it < n_u; it++) free(u_indx[it]);
  free(u_indx);
  free(W90_u); 
  for (it = 0; it < n_j; it++) free(j_indx[it]);
  free(j_indx);
  free(W90_j); 
  if (strcmp(StdI->model, "spin") == 0) free(Uspin);

}/*void StdFace_Wannier90*/
