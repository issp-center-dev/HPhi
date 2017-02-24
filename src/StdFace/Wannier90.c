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
#include "StdFace_vals.h"
#include "StdFace_ModelUtil.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>

/*
* Read Wannier90 hamiltonian file
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
static void read_W90(struct StdIntList *StdI, char *model) 
{
  FILE *fp;
  int ierr, it, jt, nt_tot;
  double dtmp[2], tmax, tabs;
  char ctmp[256], *ctmp2;
  double complex *t_tot;
  int **indx_tot;

  fprintf(stdout, "\n  @ Wannier90 Mode\n\n");

  StdFace_PrintVal_d("W90_thr", &StdI->W90_cutoff, 0.1);

  fprintf(stdout, "                       Wannier90 file = %s\n", StdI->W90_file);
  
  fp = fopen(StdI->W90_file, "r");
  ctmp2 = fgets(ctmp, 256, fp);
  ierr = fscanf(fp, "%d", &StdI->NsiteUC);
  ierr = fscanf(fp, "%d", &nt_tot);
  for (it = 0; it < nt_tot; it++) {
    ierr = fscanf(fp, "%d", &jt);
  }
  nt_tot *= StdI->NsiteUC * StdI->NsiteUC;
  t_tot = (double complex *)malloc(sizeof(double complex) * nt_tot);
  indx_tot = (int **)malloc(sizeof(int*) * nt_tot);
  for (it = 0; it < nt_tot; it++) indx_tot[it] = (int *)malloc(sizeof(int) * 5);

  tmax = 0.0;
  for (it = 0; it < nt_tot; it++) {
    ierr = fscanf(fp, "%d%d%d%d%d%lf%lf",
      &indx_tot[it][0], &indx_tot[it][1], &indx_tot[it][2], &indx_tot[it][3], &indx_tot[it][4],
      &dtmp[0], &dtmp[1]);
    t_tot[it] = dtmp[0] + I * dtmp[1];
    tabs = cabs(t_tot[it]);
    if (tmax < tabs) tmax = tabs;
    /*
    Inversion symmetry
    */
    for (jt = 0; jt < it; jt++) {
      if (
        indx_tot[it][0] == -indx_tot[jt][0] &&
        indx_tot[it][1] == -indx_tot[jt][1] &&
        indx_tot[it][2] == -indx_tot[jt][2] &&
        indx_tot[it][3] == indx_tot[jt][4] &&
        indx_tot[it][4] == indx_tot[jt][3]
        )
        t_tot[it] = 0.0;
    }/*for (jt = 0; jt < it; jt++)*/
  }/*for (ii = 0; ii < nt_tot; ii++)*/
  fclose(fp);

  fprintf(stdout, "              Total number of Hopping = %d\n", nt_tot);
  fprintf(stdout, "                      Maximum Hopping = %f\n", tmax);
  fprintf(stdout, "                Threshold for Hopping = %f\n", tmax * StdI->W90_cutoff);
  /*
    Cut-off of Hopping
  */
  /*  Query  */
  StdI->W90_nt = 0;
  for (it = 0; it < nt_tot; it++) {
    if (tmax * StdI->W90_cutoff < cabs(t_tot[it])) StdI->W90_nt += 1;
  }/*for (it = 0; it < nt_tot; it++))*/
  fprintf(stdout, "    Total number of EFFECTIVE Hopping = %d\n", StdI->W90_nt);

  StdI->W90_t = (double complex *)malloc(sizeof(double complex) * StdI->W90_nt);
  StdI->W90_indx = (int **)malloc(sizeof(int*) * StdI->W90_nt);
  for (it = 0; it < StdI->W90_nt; it++) StdI->W90_indx[it] = (int *)malloc(sizeof(int) * 5);

  /*  Store  */
  StdI->W90_nt = 0;
  for (it = 0; it < nt_tot; it++) {
    if (tmax * StdI->W90_cutoff< cabs(t_tot[it])) {
      for (jt = 0; jt < 3; jt++) StdI->W90_indx[StdI->W90_nt][jt] = indx_tot[it][jt];
      for (jt = 3; jt < 5; jt++) StdI->W90_indx[StdI->W90_nt][jt] = indx_tot[it][jt] - 1;
      StdI->W90_t[StdI->W90_nt] = t_tot[it];
    }
  }/*for (it = 0; it < nt_tot; it++))*/

  for (it = 0; it < nt_tot; it++) free(indx_tot[it]);
  free(indx_tot);
  free(t_tot);

}/*static void read_W90(struct StdIntList *StdI, char *model)*/
/*
 * Read Geometry file for wannier90
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
static void geometry_W90(struct StdIntList *StdI) {
  int isite, ii, jj, ierr;
  FILE *fp;
  double vcell, recipr[3][3], tau0[3];
  
  if (strcmp(StdI->W90_geom, "****") == 0) {
    for (ii = 0; ii < 3; ii++) {
      for (jj = 0; jj < 3; jj++) {
        StdI->direct[ii][jj] = 0.0;
      }
      StdI->direct[ii][ii] = 1.0;
    }
    for (isite = 0; isite < StdI->NsiteUC; isite++)
      for (ii = 0; ii < 3; ii++) StdI->tau[isite][ii] = 0.0;
  }/*if (strcmp(StdI->StdI->W90_geom, "****") == 0)*/
  else {
    fprintf(stdout, "                  Wannier90 Geometry file = %s\n", StdI->W90_geom);

    fp = fopen(StdI->W90_geom, "r");
    for (ii = 0; ii < 3; ii++)
      ierr = fscanf(fp, "%lf%lf%lf", &StdI->direct[ii][0], &StdI->direct[ii][1], &StdI->direct[ii][2]);
    for (isite = 0; isite < StdI->NsiteUC; isite++)
      ierr = fscanf(fp, "%lf%lf%lf", &StdI->tau[isite][0], &StdI->tau[isite][1], &StdI->tau[isite][2]);
    fclose(fp);
    /*
    Reciprocal lattice vector
    */
    vcell = 0.0;
    for (ii = 0; ii < 3; ii++) {
      vcell += StdI->direct[0][ii]
             * StdI->direct[1][(ii + 1) % 3]
             * StdI->direct[2][(ii + 2) % 3]
             - StdI->direct[0][ii]
             * StdI->direct[1][(ii + 2) % 3]
             * StdI->direct[2][(ii + 1) % 3];
    }

    for (ii = 0; ii < 3; ii++) {
      for (jj = 0; jj < 3; jj++) {
        recipr[ii][jj] = StdI->direct[(ii + 1) % 3][(jj + 1) % 3] * StdI->direct[(ii + 2) % 3][(jj + 2) % 3]
                       - StdI->direct[(ii + 1) % 3][(jj + 2) % 3] * StdI->direct[(ii + 2) % 3][(jj + 1) % 3];
        recipr[ii][jj] /= vcell;
      }/*for (jj = 0; jj < 3; jj++)*/
    }/*for (jj = 0; jj < 3; jj++)*/

    for (isite = 0; isite < StdI->NsiteUC; isite++){
      for (ii = 0; ii < 3; ii++) {
        tau0[ii] = 0.0;
        for (jj = 0; jj < 3; jj++) tau0[ii] += recipr[ii][jj] * StdI->tau[isite][jj];
      }
      for (ii = 0; ii < 3; ii++) StdI->tau[isite][ii] = tau0[ii];
    }/*for (isite = 0; isite < StdI->NsiteUC; isite++)*/

  }/*if (strcmp(StdI->StdI->W90_geom, "****") != 0)*/

  printf("      Direct lattice vectors:\n");
  for (ii = 0; ii < 3; ii++) printf("      %10.5f %10.5f %10.5f\n",
    StdI->direct[ii][0], StdI->direct[ii][1], StdI->direct[ii][2]);
  printf("      Wannier centres:\n");
  for (isite = 0; isite < StdI->NsiteUC; isite++) printf("      %10.5f %10.5f %10.5f\n",
    StdI->tau[isite][0], StdI->tau[isite][1], StdI->tau[isite][2]);

}/*static void geometry_W90(struct StdIntList *StdI) */
/**
 *
 * Setup a Hamiltonian for the Wannier90 *_hr.dat
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_Wannier90(struct StdIntList *StdI, char *model)
{
  int isite, jsite;
  int iL, iW, iH, kCell, it, ii;
  double Jtmp[3][3] = {0.0};
  FILE *fp;
  double complex Cphase;

  fprintf(stdout, "\n");
  fprintf(stdout, "#######  Parameter Summary  #######\n");
  fprintf(stdout, "\n");
  /*
   Initialize Cell
  */
  read_W90(StdI, model);
  StdFace_InitSite(StdI, fp, 3);
  geometry_W90(StdI);
  /**/
  StdFace_PrintVal_d("phase0", &StdI->phase[0], 0.0);
  StdFace_PrintVal_d("phase1", &StdI->phase[1], 0.0);
  StdFace_PrintVal_d("phase2", &StdI->phase[2], 0.0);
  /**/
  fprintf(stdout, "\n  @ Hamiltonian \n\n");
  StdFace_NotUsed_d("K", StdI->K);
  /**/
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdFace_NotUsed_i("2S", StdI->S2);
    StdFace_PrintVal_d("h", &StdI->h, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 1.0);
  }/*if (strcmp(StdI->model, "spin") == 0 )*/
  else if (strcmp(StdI->model, "hubbard") == 0) {
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 0.0);
  }
  else{
    StdFace_exit(-1);
  }/*if (model != "spin")*/
  fprintf(stdout, "\n  @ Numerical conditions\n\n");
  /*
   Local Spin
  */
  StdI->nsite = StdI->NsiteUC * StdI->NCell;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  /**/
  if(strcmp(StdI->model, "spin") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++) StdI->locspinflag[isite] = StdI->S2;
  else if(strcmp(StdI->model, "hubbard") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++) StdI->locspinflag[isite] = 0;
  /*
   The number of Transfer & Interaction
  */
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdI->ntrans = StdI->nsite * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
    StdI->nintr = StdI->NCell * (StdI->NsiteUC/*D*/ + 12/*J*/ + 0/*J'*/ + 0/*J''*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else {
    StdI->ntrans = StdI->NCell * 2/*spin*/ * (StdI->NsiteUC/*mu*/ + 24/*t*/ + 0/*t'*/ + 0/*t''*/);
    StdI->nintr = StdI->NCell * (StdI->NsiteUC/*U*/ + 4 * (12/*V*/ + 0/*V'*/ + 0/*V''*/));

    if (strcmp(StdI->model, "kondo") == 0 )  StdI->nintr += 
      StdI->nsite / 2 * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  /**/
  StdFace_MallocInteractions(StdI);
  /*
   Set Transfer & Interaction
  */
  StdI->ntrans = 0;
  StdI->nintr = 0;
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
        StdFace_Hopping(StdI, StdI->mu, isite, isite, 0);
        StdI->Cintra[StdI->NCintra] = StdI->U; StdI->CintraIndx[StdI->NCintra][0] = isite; StdI->NCintra += 1;
      }
    }/*if (strcmp(StdI->model, "spin") != 0 )*/
    /**/
    for (it = 0; it < StdI->W90_nt; it++) {
      if (StdI->W90_indx[it][0] == 0 && StdI->W90_indx[it][1] == 0 && StdI->W90_indx[it][2] == 0
        && StdI->W90_indx[it][3] == StdI->W90_indx[it][4])
      {
        if (strcmp(StdI->model, "hubbard") == 0) {
          isite = StdI->NsiteUC*kCell + StdI->W90_indx[it][3];
          StdFace_Hopping(StdI, StdI->W90_t[it], isite, isite, 0);
        }/*if (strcmp(StdI->model, "hubbrad") == 0 )*/
      }/*Local term*/
      else {
        StdFace_FindSite(StdI, iW, iL, iH,
          StdI->W90_indx[it][0], StdI->W90_indx[it][1], StdI->W90_indx[it][2],
          StdI->W90_indx[it][3], StdI->W90_indx[it][4], &isite, &jsite, &Cphase);
        if (strcmp(StdI->model, "spin") == 0) {
          for (ii = 0; ii < 3; ii++) Jtmp[ii][ii] = StdI->W90_t[it] * conj(StdI->W90_t[it]) / StdI->U;
          StdFace_GeneralJ(StdI, Jtmp, StdI->S2, StdI->S2, isite, jsite);
        }/*if (strcmp(StdI->model, "spin") == 0 )*/
        else {
          StdFace_Hopping(StdI, Cphase * StdI->W90_t[it], isite, jsite, 1);
        }
      }/*Non-local term*/
    }/*for (it = 0; it < StdI->W90_nt; it++)*/
  }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/

  StdFace_PrintXSF(StdI);
  StdFace_PrintGeometry(StdI);

  for (it = 0; it < StdI->W90_nt; it++) free(StdI->W90_indx[it]);
  free(StdI->W90_indx);
  free(StdI->W90_t);

}/*void StdFace_Pyrochlore*/
