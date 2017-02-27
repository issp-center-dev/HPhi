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
* Read Geometry file for wannier90
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
static void geometry_W90(struct StdIntList *StdI, int *wan_num) {
  int isite, ii, ierr;
  FILE *fp;

  fprintf(stdout, "                Wannier90 Geometry file = %s\n", StdI->W90_geom);

  fp = fopen(StdI->W90_geom, "r");
  /*
   Direct lattice vector
  */
  for (ii = 0; ii < 3; ii++)
    ierr = fscanf(fp, "%lf%lf%lf", &StdI->direct[ii][0], &StdI->direct[ii][1], &StdI->direct[ii][2]);
  /*
   Site position
  */
  for (isite = 0; isite < StdI->NsiteUC; isite++) free(StdI->tau[isite]);
  free(StdI->tau);
  ierr = fscanf(fp, "%d", &StdI->NsiteUC);
  fprintf(stdout, "              Number of EFFECTIVE Sites = %d\n", StdI->NsiteUC);

  StdI->tau = (double **)malloc(sizeof(double*) * StdI->NsiteUC);
  for (ii = 0; ii < StdI->NsiteUC; ii++) StdI->tau[ii] = (double *)malloc(sizeof(double) * 3);

  for (isite = 0; isite < StdI->NsiteUC; isite++)
    ierr = fscanf(fp, "%d%lf%lf%lf", &wan_num[isite],
      &StdI->tau[isite][0], &StdI->tau[isite][1], &StdI->tau[isite][2]);
  fclose(fp);

  printf("      Direct lattice vectors:\n");
  for (ii = 0; ii < 3; ii++) printf("      %10.5f %10.5f %10.5f\n",
    StdI->direct[ii][0], StdI->direct[ii][1], StdI->direct[ii][2]);
  printf("      Wannier centres:\n");
  for (isite = 0; isite < StdI->NsiteUC; isite++) printf("      %10.5f %10.5f %10.5f\n",
    StdI->tau[isite][0], StdI->tau[isite][1], StdI->tau[isite][2]);

}/*static void geometry_W90(struct StdIntList *StdI) */
/*
* Read Wannier90 hamiltonian file
*
* @author Mitsuaki Kawamura (The University of Tokyo)
*/
static void read_W90(struct StdIntList *StdI, char *model) 
{
  FILE *fp;
  int ierr, nWan, nWSC, iWSC, jWSC, iWan, jWan, iWan0, jWan0, ii;
  double dtmp[2], tmax, tabs;
  char ctmp[256], *ctmp2;
  double complex ***t_tot, **t0;
  int **indx_tot, *wan_num;

  fprintf(stdout, "\n  @ Wannier90 Mode\n\n");

  StdFace_PrintVal_d("W90_thr", &StdI->W90_cutoff, 0.1);

  fprintf(stdout, "                         Wannier90 file = %s\n", StdI->W90_hr);
  
  fp = fopen(StdI->W90_hr, "r");
  ctmp2 = fgets(ctmp, 256, fp);
  ierr = fscanf(fp, "%d", &nWan);
  ierr = fscanf(fp, "%d", &nWSC);
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    ierr = fscanf(fp, "%d", &ii);
  }
  fprintf(stdout, "                      Number of Wannier = %d\n", nWan);
  fprintf(stdout, "            Number of Wigner-Seitz Cell = %d\n", nWSC);

  t_tot = (double complex ***)malloc(sizeof(double complex **) * nWSC);
  indx_tot = (int **)malloc(sizeof(int*) * nWSC);
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    t_tot[iWSC] = (double complex **)malloc(sizeof(double complex *) * nWan);
    indx_tot[iWSC] = (int *)malloc(sizeof(int) * 3);
    for (iWan = 0; iWan < nWan; iWan++) {
      t_tot[iWSC][iWan] = (double complex *)malloc(sizeof(double complex) * nWan);
    }
  }
  t0 = (double complex **)malloc(sizeof(double complex *) * nWan);
  for (iWan = 0; iWan < nWan; iWan++)
    t0[iWan] = (double complex *)malloc(sizeof(double complex) * nWan);
  wan_num = (int *)malloc(sizeof(int) * nWan);

  geometry_W90(StdI, wan_num);

  tmax = 0.0;
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < nWan; iWan++) {
      for (jWan = 0; jWan < nWan; jWan++) {
        ierr = fscanf(fp, "%d%d%d%d%d%lf%lf",
          &indx_tot[iWSC][0], &indx_tot[iWSC][1], &indx_tot[iWSC][2], &iWan0, &jWan0,
          &dtmp[0], &dtmp[1]);
        t0[iWan0 - 1][jWan0 - 1] = dtmp[0] + I * dtmp[1];
      }
    }
    for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
      for (jWan = 0; jWan < StdI->NsiteUC; jWan++) {
        t_tot[iWSC][iWan][jWan] = t0[wan_num[iWan]][wan_num[jWan]];
      }
    }
    /*
    Inversion symmetry
    */
    for (jWSC = 0; jWSC < iWSC; jWSC++) {
      if (
        indx_tot[iWSC][0] == -indx_tot[jWSC][0] &&
        indx_tot[iWSC][1] == -indx_tot[jWSC][1] &&
        indx_tot[iWSC][2] == -indx_tot[jWSC][2]
        )
        for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
          for (jWan = 0; jWan < StdI->NsiteUC; jWan++) {
            t_tot[iWSC][iWan][jWan] = 0.0;
          }
        }
    }/*for (jWSC = 0; jWSC < iWSC; jWSC++)*/
    if (indx_tot[iWSC][0] == 0 &&
        indx_tot[iWSC][1] == 0 &&
        indx_tot[iWSC][2] == 0) 
      for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
        for (jWan = 0; jWan < iWan; jWan++) {
          t_tot[iWSC][iWan][jWan] = 0.0;
        }
      }
    /*
     Max t
    */
    for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
      for (jWan = 0; jWan < StdI->NsiteUC; jWan++) {
        tabs = cabs(t_tot[iWSC][iWan][jWan]);
        if (tmax < tabs) tmax = tabs;
      }
    }
  }/*for (iWSC = 0; iWSC < nWSC; iWSC++)*/
  fclose(fp);

  fprintf(stdout, "                        Maximum Hopping = %f\n", tmax);
  fprintf(stdout, "                  Threshold for Hopping = %f\n", tmax * StdI->W90_cutoff);
  /*
    Cut-off of Hopping
  */
  /*  Query  */
  StdI->W90_nt = 0;
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
      for (jWan = 0; jWan < StdI->NsiteUC; jWan++) {
        if (tmax * StdI->W90_cutoff < cabs(t_tot[iWSC][iWan][jWan])) StdI->W90_nt += 1;
      }
    }
  }
  fprintf(stdout, "      Total number of EFFECTIVE Hopping = %d\n", StdI->W90_nt);

  StdI->W90_t = (double complex *)malloc(sizeof(double complex) * StdI->W90_nt);
  StdI->W90_indx = (int **)malloc(sizeof(int*) * StdI->W90_nt);
  for (ii = 0; ii < StdI->W90_nt; ii++) StdI->W90_indx[ii] = (int *)malloc(sizeof(int) * 5);

  /*  Store  */
  fprintf(stdout, "      EFFECTIVE Hoppings:\n");
  StdI->W90_nt = 0;
  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < StdI->NsiteUC; iWan++) {
      for (jWan = 0; jWan < StdI->NsiteUC; jWan++) {
        if (tmax * StdI->W90_cutoff < cabs(t_tot[iWSC][iWan][jWan])) {
          for (ii = 0; ii < 3; ii++) StdI->W90_indx[StdI->W90_nt][ii] = indx_tot[iWSC][ii];
          StdI->W90_indx[StdI->W90_nt][3] = iWan;
          StdI->W90_indx[StdI->W90_nt][4] = jWan;
          StdI->W90_t[StdI->W90_nt] = t_tot[iWSC][iWan][jWan];
          fprintf(stdout, "        %5d%5d%5d%5d%5d%12.6f%12.6f\n", 
            StdI->W90_indx[StdI->W90_nt][0], StdI->W90_indx[StdI->W90_nt][1], StdI->W90_indx[StdI->W90_nt][2], 
            StdI->W90_indx[StdI->W90_nt][3], StdI->W90_indx[StdI->W90_nt][4], 
            creal(StdI->W90_t[StdI->W90_nt]), cimag(StdI->W90_t[StdI->W90_nt]));
          StdI->W90_nt += 1;
        }
      }
    }
  }

  for (iWSC = 0; iWSC < nWSC; iWSC++) {
    for (iWan = 0; iWan < nWan; iWan++) {
      free(t_tot[iWSC][iWan]);
    }
    free(t_tot[iWSC]);
    free(indx_tot[iWSC]);
  }
  free(t_tot);
  free(indx_tot);
  for (iWan = 0; iWan < nWan; iWan++) free(t0[iWan]);
  free(t0);
  free(wan_num);
}/*static void read_W90(struct StdIntList *StdI, char *model)*/
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
  StdI->NsiteUC = 1;
  StdFace_InitSite(StdI, fp, 3);
  StdFace_PrintVal_d("phase0", &StdI->phase[0], 0.0);
  StdFace_PrintVal_d("phase1", &StdI->phase[1], 0.0);
  StdFace_PrintVal_d("phase2", &StdI->phase[2], 0.0);
  /**/
  read_W90(StdI, model);
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
    StdI->nintr = StdI->NCell * (StdI->NsiteUC/*D*/ + StdI->W90_nt/*J*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else if (strcmp(StdI->model, "hubbard") == 0) {
    StdI->ntrans = StdI->NCell * 2/*spin*/ * (StdI->NsiteUC/*mu*/ + StdI->W90_nt * 2/*t*/);
    StdI->nintr = StdI->NCell * StdI->NsiteUC/*U*/;
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
