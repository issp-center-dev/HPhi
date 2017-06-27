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
#include "input.h"
#include "FileIO.h"
#include "wrapperMPI.h"

int inputHam(struct BindStruct *X){
  //Input Ham
  long int i=0;
  long int ham_i=0;
  long int ham_j=0;
  long int imax = X->Check.idim_max;
  long int ihermite=0;
  long int itmp;
  double dHam_re, dHam_im;
  char ctmp[256], ctmp2[256];

  FILE *fp;
  char sdt[D_FileNameMax];

  sprintf(sdt,cFileNamePhys_FullDiag_Ham, X->Def.CDataFileHead);
  if(childfopenMPI(sdt,"r",&fp)!=0){
    return -1;
  }

  //skip: header
  fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
  //skip: read imax, imax, ihermite
  fgetsMPI(ctmp, sizeof(ctmp) / sizeof(char), fp);
  sscanf(ctmp, "%ld %ld %ld\n", &itmp, &itmp, &ihermite);
  if(itmp != imax){
    fprintf(stdoutMPI, "Error: The dimension of input Hamiltonian is wrong: input=%ld, idim=%ld.\n", itmp, imax);
    return -1;
  }
  for(i=1; i<= ihermite; i++){
    fgetsMPI(ctmp2, sizeof(ctmp2) / sizeof(char), fp);
    sscanf(ctmp2, "%ld %ld %lf %lf\n",
           &ham_i, &ham_j, &dHam_re, &dHam_im);
    Ham[ham_i][ham_j]=dHam_re+I*dHam_im;
    Ham[ham_j][ham_i]=conj(Ham[ham_i][ham_j]);
  }
 fclose(fp);
 return 0;
}
