/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 Takahiro Misawa, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Youhei Yamaji, Synge Todo, Naoki Kawashima */

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
#include "phys.h"
#include "wrapperMPI.h"

/** 
 * 
 * 
 * @param X 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
void phys(struct BindStruct *X){

    long unsigned int i,j,i_max;
    double tmp_N,tmp_Sz;

    i_max=X->Check.idim_max;
    for(i=0;i<i_max;i++){
      for(j=0;j<i_max;j++){
	v0[j+1] =  L_vec[i][j];
      }
      X->Phys.eigen_num=i;

      if(!expec_energy(X)==0){
	fprintf(stdoutMPI, "Error: calc expec_energy.\n");
	exitMPI(-1);
      }
      if(!expec_cisajs(X,v1)==0){
	fprintf(stdoutMPI, "Error: calc OneBodyG.\n");
	exitMPI(-1);
      }
      if(!expec_cisajscktaltdc(X, v1)==0){
	fprintf(stdoutMPI, "Error: calc TwoBodyG.\n");
	exitMPI(-1);
      }
      
      if(!expec_totalspin(X, v1)==0){
	fprintf(stdoutMPI, "Error: calc TotalSpin.\n");
	exitMPI(-1);
      }
      tmp_N  = X->Phys.num_up + X->Phys.num_down;
      tmp_Sz = X->Phys.num_up - X->Phys.num_down;
      fprintf(stdoutMPI, "i=%5ld Energy=%10lf N=%10lf Sz=%10lf S2=%10lf Doublon=%10lf \n",i,X->Phys.energy,tmp_N,tmp_Sz,X->Phys.s2,X->Phys.doublon);
      X->Phys.all_energy[i]   = X->Phys.energy;
      X->Phys.all_doublon[i]  = X->Phys.doublon;
      X->Phys.all_s2[i]       = X->Phys.s2;
      X->Phys.all_num_up[i]   = X->Phys.num_up;
      X->Phys.all_num_down[i] = X->Phys.num_down;
    }
}
