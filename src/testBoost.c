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
/*-------------------------------------------------------------*/

#include <sz.h>
#include <HPhiTrans.h>
#include <output_list.h>
#include <diagonalcalc.h>
#include <CalcByLanczos.h>
#include <CalcByFullDiag.h>
#include <CalcByTPQ.h>
#include <check.h>
#include "Common.h"
#include "readdef.h"
#include "StdFace_main.h"
#include "wrapperMPI.h"
#include "mfmemory.h"
#include "xsetmem.h"


/** 
 * @test program for HPhi Boost mode
 * 
 * @param argc 
 * @param argv 
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * @author Youhei Yamaji (The University of Tokyo)
 * 
 * @return 
 */

int main(){
  
  int L;
  int j, k;
  double tmpd;
  double complex **arrayJ, **matJ;

  char *filename = "test_input";
  FILE *fp;

  if((fp = fopen(filename, "r")) == NULL){
    fprintf(stderr, "failed to open a file %s\n", filename);
    exit(EXIT_FAILURE);
  }
  

  L = (int)pow(2,6); 
  printf("%d\n", L);

  c_malloc2(arrayJ,3,3);
  c_malloc2(matJ,L,L);

  for(j=0; j<3; j++){
    for(k=0; k<3; k++){
      arrayJ[j][k] = 0.0 + 0.0*I;
    }
  }
 
  arrayJ[0][0] = 1.0 + 0.0*I; 
  arrayJ[0][1] =     - 0.5*I; 
  arrayJ[1][0] =     + 0.5*I;

  fscanf(fp, "%lf", &tmpd);
  arrayJ[0][0] = arrayJ[0][0]*tmpd;

  fclose(fp);

  printf("%f %f\n", creal(arrayJ[0][0]), cimag(arrayJ[0][0]));
  printf("%f %f\n", creal(arrayJ[0][1]), cimag(arrayJ[0][1]));
}
