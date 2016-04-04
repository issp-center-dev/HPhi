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
#include "matrixlapack.h"
#include <stdlib.h>
#include <omp.h>

void zgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double complex *ALPHA, double complex *matJL, int *LDA, double complex *arrayz, int *LDB, double complex *BETA, double complex *arrayx, int *LDC);

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

  double complex *arrayx, *arrayy, *arrayz, *arrayw;

  c_malloc1(arrayx,2);
  c_malloc1(arrayz,2);

//  arrayx[0]=0.0;
//  arrayx[1]=0.0;
  arrayz[0]=0.0;
  arrayz[1]=0.0;

  #pragma omp parallel default(none) private(arrayy,j) \
  firstprivate(arrayz) //	firstprivate(arrayz,arrayw)
  {
    c_malloc1(arrayy,2);
    //c_malloc1(arrayw,2);

    arrayy[0]=0.0;
    arrayy[1]=0.0;
    //arrayw[0]=0.0;
    //arrayw[1]=0.0;

    for(j=0;j<2;j++){
      arrayx[j] += 1.0;
      arrayy[j] += 1.0;
      arrayz[j] += 1.0;
      //arrayw[j] += 1.0;
    } 

//    for(j=0;j<2;j++){printf("arrayx %d %lf\n",j,creal(arrayx[j]));}
    for(j=0;j<2;j++){printf("arrayy %d %lf\n",j,creal(arrayy[j]));}
    for(j=0;j<2;j++){printf("arrayz %d %lf\n",j,creal(arrayz[j]));}
    //for(j=0;j<2;j++){printf("arrayw %d %lf\n",j,creal(arrayw[j]));}

    c_free1(arrayy,2);
    //c_free1(arrayw,2);
  }

  
  c_free1(arrayx,2);
  c_free1(arrayz,2);
  
}
