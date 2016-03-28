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

#ifdef MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include "Common.h"
#include "mfmemory.h"
#include "xsetmem.h"
#include "wrapperMPI.h"
#include "defmodelBoost.h"
/**
 *
 * Creating list and exchange couplings for Boost mode  
 *
 * @author Youhei Yamaji (The University of Tokyo)
 */
void defmodelBoost(
  long unsigned int W0,
  long unsigned int R0,
  long unsigned int num_pivot,
  long unsigned int ishift_nspin,
  int  **list_6spin_star,
  int ***list_6spin_pair,
  long unsigned int model_num,
      /* 1: Kitaev on standard lattice  
         2: Kagome on standard lattice ?  
         3: Kagome on a cluster equivalent to cluster (E) in JPSJ 79, 053707 (2010)   
      */
  double complex ***arrayJ,
  double complex *vecB
)
{
  char *filename = "arrayJ";
  FILE *fp;
  char *filenameb = "vecB";
  FILE *fpb;
  long unsigned int j, k, ell, NJ; 
  long unsigned int i0,i1,i2; 
  double ReJex, ImJex;
  double ReBx, ImBx, ReBy, ImBy, ReBz, ImBz;
  double complex *tmpJ;

  c_malloc1(tmpJ, 27);

  if((fp = fopen(filename, "r")) == NULL){
    fprintf(stderr, "\n ###Boost### failed to open a file %s\n", filename);
    exitMPI(EXIT_FAILURE);
  }
  if((fpb = fopen(filenameb, "r")) == NULL){
    fprintf(stderr, "\n ###Boost### failed to open a file %s\n", filenameb);
    exitMPI(EXIT_FAILURE);
  }

  for(j=0; j < 3; j++){
    for(k=0; k < 3; k++){
      for(ell=0; ell < 3; ell++){
        arrayJ[j][k][ell] = 0.0;
        tmpJ[j+3*k+9*ell];
      }
    }
  }
  if(myrank==0){
    fscanf(fp, "%ld\n", &NJ);
    for(j=0;j<NJ;j++){
      fscanf(fp, "%ld %ld %ld %lf %lf\n", &i0, &i1, &i2, &ReJex, &ImJex);
//      printf(" ###Boost### ReJex %lf\n",ReJex);
      arrayJ[i0][i1][i2] = ReJex + I*ImJex;
    }
    for(j=0; j < 3; j++){
      for(k=0; k < 3; k++){
        for(ell=0; ell < 3; ell++){
          tmpJ[j+3*k+9*ell] = arrayJ[j][k][ell];
        }
      }
    }
    fscanf(fpb, "%lf %lf %lf %lf %lf %lf\n", &ReBx, &ImBx, &ReBy, &ImBy, &ReBz, &ImBz);
    vecB[0] = ReBx + I*ImBx;
    vecB[1] = ReBy + I*ImBy;
    vecB[2] = ReBz + I*ImBz;
  }
  fclose(fp);
  fclose(fpb);
  
#ifdef MPI
  MPI_Bcast(tmpJ, 27, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(vecB, 3, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
#endif
  for(j=0; j < 3; j++){
    for(k=0; k < 3; k++){
      for(ell=0; ell < 3; ell++){
        arrayJ[j][k][ell] = tmpJ[j+3*k+9*ell]; 
      }
    }
  }
//  printf(" ###Boost### Bcast arrayJ %lf\n",creal(arrayJ[0][0][0]));
  
  if(model_num==1){
// Kitaev on standard lattice N=W0*R0
// W0: 6 (fixed) 
// R0: R0>1
// num_pivot: 2 (fixed)
// ishift_nspin: 3 (fixed)
    
/* definition of model */
/*
  for(j=0; j < 3; j++){
    for(k=0; k < 3; k++){
      for(ell=0; ell < 3; ell++){
        arrayJ[j][k][ell] = 0.0;
      }
    }
  }
  arrayJ[0][0][0] = -1.0; //type=1 Jxx
  arrayJ[1][1][1] = -1.0; //type=2 Jyy
  arrayJ[2][2][2] = -1.0; //type=3 Jzz
*/
  for(j=0; j < 3; j++){
    for(k=0; k < 3; k++){
      for(ell=0; ell < 3; ell++){
        arrayJ[j][k][ell] = arrayJ[j][k][ell]*0.25;
      }
    }
    vecB[j] = vecB[j]*0.5;
  }

  for(j=0; j < (R0*num_pivot); j++){
    for(ell=0; ell < 7; ell++){
      for(k=0; k < 15; k++){
        list_6spin_pair[j][ell][k]=0;
      }
    } 
  }
  for(j=0; j < (R0*num_pivot); j++){
    for(ell=0; ell < 7; ell++){
      list_6spin_star[j][ell]=0;
    }
  }


  /* define list_6spin */ 
  for(j=0; j < R0; j++){

    list_6spin_star[2*j][0]=5; // num of J
    list_6spin_star[2*j][1]=1;
    list_6spin_star[2*j][2]=1;
    list_6spin_star[2*j][3]=1;
    list_6spin_star[2*j][4]=2;
    list_6spin_star[2*j][5]=1;
    list_6spin_star[2*j][6]=1; // flag

    list_6spin_pair[2*j][0][0]=0; //(1,1,1+2*j)=0 
    list_6spin_pair[2*j][1][0]=1; //(2,1,1+2*j)=1
    list_6spin_pair[2*j][2][0]=2; //(3,1,1+2*j)=2
    list_6spin_pair[2*j][3][0]=3; //(4,1,1+2*j)=3
    list_6spin_pair[2*j][4][0]=4; //(5,1,1+2*j)=4
    list_6spin_pair[2*j][5][0]=5; //(6,1,1+2*j)=5
    list_6spin_pair[2*j][6][0]=3; //(7,1,1+2*j)=3 ! type of J
    list_6spin_pair[2*j][0][1]=1; //(1,2,1+2*j)=1 
    list_6spin_pair[2*j][1][1]=2; //(2,2,1+2*j)=2
    list_6spin_pair[2*j][2][1]=0; //(3,2,1+2*j)=0
    list_6spin_pair[2*j][3][1]=3; //(4,2,1+2*j)=3
    list_6spin_pair[2*j][4][1]=4; //(5,2,1+2*j)=4
    list_6spin_pair[2*j][5][1]=5; //(6,2,1+2*j)=5
    list_6spin_pair[2*j][6][1]=1; //(7,2,1+2*j)=1 ! type of J
    list_6spin_pair[2*j][0][2]=2; //(1,3,1+2*j)=2 
    list_6spin_pair[2*j][1][2]=3; //(2,3,1+2*j)=3
    list_6spin_pair[2*j][2][2]=0; //(3,3,1+2*j)=0
    list_6spin_pair[2*j][3][2]=1; //(4,3,1+2*j)=1
    list_6spin_pair[2*j][4][2]=4; //(5,3,1+2*j)=4
    list_6spin_pair[2*j][5][2]=5; //(6,3,1+2*j)=5
    list_6spin_pair[2*j][6][2]=3; //(7,3,1+2*j)=3 ! type of J
    list_6spin_pair[2*j][0][3]=0; //(1,4,1+2*j)=0 
    list_6spin_pair[2*j][1][3]=4; //(2,4,1+2*j)=4
    list_6spin_pair[2*j][2][3]=1; //(3,4,1+2*j)=1
    list_6spin_pair[2*j][3][3]=2; //(4,4,1+2*j)=2
    list_6spin_pair[2*j][4][3]=3; //(5,4,1+2*j)=3
    list_6spin_pair[2*j][5][3]=5; //(6,4,1+2*j)=5
    list_6spin_pair[2*j][6][3]=1; //(7,4,1+2*j)=1 ! type of J
    list_6spin_pair[2*j][0][4]=1; //(1,5,1+2*j)=1 
    list_6spin_pair[2*j][1][4]=5; //(2,5,1+2*j)=5
    list_6spin_pair[2*j][2][4]=0; //(3,5,1+2*j)=0
    list_6spin_pair[2*j][3][4]=2; //(4,5,1+2*j)=2
    list_6spin_pair[2*j][4][4]=3; //(5,5,1+2*j)=3
    list_6spin_pair[2*j][5][4]=4; //(6,5,1+2*j)=4
    list_6spin_pair[2*j][6][4]=2; //(7,5,1+2*j)=2 ! type of J


    list_6spin_star[(2*j+1)][0]=4; //(0,2+2*j)=4 ! num of J
    list_6spin_star[(2*j+1)][1]=1; //(1,2+2*j)=1
    list_6spin_star[(2*j+1)][2]=1; //(2,2+2*j)=1
    list_6spin_star[(2*j+1)][3]=1; //(3,2+2*j)=1
    list_6spin_star[(2*j+1)][4]=2; //(4,2+2*j)=2
    list_6spin_star[(2*j+1)][5]=2; //(5,2+2*j)=2
    list_6spin_star[(2*j+1)][6]=1; //(6,2+2*j)=1 ! flag

    list_6spin_pair[(2*j+1)][0][0]=0; //(1,1,2+2*j)=0 
    list_6spin_pair[(2*j+1)][1][0]=1; //(2,1,2+2*j)=1
    list_6spin_pair[(2*j+1)][2][0]=2; //(3,1,2+2*j)=2
    list_6spin_pair[(2*j+1)][3][0]=3; //(4,1,2+2*j)=3
    list_6spin_pair[(2*j+1)][4][0]=4; //(5,1,2+2*j)=4
    list_6spin_pair[(2*j+1)][5][0]=5; //(6,1,2+2*j)=5
    list_6spin_pair[(2*j+1)][6][0]=1; //(7,1,2+2*j)=1 ! type of J
    list_6spin_pair[(2*j+1)][0][1]=1; //(1,2,2+2*j)=1 
    list_6spin_pair[(2*j+1)][1][1]=2; //(2,2,2+2*j)=2
    list_6spin_pair[(2*j+1)][2][1]=0; //(3,2,2+2*j)=0
    list_6spin_pair[(2*j+1)][3][1]=3; //(4,2,2+2*j)=3
    list_6spin_pair[(2*j+1)][4][1]=4; //(5,2,2+2*j)=4
    list_6spin_pair[(2*j+1)][5][1]=5; //(6,2,2+2*j)=5
    list_6spin_pair[(2*j+1)][6][1]=3; //(7,2,2+2*j)=3 ! type of J
    list_6spin_pair[(2*j+1)][0][2]=0; //(1,3,2+2*j)=0 
    list_6spin_pair[(2*j+1)][1][2]=4; //(2,3,2+2*j)=4
    list_6spin_pair[(2*j+1)][2][2]=1; //(3,3,2+2*j)=1
    list_6spin_pair[(2*j+1)][3][2]=2; //(4,3,2+2*j)=2
    list_6spin_pair[(2*j+1)][4][2]=3; //(5,3,2+2*j)=3
    list_6spin_pair[(2*j+1)][5][2]=5; //(6,3,2+2*j)=5
    list_6spin_pair[(2*j+1)][6][2]=2; //(7,3,2+2*j)=2 ! type of J
    list_6spin_pair[(2*j+1)][0][3]=2; //(1,4,2+2*j)=2 
    list_6spin_pair[(2*j+1)][1][3]=5; //(2,4,2+2*j)=5
    list_6spin_pair[(2*j+1)][2][3]=0; //(3,4,2+2*j)=0
    list_6spin_pair[(2*j+1)][3][3]=1; //(4,4,2+2*j)=1
    list_6spin_pair[(2*j+1)][4][3]=3; //(5,4,2+2*j)=3
    list_6spin_pair[(2*j+1)][5][3]=4; //(6,4,2+2*j)=4
    list_6spin_pair[(2*j+1)][6][3]=2; //(7,4,2+2*j)=2 ! type of J
  }/* define list_6spin */ 
/* definition of model */

  } else if(model_num==2){
// Kagome on standard lattice N=W0*R0
// W0: 9 (fixed) 
// R0: R0>1
// num_pivot: 4 (fixed)
// ishift_nspin: 3 (fixed)
    
/* definition of model */
/*
  for(j=0; j < 3; j++){
    for(k=0; k < 3; k++){
      for(ell=0; ell < 3; ell++){
        arrayJ[j][k][ell] = 0.0;
      }
    }
  }
  arrayJ[0][0][0] = 1.0;
  arrayJ[0][1][1] = 1.0;
  arrayJ[0][2][2] = 1.0;
*/
  for(j=0; j < 3; j++){
    for(k=0; k < 3; k++){
      for(ell=0; ell < 3; ell++){
        arrayJ[j][k][ell] = arrayJ[j][k][ell]*0.25;
      }
    }
    vecB[j] = vecB[j]*0.5;
  }
  for(j=0; j < (R0*num_pivot); j++){
    for(ell=0; ell < 7; ell++){
      for(k=0; k < 15; k++){
        list_6spin_pair[j][ell][k]=0;
      }
    } 
  }
  for(j=0; j < (R0*num_pivot); j++){
    for(ell=0; ell < 7; ell++){
      list_6spin_star[j][ell]=0;
    }
  }
  
  for(j=0; j<R0; j++){
    list_6spin_star[4*j][0]=1; // num of J
    list_6spin_star[4*j][1]=1;
    list_6spin_star[4*j][2]=1;
    list_6spin_star[4*j][3]=1;
    list_6spin_star[4*j][4]=4;
    list_6spin_star[4*j][5]=2;
    list_6spin_star[4*j][6]=-1; // flag
    list_6spin_pair[4*j][0][0]=0; //(1,1,1+2*j)=0 
    list_6spin_pair[4*j][1][0]=4; //(2,1,1+2*j)=1
    list_6spin_pair[4*j][2][0]=1; //(3,1,1+2*j)=2
    list_6spin_pair[4*j][3][0]=2; //(4,1,1+2*j)=3
    list_6spin_pair[4*j][4][0]=3; //(5,1,1+2*j)=4
    list_6spin_pair[4*j][5][0]=5; //(6,1,1+2*j)=5
    list_6spin_pair[4*j][6][0]=1; //(7,1,1+2*j)=3 ! type of J

    list_6spin_star[4*j+1][0]=6; // num of J
    list_6spin_star[4*j+1][1]=1;
    list_6spin_star[4*j+1][2]=1;
    list_6spin_star[4*j+1][3]=1;
    list_6spin_star[4*j+1][4]=6;
    list_6spin_star[4*j+1][5]=7;
    list_6spin_star[4*j+1][6]=1; // flag
    list_6spin_pair[4*j+1][0][0]=0;
    list_6spin_pair[4*j+1][1][0]=1;
    list_6spin_pair[4*j+1][2][0]=2;
    list_6spin_pair[4*j+1][3][0]=3;
    list_6spin_pair[4*j+1][4][0]=4;
    list_6spin_pair[4*j+1][5][0]=5;
    list_6spin_pair[4*j+1][6][0]=1; // type of J
    list_6spin_pair[4*j+1][0][1]=1;
    list_6spin_pair[4*j+1][1][1]=2;
    list_6spin_pair[4*j+1][2][1]=0;
    list_6spin_pair[4*j+1][3][1]=3;
    list_6spin_pair[4*j+1][4][1]=4;
    list_6spin_pair[4*j+1][5][1]=5;
    list_6spin_pair[4*j+1][6][1]=1; // type of J
    list_6spin_pair[4*j+1][0][2]=0;
    list_6spin_pair[4*j+1][1][2]=2;
    list_6spin_pair[4*j+1][2][2]=1;
    list_6spin_pair[4*j+1][3][2]=3;
    list_6spin_pair[4*j+1][4][2]=4;
    list_6spin_pair[4*j+1][5][2]=5;
    list_6spin_pair[4*j+1][6][2]=1; // type of J
    list_6spin_pair[4*j+1][0][3]=1;
    list_6spin_pair[4*j+1][1][3]=3;
    list_6spin_pair[4*j+1][2][3]=0;
    list_6spin_pair[4*j+1][3][3]=2;
    list_6spin_pair[4*j+1][4][3]=4;
    list_6spin_pair[4*j+1][5][3]=5;
    list_6spin_pair[4*j+1][6][3]=1; // type of J
    list_6spin_pair[4*j+1][0][4]=2;
    list_6spin_pair[4*j+1][1][4]=4;
    list_6spin_pair[4*j+1][2][4]=0;
    list_6spin_pair[4*j+1][3][4]=1;
    list_6spin_pair[4*j+1][4][4]=3;
    list_6spin_pair[4*j+1][5][4]=5;
    list_6spin_pair[4*j+1][6][4]=1; // type of J
    list_6spin_pair[4*j+1][0][5]=2;
    list_6spin_pair[4*j+1][1][5]=5;
    list_6spin_pair[4*j+1][2][5]=0;
    list_6spin_pair[4*j+1][3][5]=1;
    list_6spin_pair[4*j+1][4][5]=3;
    list_6spin_pair[4*j+1][5][5]=4;
    list_6spin_pair[4*j+1][6][5]=1; // type of J
    
    list_6spin_star[4*j+2][0]=6; // num of J
    list_6spin_star[4*j+2][1]=1;
    list_6spin_star[4*j+2][2]=1;
    list_6spin_star[4*j+2][3]=1;
    list_6spin_star[4*j+2][4]=4;
    list_6spin_star[4*j+2][5]=2;
    list_6spin_star[4*j+2][6]=1; // flag
    list_6spin_pair[4*j+2][0][0]=0;
    list_6spin_pair[4*j+2][1][0]=1;
    list_6spin_pair[4*j+2][2][0]=2;
    list_6spin_pair[4*j+2][3][0]=3;
    list_6spin_pair[4*j+2][4][0]=4;
    list_6spin_pair[4*j+2][5][0]=5;
    list_6spin_pair[4*j+2][6][0]=1; // type of J
    list_6spin_pair[4*j+2][0][1]=1;
    list_6spin_pair[4*j+2][1][1]=2;
    list_6spin_pair[4*j+2][2][1]=0;
    list_6spin_pair[4*j+2][3][1]=3;
    list_6spin_pair[4*j+2][4][1]=4;
    list_6spin_pair[4*j+2][5][1]=5;
    list_6spin_pair[4*j+2][6][1]=1; // type of J
    list_6spin_pair[4*j+2][0][2]=0;
    list_6spin_pair[4*j+2][1][2]=2;
    list_6spin_pair[4*j+2][2][2]=1;
    list_6spin_pair[4*j+2][3][2]=3;
    list_6spin_pair[4*j+2][4][2]=4;
    list_6spin_pair[4*j+2][5][2]=5;
    list_6spin_pair[4*j+2][6][2]=1; // type of J
    list_6spin_pair[4*j+2][0][3]=1;
    list_6spin_pair[4*j+2][1][3]=3;
    list_6spin_pair[4*j+2][2][3]=0;
    list_6spin_pair[4*j+2][3][3]=2;
    list_6spin_pair[4*j+2][4][3]=4;
    list_6spin_pair[4*j+2][5][3]=5;
    list_6spin_pair[4*j+2][6][3]=1; // type of J
    list_6spin_pair[4*j+2][0][4]=2;
    list_6spin_pair[4*j+2][1][4]=5;
    list_6spin_pair[4*j+2][2][4]=0;
    list_6spin_pair[4*j+2][3][4]=1;
    list_6spin_pair[4*j+2][4][4]=3;
    list_6spin_pair[4*j+2][5][4]=4;
    list_6spin_pair[4*j+2][6][4]=1; // type of J
    list_6spin_pair[4*j+2][0][5]=2;
    list_6spin_pair[4*j+2][1][5]=4;
    list_6spin_pair[4*j+2][2][5]=0;
    list_6spin_pair[4*j+2][3][5]=1;
    list_6spin_pair[4*j+2][4][5]=3;
    list_6spin_pair[4*j+2][5][5]=5;
    list_6spin_pair[4*j+2][6][5]=1; // type of J

    list_6spin_star[4*j+3][0]=5; // num of J
    list_6spin_star[4*j+3][1]=1;
    list_6spin_star[4*j+3][2]=1;
    list_6spin_star[4*j+3][3]=1;
    list_6spin_star[4*j+3][4]=4;
    list_6spin_star[4*j+3][5]=2;
    list_6spin_star[4*j+3][6]=1; // flag
    list_6spin_pair[4*j+3][0][0]=0;
    list_6spin_pair[4*j+3][1][0]=1;
    list_6spin_pair[4*j+3][2][0]=2;
    list_6spin_pair[4*j+3][3][0]=3;
    list_6spin_pair[4*j+3][4][0]=4;
    list_6spin_pair[4*j+3][5][0]=5;
    list_6spin_pair[4*j+3][6][0]=1; // type of J
    list_6spin_pair[4*j+3][0][1]=1;
    list_6spin_pair[4*j+3][1][1]=2;
    list_6spin_pair[4*j+3][2][1]=0;
    list_6spin_pair[4*j+3][3][1]=3;
    list_6spin_pair[4*j+3][4][1]=4;
    list_6spin_pair[4*j+3][5][1]=5;
    list_6spin_pair[4*j+3][6][1]=1; // type of J
    list_6spin_pair[4*j+3][0][2]=0;
    list_6spin_pair[4*j+3][1][2]=2;
    list_6spin_pair[4*j+3][2][2]=1;
    list_6spin_pair[4*j+3][3][2]=3;
    list_6spin_pair[4*j+3][4][2]=4;
    list_6spin_pair[4*j+3][5][2]=5;
    list_6spin_pair[4*j+3][6][2]=1; // type of J
    list_6spin_pair[4*j+3][0][3]=2;
    list_6spin_pair[4*j+3][1][3]=5;
    list_6spin_pair[4*j+3][2][3]=0;
    list_6spin_pair[4*j+3][3][3]=1;
    list_6spin_pair[4*j+3][4][3]=3;
    list_6spin_pair[4*j+3][5][3]=4;
    list_6spin_pair[4*j+3][6][3]=1; // type of J
    list_6spin_pair[4*j+3][0][4]=2;
    list_6spin_pair[4*j+3][1][4]=4;
    list_6spin_pair[4*j+3][2][4]=0;
    list_6spin_pair[4*j+3][3][4]=1;
    list_6spin_pair[4*j+3][4][4]=3;
    list_6spin_pair[4*j+3][5][4]=5;
    list_6spin_pair[4*j+3][6][4]=1; // type of J
  }
/* definition of model */

  } else if(model_num==3){
// Kagome on a 36 site cluster (cluster (E) in Nakano-Sakai JPSJ 79, 053707 (2010))     
// W0: 18 (fixed) 
// R0: 2 (fixed)
// num_pivot: 8 (fixed)
// ishift_nspin: 3 (fixed)

/* definition of model */
/*
  for(j=0; j < 3; j++){
    for(k=0; k < 3; k++){
      for(ell=0; ell < 3; ell++){
        arrayJ[j][k][ell] = 0.0;
      }
    }
  }
  arrayJ[0][0][0] = 1.0;
  arrayJ[0][1][1] = 1.0;
  arrayJ[0][2][2] = 1.0;
*/
  for(j=0; j < 3; j++){
    for(k=0; k < 3; k++){
      for(ell=0; ell < 3; ell++){
        arrayJ[j][k][ell] = arrayJ[j][k][ell]*0.25;
      }
    }
    vecB[j] = vecB[j]*0.5;
  }
  for(j=0; j < (R0*num_pivot); j++){
    for(ell=0; ell < 7; ell++){
      for(k=0; k < 15; k++){
        list_6spin_pair[j][ell][k]=0;
      }
    } 
  }
  for(j=0; j < (R0*num_pivot); j++){
    for(ell=0; ell < 7; ell++){
      list_6spin_star[j][ell]=0;
    }
  }
  
  for(j=0; j<R0; j++){
    list_6spin_star[8*j][0]=6; // num of J
    list_6spin_star[8*j][1]=1;
    list_6spin_star[8*j][2]=1;
    list_6spin_star[8*j][3]=1;
    list_6spin_star[8*j][4]=6;
    list_6spin_star[8*j][5]=8;
    list_6spin_star[8*j][6]=1; // flag
    list_6spin_pair[8*j][0][0]=0; 
    list_6spin_pair[8*j][1][0]=1;
    list_6spin_pair[8*j][2][0]=2;
    list_6spin_pair[8*j][3][0]=3;
    list_6spin_pair[8*j][4][0]=4;
    list_6spin_pair[8*j][5][0]=5;
    list_6spin_pair[8*j][6][0]=1; // type of J
    list_6spin_pair[8*j][0][1]=1; 
    list_6spin_pair[8*j][1][1]=2;
    list_6spin_pair[8*j][2][1]=0;
    list_6spin_pair[8*j][3][1]=3;
    list_6spin_pair[8*j][4][1]=4;
    list_6spin_pair[8*j][5][1]=5;
    list_6spin_pair[8*j][6][1]=1; // type of J
    list_6spin_pair[8*j][0][2]=0; 
    list_6spin_pair[8*j][1][2]=2;
    list_6spin_pair[8*j][2][2]=1;
    list_6spin_pair[8*j][3][2]=3;
    list_6spin_pair[8*j][4][2]=4;
    list_6spin_pair[8*j][5][2]=5;
    list_6spin_pair[8*j][6][2]=1; // type of J
    list_6spin_pair[8*j][0][3]=1; 
    list_6spin_pair[8*j][1][3]=3;
    list_6spin_pair[8*j][2][3]=0;
    list_6spin_pair[8*j][3][3]=2;
    list_6spin_pair[8*j][4][3]=4;
    list_6spin_pair[8*j][5][3]=5;
    list_6spin_pair[8*j][6][3]=1; // type of J
    list_6spin_pair[8*j][0][4]=2; 
    list_6spin_pair[8*j][1][4]=4;
    list_6spin_pair[8*j][2][4]=0;
    list_6spin_pair[8*j][3][4]=1;
    list_6spin_pair[8*j][4][4]=3;
    list_6spin_pair[8*j][5][4]=5;
    list_6spin_pair[8*j][6][4]=1; // type of J
    list_6spin_pair[8*j][0][5]=0; 
    list_6spin_pair[8*j][1][5]=5;
    list_6spin_pair[8*j][2][5]=1;
    list_6spin_pair[8*j][3][5]=2;
    list_6spin_pair[8*j][4][5]=3;
    list_6spin_pair[8*j][5][5]=4;
    list_6spin_pair[8*j][6][5]=1; // type of J

    list_6spin_star[8*j+1][0]=6; // num of J
    list_6spin_star[8*j+1][1]=1;
    list_6spin_star[8*j+1][2]=1;
    list_6spin_star[8*j+1][3]=1;
    list_6spin_star[8*j+1][4]=4;
    list_6spin_star[8*j+1][5]=2;
    list_6spin_star[8*j+1][6]=1; // flag
    list_6spin_pair[8*j+1][0][0]=0; 
    list_6spin_pair[8*j+1][1][0]=1;
    list_6spin_pair[8*j+1][2][0]=2;
    list_6spin_pair[8*j+1][3][0]=3;
    list_6spin_pair[8*j+1][4][0]=4;
    list_6spin_pair[8*j+1][5][0]=5;
    list_6spin_pair[8*j+1][6][0]=1; // type of J
    list_6spin_pair[8*j+1][0][1]=1; 
    list_6spin_pair[8*j+1][1][1]=2;
    list_6spin_pair[8*j+1][2][1]=0;
    list_6spin_pair[8*j+1][3][1]=3;
    list_6spin_pair[8*j+1][4][1]=4;
    list_6spin_pair[8*j+1][5][1]=5;
    list_6spin_pair[8*j+1][6][1]=1; // type of J
    list_6spin_pair[8*j+1][0][2]=0; 
    list_6spin_pair[8*j+1][1][2]=2;
    list_6spin_pair[8*j+1][2][2]=1;
    list_6spin_pair[8*j+1][3][2]=3;
    list_6spin_pair[8*j+1][4][2]=4;
    list_6spin_pair[8*j+1][5][2]=5;
    list_6spin_pair[8*j+1][6][2]=1; // type of J
    list_6spin_pair[8*j+1][0][3]=1; 
    list_6spin_pair[8*j+1][1][3]=3;
    list_6spin_pair[8*j+1][2][3]=0;
    list_6spin_pair[8*j+1][3][3]=2;
    list_6spin_pair[8*j+1][4][3]=4;
    list_6spin_pair[8*j+1][5][3]=5;
    list_6spin_pair[8*j+1][6][3]=1; // type of J
    list_6spin_pair[8*j+1][0][4]=2; 
    list_6spin_pair[8*j+1][1][4]=4;
    list_6spin_pair[8*j+1][2][4]=0;
    list_6spin_pair[8*j+1][3][4]=1;
    list_6spin_pair[8*j+1][4][4]=3;
    list_6spin_pair[8*j+1][5][4]=5;
    list_6spin_pair[8*j+1][6][4]=1; // type of J
    list_6spin_pair[8*j+1][0][5]=2; 
    list_6spin_pair[8*j+1][1][5]=5;
    list_6spin_pair[8*j+1][2][5]=0;
    list_6spin_pair[8*j+1][3][5]=1;
    list_6spin_pair[8*j+1][4][5]=3;
    list_6spin_pair[8*j+1][5][5]=4;
    list_6spin_pair[8*j+1][6][5]=1; // type of J

    list_6spin_star[8*j+2][0]=5; // num of J
    list_6spin_star[8*j+2][1]=1;
    list_6spin_star[8*j+2][2]=1;
    list_6spin_star[8*j+2][3]=1;
    list_6spin_star[8*j+2][4]=4;
    list_6spin_star[8*j+2][5]=2;
    list_6spin_star[8*j+2][6]=-1; // flag
    list_6spin_pair[8*j+2][0][0]=0; 
    list_6spin_pair[8*j+2][1][0]=1;
    list_6spin_pair[8*j+2][2][0]=2;
    list_6spin_pair[8*j+2][3][0]=3;
    list_6spin_pair[8*j+2][4][0]=4;
    list_6spin_pair[8*j+2][5][0]=5;
    list_6spin_pair[8*j+2][6][0]=1; // type of J
    list_6spin_pair[8*j+2][0][1]=1; 
    list_6spin_pair[8*j+2][1][1]=2;
    list_6spin_pair[8*j+2][2][1]=0;
    list_6spin_pair[8*j+2][3][1]=3;
    list_6spin_pair[8*j+2][4][1]=4;
    list_6spin_pair[8*j+2][5][1]=5;
    list_6spin_pair[8*j+2][6][1]=1; // type of J
    list_6spin_pair[8*j+2][0][2]=0; 
    list_6spin_pair[8*j+2][1][2]=2;
    list_6spin_pair[8*j+2][2][2]=1;
    list_6spin_pair[8*j+2][3][2]=3;
    list_6spin_pair[8*j+2][4][2]=4;
    list_6spin_pair[8*j+2][5][2]=5;
    list_6spin_pair[8*j+2][6][2]=1; // type of J
    list_6spin_pair[8*j+2][0][3]=2; 
    list_6spin_pair[8*j+2][1][3]=4;
    list_6spin_pair[8*j+2][2][3]=0;
    list_6spin_pair[8*j+2][3][3]=1;
    list_6spin_pair[8*j+2][4][3]=3;
    list_6spin_pair[8*j+2][5][3]=5;
    list_6spin_pair[8*j+2][6][3]=1; // type of J
    list_6spin_pair[8*j+2][0][4]=2; 
    list_6spin_pair[8*j+2][1][4]=5;
    list_6spin_pair[8*j+2][2][4]=0;
    list_6spin_pair[8*j+2][3][4]=1;
    list_6spin_pair[8*j+2][4][4]=3;
    list_6spin_pair[8*j+2][5][4]=4;
    list_6spin_pair[8*j+2][6][4]=1; // type of J

    list_6spin_star[8*j+3][0]=1; // num of J
    list_6spin_star[8*j+3][1]=1;
    list_6spin_star[8*j+3][2]=1;
    list_6spin_star[8*j+3][3]=1;
    list_6spin_star[8*j+3][4]=6;
    list_6spin_star[8*j+3][5]=3;
    list_6spin_star[8*j+3][6]=1; // flag
    list_6spin_pair[8*j+3][0][0]=1; 
    list_6spin_pair[8*j+3][1][0]=5;
    list_6spin_pair[8*j+3][2][0]=0;
    list_6spin_pair[8*j+3][3][0]=2;
    list_6spin_pair[8*j+3][4][0]=3;
    list_6spin_pair[8*j+3][5][0]=4;
    list_6spin_pair[8*j+3][6][0]=1; // type of J


    list_6spin_star[8*j+4][0]=6; // num of J
    list_6spin_star[8*j+4][1]=1;
    list_6spin_star[8*j+4][2]=1;
    list_6spin_star[8*j+4][3]=1;
    list_6spin_star[8*j+4][4]=7;
    list_6spin_star[8*j+4][5]=2;
    list_6spin_star[8*j+4][6]=1; // flag
    list_6spin_pair[8*j+4][0][0]=0; 
    list_6spin_pair[8*j+4][1][0]=1;
    list_6spin_pair[8*j+4][2][0]=2;
    list_6spin_pair[8*j+4][3][0]=3;
    list_6spin_pair[8*j+4][4][0]=4;
    list_6spin_pair[8*j+4][5][0]=5;
    list_6spin_pair[8*j+4][6][0]=1; // type of J
    list_6spin_pair[8*j+4][0][1]=1; 
    list_6spin_pair[8*j+4][1][1]=2;
    list_6spin_pair[8*j+4][2][1]=0;
    list_6spin_pair[8*j+4][3][1]=3;
    list_6spin_pair[8*j+4][4][1]=4;
    list_6spin_pair[8*j+4][5][1]=5;
    list_6spin_pair[8*j+4][6][1]=1; // type of J
    list_6spin_pair[8*j+4][0][2]=0; 
    list_6spin_pair[8*j+4][1][2]=2;
    list_6spin_pair[8*j+4][2][2]=1;
    list_6spin_pair[8*j+4][3][2]=3;
    list_6spin_pair[8*j+4][4][2]=4;
    list_6spin_pair[8*j+4][5][2]=5;
    list_6spin_pair[8*j+4][6][2]=1; // type of J
    list_6spin_pair[8*j+4][0][3]=1; 
    list_6spin_pair[8*j+4][1][3]=3;
    list_6spin_pair[8*j+4][2][3]=0;
    list_6spin_pair[8*j+4][3][3]=2;
    list_6spin_pair[8*j+4][4][3]=4;
    list_6spin_pair[8*j+4][5][3]=5;
    list_6spin_pair[8*j+4][6][3]=1; // type of J
    list_6spin_pair[8*j+4][0][4]=2; 
    list_6spin_pair[8*j+4][1][4]=4;
    list_6spin_pair[8*j+4][2][4]=0;
    list_6spin_pair[8*j+4][3][4]=1;
    list_6spin_pair[8*j+4][4][4]=3;
    list_6spin_pair[8*j+4][5][4]=5;
    list_6spin_pair[8*j+4][6][4]=1; // type of J
    list_6spin_pair[8*j+4][0][5]=2; 
    list_6spin_pair[8*j+4][1][5]=5;
    list_6spin_pair[8*j+4][2][5]=0;
    list_6spin_pair[8*j+4][3][5]=1;
    list_6spin_pair[8*j+4][4][5]=3;
    list_6spin_pair[8*j+4][5][5]=4;
    list_6spin_pair[8*j+4][6][5]=1; // type of J

    list_6spin_star[8*j+5][0]=6; // num of J
    list_6spin_star[8*j+5][1]=1;
    list_6spin_star[8*j+5][2]=1;
    list_6spin_star[8*j+5][3]=1;
    list_6spin_star[8*j+5][4]=7;
    list_6spin_star[8*j+5][5]=2;
    list_6spin_star[8*j+5][6]=1; // flag
    list_6spin_pair[8*j+5][0][0]=0; 
    list_6spin_pair[8*j+5][1][0]=1;
    list_6spin_pair[8*j+5][2][0]=2;
    list_6spin_pair[8*j+5][3][0]=3;
    list_6spin_pair[8*j+5][4][0]=4;
    list_6spin_pair[8*j+5][5][0]=5;
    list_6spin_pair[8*j+5][6][0]=1; // type of J
    list_6spin_pair[8*j+5][0][1]=1; 
    list_6spin_pair[8*j+5][1][1]=2;
    list_6spin_pair[8*j+5][2][1]=0;
    list_6spin_pair[8*j+5][3][1]=3;
    list_6spin_pair[8*j+5][4][1]=4;
    list_6spin_pair[8*j+5][5][1]=5;
    list_6spin_pair[8*j+5][6][1]=1; // type of J
    list_6spin_pair[8*j+5][0][2]=0; 
    list_6spin_pair[8*j+5][1][2]=2;
    list_6spin_pair[8*j+5][2][2]=1;
    list_6spin_pair[8*j+5][3][2]=3;
    list_6spin_pair[8*j+5][4][2]=4;
    list_6spin_pair[8*j+5][5][2]=5;
    list_6spin_pair[8*j+5][6][2]=1; // type of J
    list_6spin_pair[8*j+5][0][3]=1; 
    list_6spin_pair[8*j+5][1][3]=3;
    list_6spin_pair[8*j+5][2][3]=0;
    list_6spin_pair[8*j+5][3][3]=2;
    list_6spin_pair[8*j+5][4][3]=4;
    list_6spin_pair[8*j+5][5][3]=5;
    list_6spin_pair[8*j+5][6][3]=1; // type of J
    list_6spin_pair[8*j+5][0][4]=2; 
    list_6spin_pair[8*j+5][1][4]=4;
    list_6spin_pair[8*j+5][2][4]=0;
    list_6spin_pair[8*j+5][3][4]=1;
    list_6spin_pair[8*j+5][4][4]=3;
    list_6spin_pair[8*j+5][5][4]=5;
    list_6spin_pair[8*j+5][6][4]=1; // type of J
    list_6spin_pair[8*j+5][0][5]=2; 
    list_6spin_pair[8*j+5][1][5]=5;
    list_6spin_pair[8*j+5][2][5]=0;
    list_6spin_pair[8*j+5][3][5]=1;
    list_6spin_pair[8*j+5][4][5]=3;
    list_6spin_pair[8*j+5][5][5]=4;
    list_6spin_pair[8*j+5][6][5]=1; // type of J

    list_6spin_star[8*j+6][0]=5; // num of J
    list_6spin_star[8*j+6][1]=1;
    list_6spin_star[8*j+6][2]=1;
    list_6spin_star[8*j+6][3]=1;
    list_6spin_star[8*j+6][4]=2;
    list_6spin_star[8*j+6][5]=5;
    list_6spin_star[8*j+6][6]=-1; // flag
    list_6spin_pair[8*j+6][0][0]=0; 
    list_6spin_pair[8*j+6][1][0]=1;
    list_6spin_pair[8*j+6][2][0]=2;
    list_6spin_pair[8*j+6][3][0]=3;
    list_6spin_pair[8*j+6][4][0]=4;
    list_6spin_pair[8*j+6][5][0]=5;
    list_6spin_pair[8*j+6][6][0]=1; // type of J
    list_6spin_pair[8*j+6][0][1]=1; 
    list_6spin_pair[8*j+6][1][1]=2;
    list_6spin_pair[8*j+6][2][1]=0;
    list_6spin_pair[8*j+6][3][1]=3;
    list_6spin_pair[8*j+6][4][1]=4;
    list_6spin_pair[8*j+6][5][1]=5;
    list_6spin_pair[8*j+6][6][1]=1; // type of J
    list_6spin_pair[8*j+6][0][2]=0; 
    list_6spin_pair[8*j+6][1][2]=2;
    list_6spin_pair[8*j+6][2][2]=1;
    list_6spin_pair[8*j+6][3][2]=3;
    list_6spin_pair[8*j+6][4][2]=4;
    list_6spin_pair[8*j+6][5][2]=5;
    list_6spin_pair[8*j+6][6][2]=1; // type of J
    list_6spin_pair[8*j+6][0][3]=1; 
    list_6spin_pair[8*j+6][1][3]=4;
    list_6spin_pair[8*j+6][2][3]=0;
    list_6spin_pair[8*j+6][3][3]=2;
    list_6spin_pair[8*j+6][4][3]=3;
    list_6spin_pair[8*j+6][5][3]=5;
    list_6spin_pair[8*j+6][6][3]=1; // type of J
    list_6spin_pair[8*j+6][0][4]=2; 
    list_6spin_pair[8*j+6][1][4]=5;
    list_6spin_pair[8*j+6][2][4]=0;
    list_6spin_pair[8*j+6][3][4]=1;
    list_6spin_pair[8*j+6][4][4]=3;
    list_6spin_pair[8*j+6][5][4]=4;
    list_6spin_pair[8*j+6][6][4]=1; // type of J

    list_6spin_star[8*j+7][0]=1; // num of J
    list_6spin_star[8*j+7][1]=1;
    list_6spin_star[8*j+7][2]=1;
    list_6spin_star[8*j+7][3]=1;
    list_6spin_star[8*j+7][4]=2;
    list_6spin_star[8*j+7][5]=7;
    list_6spin_star[8*j+7][6]=1; // flag
    list_6spin_pair[8*j+7][0][0]=1; 
    list_6spin_pair[8*j+7][1][0]=5;
    list_6spin_pair[8*j+7][2][0]=0;
    list_6spin_pair[8*j+7][3][0]=2;
    list_6spin_pair[8*j+7][4][0]=3;
    list_6spin_pair[8*j+7][5][0]=4;
    list_6spin_pair[8*j+7][6][0]=1; // type of J
  }
/* definition of model */

  } 

  c_free1(tmpJ, 27);
  /*
    Print Input file for the expartmode + Boost
  */
  fp = fopenMPI("Boost.def", "w");
  fprintf(fp, "%d\n", 1);/**/

}
