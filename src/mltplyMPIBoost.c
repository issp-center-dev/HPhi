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

//Define Mode for mltply
// complex version

#ifdef MPI
#include "mpi.h"
#endif
#include "Common.h"
#include "mfmemory.h"
#include "xsetmem.h"
#include "mltply.h"
#include "bitcalc.h"
#include "wrapperMPI.h"
#include "mltplyMPI.h"


/**
 *
 * Exchange term in Spin model
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 * @author Youhei Yamaji (The University of Tokyo)
 */
void child_general_int_spin_MPIBoost(
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/,
  double complex *tmp_v2 /**< [inout] bufffer*/,
  double complex *tmp_v3 /**< [inout] bufffer*/
  )
{
#ifdef MPI
  
  double complex dam_pr = 0;
  MPI_Status statusMPI;

/* read def for star and pair  */
  char *filename1 = "starBoost";
  char *filename2 = "pairBoost";
  FILE *fp1;
  FILE *fp2;

  int ierr;
  int INFO;
  char TRANSA, TRANSB;
  int M, N, K, LDA, LDB, LDC;
  double complex ALPHA, BETA;

  long unsigned int i_max;
  long unsigned int j, k, ell, iloop;
  long unsigned int i1, i2;
  long unsigned int iomp;
  long unsigned int ell4, ell5, ell6, m0, Ipart1;
  long unsigned int W0, R0, num_pivot, ishift_nspin;
  long unsigned int mi, mj, mri, mrj, mrk, mrl, indj;
  long unsigned int ellrl, ellrk, ellrj, ellri, elli1, elli2, ellj1, ellj2;
  long unsigned int iSS1, iSS2, iSSL1, iSSL2;
  double complex ***arrayJ;
  double complex **vecJ;
  double complex **matJ, **matJ2;
  double complex *matJL;
  double complex *arrayz, *arrayx;
  int  **list_6spin_star;
  int ***list_6spin_pair;
  long unsigned int ishift1, ishift2, ishift3, ishift4, ishift5, pivot_flag, num_J_star; 
  long unsigned int pow1, pow2, pow3, pow4, pow5, pow11, pow21, pow31, pow41, pow51; 

  i_max = X->Check.idim_max;

  if((fp1 = fopen(filename1, "r")) == NULL){
    fprintf(stderr, "\n ###Boost### failed to open a file %s\n", filename1);
    exit(EXIT_FAILURE);
  }
  
  fscanf(fp1, "%d %d %d %d\n", &W0, &R0, &num_pivot, &ishift_nspin);

  fclose(fp1);
  
  if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode subroutine %d %d %d %d \n\n", W0, R0, num_pivot, ishift_nspin);}

  c_malloc3(arrayJ, 3, 3, 3); 
  c_malloc2(vecJ, 3, 3); 
  c_malloc2(matJ, 4, 4); 
  c_malloc2(matJ2, 4, 4); 
  c_malloc1(matJL, (int)(64*64)); 
  i_malloc2(list_6spin_star, (int)(R0*num_pivot), 7); 
  i_malloc3(list_6spin_pair, (int)(R0*num_pivot), 7, 21); 

  if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode subroutine list allocated \n\n");}

  for(j=0; j < 3; j++){
    for(k=0; k < 3; k++){
      for(ell=0; ell < 3; ell++){
        arrayJ[j][k][ell] = 0.0;
      }
    }
  }
  if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode subroutine J zero clear \n\n");}
  arrayJ[0][0][0] = -1.0; //type=1 Jxx
  arrayJ[1][1][1] = -1.0; //type=2 Jyy
  arrayJ[2][2][2] = -1.0; //type=3 Jzz
  for(j=0; j < 3; j++){
    for(k=0; k < 3; k++){
      for(ell=0; ell < 3; ell++){
        arrayJ[j][k][ell] *= 0.25;
      }
    }
  }
  if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode subroutine J multiplied 1/4 \n\n");}

  for(j=0; j < (int)(R0*num_pivot); j++){
    for(ell=0; ell < 7; ell++){
      for(k=0; k < 21; k++){
        list_6spin_pair[j][ell][k]=0;
      }
    } 
  }
  for(j=0; j < (int)(R0*num_pivot); j++){
    for(ell=0; ell < 7; ell++){
      list_6spin_star[j][ell]=0;
    }
  }

  if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode subroutine list set \n\n");}

  /* define list_6spin */ 
  for(j=0; j < R0; j++){

    
    list_6spin_star[(int)2*j][0]=5; // num of J
    list_6spin_star[(int)2*j][1]=1;
    list_6spin_star[(int)2*j][2]=1;
    list_6spin_star[(int)2*j][3]=1;
    list_6spin_star[(int)2*j][4]=2;
    list_6spin_star[(int)2*j][5]=1;
    list_6spin_star[(int)2*j][6]=1; // flag

    list_6spin_pair[(int)2*j][0][0]=0; //(1,1,1+2*j)=0 
    list_6spin_pair[(int)2*j][1][0]=1; //(2,1,1+2*j)=1
    list_6spin_pair[(int)2*j][2][0]=2; //(3,1,1+2*j)=2
    list_6spin_pair[(int)2*j][3][0]=3; //(4,1,1+2*j)=3
    list_6spin_pair[(int)2*j][4][0]=4; //(5,1,1+2*j)=4
    list_6spin_pair[(int)2*j][5][0]=5; //(6,1,1+2*j)=5
    list_6spin_pair[(int)2*j][6][0]=3; //(7,1,1+2*j)=3 ! type of J
    list_6spin_pair[(int)2*j][0][1]=1; //(1,2,1+2*j)=1 
    list_6spin_pair[(int)2*j][1][1]=2; //(2,2,1+2*j)=2
    list_6spin_pair[(int)2*j][2][1]=0; //(3,2,1+2*j)=0
    list_6spin_pair[(int)2*j][3][1]=3; //(4,2,1+2*j)=3
    list_6spin_pair[(int)2*j][4][1]=4; //(5,2,1+2*j)=4
    list_6spin_pair[(int)2*j][5][1]=5; //(6,2,1+2*j)=5
    list_6spin_pair[(int)2*j][6][1]=1; //(7,2,1+2*j)=1 ! type of J
    list_6spin_pair[(int)2*j][0][2]=2; //(1,3,1+2*j)=2 
    list_6spin_pair[(int)2*j][1][2]=3; //(2,3,1+2*j)=3
    list_6spin_pair[(int)2*j][2][2]=0; //(3,3,1+2*j)=0
    list_6spin_pair[(int)2*j][3][2]=1; //(4,3,1+2*j)=1
    list_6spin_pair[(int)2*j][4][2]=4; //(5,3,1+2*j)=4
    list_6spin_pair[(int)2*j][5][2]=5; //(6,3,1+2*j)=5
    list_6spin_pair[(int)2*j][6][2]=3; //(7,3,1+2*j)=3 ! type of J
    list_6spin_pair[(int)2*j][0][3]=0; //(1,4,1+2*j)=0 
    list_6spin_pair[(int)2*j][1][3]=4; //(2,4,1+2*j)=4
    list_6spin_pair[(int)2*j][2][3]=1; //(3,4,1+2*j)=1
    list_6spin_pair[(int)2*j][3][3]=2; //(4,4,1+2*j)=2
    list_6spin_pair[(int)2*j][4][3]=3; //(5,4,1+2*j)=3
    list_6spin_pair[(int)2*j][5][3]=5; //(6,4,1+2*j)=5
    list_6spin_pair[(int)2*j][6][3]=1; //(7,4,1+2*j)=1 ! type of J
    list_6spin_pair[(int)2*j][0][4]=1; //(1,5,1+2*j)=1 
    list_6spin_pair[(int)2*j][1][4]=5; //(2,5,1+2*j)=5
    list_6spin_pair[(int)2*j][2][4]=0; //(3,5,1+2*j)=0
    list_6spin_pair[(int)2*j][3][4]=2; //(4,5,1+2*j)=2
    list_6spin_pair[(int)2*j][4][4]=3; //(5,5,1+2*j)=3
    list_6spin_pair[(int)2*j][5][4]=4; //(6,5,1+2*j)=4
    list_6spin_pair[(int)2*j][6][4]=2; //(7,5,1+2*j)=2 ! type of J


    list_6spin_star[(int)(2*j+1)][0]=4; //(0,2+2*j)=4 ! num of J
    list_6spin_star[(int)(2*j+1)][1]=1; //(1,2+2*j)=1
    list_6spin_star[(int)(2*j+1)][2]=1; //(2,2+2*j)=1
    list_6spin_star[(int)(2*j+1)][3]=1; //(3,2+2*j)=1
    list_6spin_star[(int)(2*j+1)][4]=2; //(4,2+2*j)=2
    list_6spin_star[(int)(2*j+1)][5]=2; //(5,2+2*j)=2
    list_6spin_star[(int)(2*j+1)][6]=1; //(6,2+2*j)=1 ! flag

    list_6spin_pair[(int)(2*j+1)][0][0]=0; //(1,1,2+2*j)=0 
    list_6spin_pair[(int)(2*j+1)][1][0]=1; //(2,1,2+2*j)=1
    list_6spin_pair[(int)(2*j+1)][2][0]=2; //(3,1,2+2*j)=2
    list_6spin_pair[(int)(2*j+1)][3][0]=3; //(4,1,2+2*j)=3
    list_6spin_pair[(int)(2*j+1)][4][0]=4; //(5,1,2+2*j)=4
    list_6spin_pair[(int)(2*j+1)][5][0]=5; //(6,1,2+2*j)=5
    list_6spin_pair[(int)(2*j+1)][6][0]=1; //(7,1,2+2*j)=1 ! type of J
    list_6spin_pair[(int)(2*j+1)][0][1]=1; //(1,2,2+2*j)=1 
    list_6spin_pair[(int)(2*j+1)][1][1]=2; //(2,2,2+2*j)=2
    list_6spin_pair[(int)(2*j+1)][2][1]=0; //(3,2,2+2*j)=0
    list_6spin_pair[(int)(2*j+1)][3][1]=3; //(4,2,2+2*j)=3
    list_6spin_pair[(int)(2*j+1)][4][1]=4; //(5,2,2+2*j)=4
    list_6spin_pair[(int)(2*j+1)][5][1]=5; //(6,2,2+2*j)=5
    list_6spin_pair[(int)(2*j+1)][6][1]=3; //(7,2,2+2*j)=3 ! type of J
    list_6spin_pair[(int)(2*j+1)][0][2]=0; //(1,3,2+2*j)=0 
    list_6spin_pair[(int)(2*j+1)][1][2]=4; //(2,3,2+2*j)=4
    list_6spin_pair[(int)(2*j+1)][2][2]=1; //(3,3,2+2*j)=1
    list_6spin_pair[(int)(2*j+1)][3][2]=2; //(4,3,2+2*j)=2
    list_6spin_pair[(int)(2*j+1)][4][2]=3; //(5,3,2+2*j)=3
    list_6spin_pair[(int)(2*j+1)][5][2]=5; //(6,3,2+2*j)=5
    list_6spin_pair[(int)(2*j+1)][6][2]=2; //(7,3,2+2*j)=2 ! type of J
    list_6spin_pair[(int)(2*j+1)][0][3]=2; //(1,4,2+2*j)=2 
    list_6spin_pair[(int)(2*j+1)][1][3]=5; //(2,4,2+2*j)=5
    list_6spin_pair[(int)(2*j+1)][2][3]=0; //(3,4,2+2*j)=0
    list_6spin_pair[(int)(2*j+1)][3][3]=1; //(4,4,2+2*j)=1
    list_6spin_pair[(int)(2*j+1)][4][3]=3; //(5,4,2+2*j)=3
    list_6spin_pair[(int)(2*j+1)][5][3]=4; //(6,4,2+2*j)=4
    list_6spin_pair[(int)(2*j+1)][6][3]=2; //(7,4,2+2*j)=2 ! type of J
  }/* define list_6spin */ 


  for(iloop=0; iloop < R0; iloop++){

    if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode subroutine %d th column\n\n",iloop);}

    for(j=iloop*num_pivot; j < (iloop+1)*num_pivot; j++){
      
      if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode subroutine %d th pivot\n\n",j);}
       
      num_J_star = list_6spin_star[j][0]; //(0,j) 
      ishift1    = list_6spin_star[j][1]; //(1,j) 
      ishift2    = list_6spin_star[j][2]; //(2,j) 
      ishift3    = list_6spin_star[j][3]; //(3,j)
      ishift4    = list_6spin_star[j][4]; //(4,j)
      ishift5    = list_6spin_star[j][5]; //(5,j)
      pivot_flag = list_6spin_star[j][6]; //(6,j)
      pow1 = pow(2,ishift1);
      pow2 = pow(2,ishift1+ishift2);
      pow3 = pow(2,ishift1+ishift2+ishift3);
      pow4 = pow(2,ishift1+ishift2+ishift3+ishift4);
      pow5 = pow(2,ishift1+ishift2+ishift3+ishift4+ishift5);
      pow11= pow(2,ishift1+1);
      pow21= pow(2,ishift1+ishift2+1);
      pow31= pow(2,ishift1+ishift2+ishift3+1);
      pow41= pow(2,ishift1+ishift2+ishift3+ishift4+1);
      pow51= pow(2,ishift1+ishift2+ishift3+ishift4+ishift5+1);

      for(k=0; k < 64*64; k++){
        matJL[k] = 0.0 + 0.0*I;
      }

      for(ell=0; ell < num_J_star; ell++){
        mi   = list_6spin_pair[j][0][ell]; //(1,ell,j)
        mj   = list_6spin_pair[j][1][ell]; //(2,ell,j)
        mri  = list_6spin_pair[j][2][ell]; //(3,ell,j)
        mrj  = list_6spin_pair[j][3][ell]; //(4,ell,j)
        mrk  = list_6spin_pair[j][4][ell]; //(5,ell,j)
        mrl  = list_6spin_pair[j][5][ell]; //(6,ell,j)
        indj = list_6spin_pair[j][6][ell]; //(6,ell,j)
        for(i1 = 0; i1 < 3; i1++){
          for(i2 = 0; i2 < 3; i2++){
            vecJ[i1][i2] = arrayJ[indj-1][i1][i2];
          }
        } 
        if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode subroutine vecJ set %d\n\n",ell);}
        //matJSS(1,1) = vecJ(3,3)
        matJ[0][0] = vecJ[2][2];
        //matJSS(1,2)= vecJ(1,1)-vecJ(2,2)-dcmplx(0.0d0,1.0d0)*vecJ(1,2)-dcmplx(0.0d0,1.0d0)*vecJ(2,1)
        matJ[0][1] = vecJ[0][0]-vecJ[1][1]-I*vecJ[0][1]-I*vecJ[1][0];
        //matJSS(1,3)= vecJ(3,1)-dcmplx(0.0d0,1.0d0)*vecJ(3,2)
        matJ[0][2] = vecJ[2][0]-I*vecJ[2][1];
        //matJSS(1,4)= vecJ(1,3)-dcmplx(0.0d0,1.0d0)*vecJ(2,3)
        matJ[0][3] = vecJ[0][2]-I*vecJ[1][2];
        //matJSS(2,1)= vecJ(1,1)-vecJ(2,2)+dcmplx(0.0d0,1.0d0)*vecJ(1,2)+dcmplx(0.0d0,1.0d0)*vecJ(2,1)
        matJ[1][0] = vecJ[0][0]-vecJ[1][1]+I*vecJ[0][1]+I*vecJ[1][0];
        //matJSS(2,2)= vecJ(3,3)
        matJ[1][1] = vecJ[2][2];
        //matJSS(2,3)=dcmplx(-1.0d0,0.0d0)*vecJ(1,3)-dcmplx(0.0d0,1.0d0)*vecJ(2,3)
        matJ[1][2] =-vecJ[0][2]-I*vecJ[1][2];
        //matJSS(2,4)=dcmplx(-1.0d0,0.0d0)*vecJ(3,1)-dcmplx(0.0d0,1.0d0)*vecJ(3,2)
        matJ[1][3] =-vecJ[2][0]-I*vecJ[2][1];
        //matJSS(3,1)= vecJ(3,1)+dcmplx(0.0d0,1.0d0)*vecJ(3,2)
        matJ[2][0] = vecJ[2][0]+I*vecJ[2][1];
        //matJSS(3,2)=dcmplx(-1.0d0,0.0d0)*vecJ(1,3)+dcmplx(0.0d0,1.0d0)*vecJ(2,3)
        matJ[2][1] =-vecJ[0][2]+I*vecJ[1][2];
        //matJSS(3,3)=dcmplx(-1.0d0,0.0d0)*vecJ(3,3)
        matJ[2][2] =-vecJ[2][2];
        //matJSS(3,4)= vecJ(1,1)+vecJ(2,2)+dcmplx(0.0d0,1.0d0)*vecJ(1,2)-dcmplx(0.0d0,1.0d0)*vecJ(2,1)
        matJ[2][3] = vecJ[0][0]+vecJ[1][1]+I*vecJ[0][1]-I*vecJ[1][0];
        //matJSS(4,1)= vecJ(1,3)+dcmplx(0.0d0,1.0d0)*vecJ(2,3)
        matJ[3][0] = vecJ[0][2]+I*vecJ[1][2];
        //matJSS(4,2)=dcmplx(-1.0d0,0.0d0)*vecJ(3,1)+dcmplx(0.0d0,1.0d0)*vecJ(3,2)
        matJ[3][1] =-vecJ[2][0]+I*vecJ[2][1];
        //matJSS(4,3)= vecJ(1,1)+vecJ(2,2)-dcmplx(0.0d0,1.0d0)*vecJ(1,2)+dcmplx(0.0d0,1.0d0)*vecJ(2,1)
        matJ[3][2] = vecJ[0][0]+vecJ[1][1]-I*vecJ[0][1]+I*vecJ[1][0];
        //matJSS(4,4)=dcmplx(-1.0d0,0.0d0)*vecJ(3,3)
        matJ[3][3] =-vecJ[2][2];
        
        matJ2[1][1] = matJ[0][0]; 
        matJ2[1][2] = matJ[0][1]; 
        matJ2[1][3] = matJ[0][2]; 
        matJ2[1][0] = matJ[0][3]; 
        matJ2[2][1] = matJ[1][0]; 
        matJ2[2][2] = matJ[1][1]; 
        matJ2[2][3] = matJ[1][2]; 
        matJ2[2][0] = matJ[1][3]; 
        matJ2[3][1] = matJ[2][0]; 
        matJ2[3][2] = matJ[2][1];
        matJ2[3][3] = matJ[2][2]; 
        matJ2[3][0] = matJ[2][3]; 
        matJ2[0][1] = matJ[3][0]; 
        matJ2[0][2] = matJ[3][1]; 
        matJ2[0][3] = matJ[3][2]; 
        matJ2[0][0] = matJ[3][3]; 

        for(ellri=0; ellri<2; ellri++){
        for(ellrj=0; ellrj<2; ellrj++){
        for(ellrk=0; ellrk<2; ellrk++){
        for(ellrl=0; ellrl<2; ellrl++){
          for(elli1=0; ellri<2; ellri++){
          for(ellj1=0; ellrj<2; ellrj++){
          for(elli2=0; ellrk<2; ellrk++){
          for(ellj2=0; ellrl<2; ellrl++){
            
            iSSL1 = elli1*pow(2,mi) + ellj1*pow(2,mj) + ellri*pow(2,mri) + ellrj*pow(2,mrj) + ellrk*pow(2,mrk) + ellrl*pow(2,mrl);
            iSSL2 = elli2*pow(2,mi) + ellj2*pow(2,mj) + ellri*pow(2,mri) + ellrj*pow(2,mrj) + ellrk*pow(2,mrk) + ellrl*pow(2,mrl);
            iSS1  = elli1 + 2*ellj1;
            iSS2  = elli2 + 2*ellj2;
            matJL[iSSL1+64*iSSL2] += matJ2[iSS1][iSS2];
          }
          }
          }
          }
        }
        }
        }
        }
        
      }/* loop for ell */

      if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode subroutine iomp\n\n");}

      iomp=i_max/pow(2,ishift1+ishift2+ishift3+ishift4+ishift5+2); 
      if(myrank==0){printf("\n\n###Boost### SpinGC Boost mode subroutine iomp %d\n\n",iomp);}
      #pragma omp parallel default(none) private(arrayz,arrayx,ell4,ell5,ell6,m0,Ipart1,TRANSA,TRANSB,M,N,K,LDA,LDB,LDC,ALPHA,BETA,INFO) \
      firstprivate(matJL,iomp) shared(ishift1,ishift2,ishift3,ishift4,ishift5,pow4,pow5,pow41,pow51,tmp_v0,tmp_v1,tmp_v3)
      {
        c_malloc1(arrayz, (int)(64*pow(2,ishift4+ishift5-1)));  
        c_malloc1(arrayx, (int)(64*pow(2,ishift4+ishift5-1)));
        #pragma omp for
        //for(ell6 = 0; ell6 < i_max/pow(2,ishift1+ishift2+ishift3+ishift4+ishift5+2); ell6++){
        for(ell6 = 0; ell6 < iomp; ell6++){
          Ipart1=pow51*2*ell6;
          for(ell5 = 0; ell5 < (int)pow(2,ishift5); ell5++){
          for(ell4 = 0; ell4 < (int)pow(2,ishift4); ell4++){
            for(m0 = 0; m0 < 16; m0++){
              arrayz[0 + (int)(m0 +64*(ell4+ell5*pow(2,ishift4-1)))] = tmp_v1[1 + (int)(m0+16*ell4          +pow41*ell5+Ipart1)];
              arrayz[16+ (int)(m0 +64*(ell4+ell5*pow(2,ishift4-1)))] = tmp_v1[1 + (int)(m0+16*ell4+pow4     +pow41*ell5+Ipart1)];
              arrayz[32+ (int)(m0 +64*(ell4+ell5*pow(2,ishift4-1)))] = tmp_v1[1 + (int)(m0+16*ell4+pow5     +pow41*ell5+Ipart1)];
              arrayz[48+ (int)(m0 +64*(ell4+ell5*pow(2,ishift4-1)))] = tmp_v1[1 + (int)(m0+16*ell4+pow4+pow5+pow41*ell5+Ipart1)];
              tmp_v3[1 + (int)(m0+16*ell4          +pow41*ell5+Ipart1)]=tmp_v1[1 + (int)(m0+16*ell4          +pow41*ell5+Ipart1)];
              tmp_v3[1 + (int)(m0+16*ell4+pow4     +pow41*ell5+Ipart1)]=tmp_v1[1 + (int)(m0+16*ell4+pow4     +pow41*ell5+Ipart1)];
              tmp_v3[1 + (int)(m0+16*ell4+pow5     +pow41*ell5+Ipart1)]=tmp_v1[1 + (int)(m0+16*ell4+pow5     +pow41*ell5+Ipart1)];
              tmp_v3[1 + (int)(m0+16*ell4+pow4+pow5+pow41*ell5+Ipart1)]=tmp_v1[1 + (int)(m0+16*ell4+pow4+pow5+pow41*ell5+Ipart1)];
              arrayx[0 + (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)))] = tmp_v0[1 + (int)(m0+16*ell4          +pow41*ell5+Ipart1)];
              arrayx[16+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)))] = tmp_v0[1 + (int)(m0+16*ell4+pow4     +pow41*ell5+Ipart1)];
              arrayx[32+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)))] = tmp_v0[1 + (int)(m0+16*ell4+pow5     +pow41*ell5+Ipart1)];
              arrayx[48+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)))] = tmp_v0[1 + (int)(m0+16*ell4+pow4+pow5+pow41*ell5+Ipart1)];
            } 
          }
          }
          for(ell5 = 0; ell5 < (int)pow(2,ishift5); ell5++){
          for(ell4 = 0; ell4 < (int)pow(2,ishift4); ell4++){
            for(m0 = 0; m0 < 16; m0++){
              arrayz[0 + (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))] = tmp_v1[1 + m0+16*ell4          +pow41*ell5+pow51+Ipart1];
              arrayz[16+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))] = tmp_v1[1 + m0+16*ell4+pow4     +pow41*ell5+pow51+Ipart1];
              arrayz[32+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))] = tmp_v1[1 + m0+16*ell4+pow5     +pow41*ell5+pow51+Ipart1];
              arrayz[48+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))] = tmp_v1[1 + m0+16*ell4+pow4+pow5+pow41*ell5+pow51+Ipart1];
              tmp_v3[1 + m0+16*ell4          +pow41*ell5+pow51+Ipart1] = tmp_v1[1 + m0+16*ell4          +pow41*ell5+pow51+Ipart1];
              tmp_v3[1 + m0+16*ell4+pow4     +pow41*ell5+pow51+Ipart1] = tmp_v1[1 + m0+16*ell4+pow4     +pow41*ell5+pow51+Ipart1];
              tmp_v3[1 + m0+16*ell4+pow5     +pow41*ell5+pow51+Ipart1] = tmp_v1[1 + m0+16*ell4+pow5     +pow41*ell5+pow51+Ipart1];
              tmp_v3[1 + m0+16*ell4+pow4+pow5+pow41*ell5+pow51+Ipart1] = tmp_v1[1 + m0+16*ell4+pow4+pow5+pow41*ell5+pow51+Ipart1];
              arrayx[0 + (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))] = tmp_v0[1 + m0+16*ell4          +pow41*ell5+pow51+Ipart1];
              arrayx[16+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))] = tmp_v0[1 + m0+16*ell4+pow4     +pow41*ell5+pow51+Ipart1];
              arrayx[32+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))] = tmp_v0[1 + m0+16*ell4+pow5     +pow41*ell5+pow51+Ipart1];
              arrayx[48+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))] = tmp_v0[1 + m0+16*ell4+pow4+pow5+pow41*ell5+pow51+Ipart1];

            }
          }
          } 
          TRANSA = 'N';
          TRANSB = 'N';
          M = 64;
          N = (int)pow(2,ishift4+ishift5-1);
          K = 64;
          ALPHA = 1.0;
          LDA = 64;
          LDB = 64;
          BETA = 1.0;
          LDC = 64;
          zgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,matJL,&LDA,arrayz,&LDB,&BETA,arrayx,&LDC,&INFO);
          for(ell5 = 0; ell5 < pow(2,ishift5); ell5++){
          for(ell4 = 0; ell4 < pow(2,ishift4); ell4++){
            for(m0 = 0; m0 < 16; m0++){
              tmp_v1[1 + m0+16*ell4          +pow41*ell5+Ipart1]       = arrayx[0 + (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)))];
              tmp_v1[1 + m0+16*ell4+pow4     +pow41*ell5+Ipart1]       = arrayx[16+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)))];
              tmp_v1[1 + m0+16*ell4+pow5     +pow41*ell5+Ipart1]       = arrayx[32+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)))];
              tmp_v1[1 + m0+16*ell4+pow4+pow5+pow41*ell5+Ipart1]       = arrayx[48+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)))];
            }
          }
          }
          for(ell5 = 0; ell5 < (int)pow(2,ishift5); ell5++){
          for(ell4 = 0; ell4 < (int)pow(2,ishift4); ell4++){
            for(m0 = 0; m0 < 16; m0++){
              tmp_v1[1 + m0+16*ell4          +pow41*ell5+pow51+Ipart1] = arrayx[0 + (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))];
              tmp_v1[1 + m0+16*ell4+pow4     +pow41*ell5+pow51+Ipart1] = arrayx[16+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))];
              tmp_v1[1 + m0+16*ell4+pow5     +pow41*ell5+pow51+Ipart1] = arrayx[32+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))];
              tmp_v1[1 + m0+16*ell4+pow4+pow5+pow41*ell5+pow51+Ipart1] = arrayx[48+ (int)(m0+64*(ell4+ell5*pow(2,ishift4-1)+pow(2,ishift4+ishift5-2)))];
            }
          }
          }
         
        }/* omp parallel for */ 
        c_free1(arrayz, (int)(64*pow(2,ishift4+ishift5-1)));  
        c_free1(arrayx, (int)(64*pow(2,ishift4+ishift5-1)));  
      }/* omp parallel */
      if(pivot_flag==1){
        iomp=i_max/pow(2,ishift_nspin);
        #pragma omp parallel for default(none) private(ell4,ell5,ell6,m0,Ipart1,TRANSA,TRANSB,M,N,K,LDA,LDB,LDC,ALPHA,BETA) \
        firstprivate(i_max,iomp) shared(ishift1,ishift2,ishift3,ishift4,ishift5,pow4,pow5,pow41,pow51,ishift_nspin,tmp_v0,tmp_v1)
        //for(ell5 = 0; ell5 < i_max/pow(2,ishift_nspin); ell5++ ){
        for(ell5 = 0; ell5 < iomp; ell5++ ){
          for(ell4 = 0; ell4 < (int)pow(2,ishift_nspin); ell4++){
            tmp_v0[1 + (int)(ell5+(i_max/pow(2,ishift_nspin))*ell4)] = tmp_v1[1 + (int)(ell4+pow(2,ishift_nspin)*ell5)];
          } 
        }
        iomp=i_max/pow(2,ishift_nspin);
        #pragma omp parallel for default(none) private(ell4,ell5) \
        firstprivate(i_max,iomp) shared(ishift_nspin,tmp_v1,tmp_v3)
        //for(ell5 = 0; ell5 < i_max/pow(2,ishift_nspin); ell5++ ){
        for(ell5 = 0; ell5 < iomp; ell5++ ){
          for(ell4 = 0; ell4 < (int)pow(2,ishift_nspin); ell4++){
            tmp_v1[1 + (int)(ell5+(i_max/pow(2,ishift_nspin))*ell4)] = tmp_v3[1 + (int)(ell4+pow(2,ishift_nspin)*ell5)];
          } 
        }
      }
      else{ 
        #pragma omp parallel for default(none) private(ell4) \
        firstprivate(i_max) shared(tmp_v0,tmp_v1,tmp_v3)
        for(ell4 = 0; ell4 < i_max; ell4++ ){
          tmp_v0[1 + ell4] = tmp_v1[1 + ell4];
          tmp_v1[1 + ell4] =  tmp_v3[1 + ell4];
        }
      }/* if pivot_flag */

    }/* loop for j */
/* ucyy */

    ierr = MPI_Alltoall(&tmp_v1[1],(int)(i_max/nproc),MPI_DOUBLE_COMPLEX,&tmp_v3[1],(int)(i_max/nproc),MPI_COMM_WORLD, &statusMPI);
    ierr = MPI_Alltoall(&tmp_v0[1],(int)(i_max/nproc),MPI_DOUBLE_COMPLEX,&tmp_v2[1],(int)(i_max/nproc),MPI_COMM_WORLD, &statusMPI);

    iomp=pow(2,W0)/nproc;
    #pragma omp parallel for default(none) private(ell4,ell5,ell6) \
    firstprivate(i_max,W0,iomp,nproc) shared(tmp_v0,tmp_v1,tmp_v2,tmp_v3)
    //for(ell4 = 0; ell4 < pow(2,W0)/nproc; ell4++ ){
    for(ell4 = 0; ell4 < iomp; ell4++ ){
      for(ell5 = 0; ell5 < nproc; ell5++ ){
        for(ell6 = 0; ell6 < (int)(i_max/pow(2,W0)); ell6++ ){
          tmp_v1[1 + (int)(ell6+ell5*i_max/pow(2,W0)+ell4*i_max/(pow(2,W0)/nproc))] = tmp_v3[1 + (int)(ell6+ell4*i_max/pow(2,W0)+ell5*i_max/nproc)];
          tmp_v0[1 + (int)(ell6+ell5*i_max/pow(2,W0)+ell4*i_max/(pow(2,W0)/nproc))] = tmp_v2[1 + (int)(ell6+ell4*i_max/pow(2,W0)+ell5*i_max/nproc)];
        }
      }   
    }

  }/* loop for iloop */

/*
  dam_pr= X_child_general_int_spin_MPIBoost
    (
     matJ, X, tmp_v0, tmp_v1);
  
  X->Large.prdct += dam_pr;
*/

#endif
}/*void child_general_int_spin_MPIBoost*/

