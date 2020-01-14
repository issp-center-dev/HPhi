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

//Define Mode for mltply
// complex version

#ifdef MPI
#include "mpi.h"
#endif
#include "Common.h"
#include "common/setmemory.h"
#include "wrapperMPI.h"

void zgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double complex *ALPHA, double complex *matJL, int *LDA, double complex *arrayz, int *LDB, double complex *BETA, double complex *arrayx, int *LDC);

/**
 *
 * Exchange term in Spin model
 *
 * @author Mitsuaki Kawamura (The University of Tokyo)
 * @author Youhei Yamaji (The University of Tokyo)
 */
void general_int_spin_MPIBoost(
  struct BindStruct *X /**< [inout]*/,
  double complex *tmp_v0 /**< [out] Result v0 = H v1*/,
  double complex *tmp_v1 /**< [in] v0 = H v1*/,
  double complex *tmp_v2 /**< [inout] bufffer*/,
  double complex *tmp_v3 /**< [inout] bufffer*/
  )
{
#ifdef MPI
  
  //double complex dam_pr = 0;
  // MPI_Status statusMPI;

  //  int ierr;
  //  int INFO;
  char TRANSA, TRANSB;
  int M, N, K, LDA, LDB, LDC;
  double complex ALPHA, BETA;  
  long unsigned int i_max;
  long unsigned int j, k, ell, iloop;
  long unsigned int i1, i2;
  long unsigned int iomp;
  long unsigned int ell4, ell5, ell6, m0, Ipart1;
  long unsigned int mi, mj, mri, mrj, mrk, mrl;
  int indj;
  long unsigned int ellrl, ellrk, ellrj, ellri, elli1, elli2, ellj1, ellj2;
  long unsigned int iSS1, iSS2, iSSL1, iSSL2;
  double complex **vecJ;
  double complex **matJ, **matJ2;
  double complex *matJL;
  double complex *matI;
  double complex **matB;
  double complex *arrayz;
  double complex *arrayx;
  double complex *arrayw;
  long unsigned int ishift1, ishift2, ishift3, ishift4, ishift5, pivot_flag, num_J_star;
  long unsigned int pow4, pow5, pow41, pow51;  
  //long unsigned int pow1, pow2, pow3, pow4, pow5, pow11, pow21, pow31, pow41, pow51; 

  i_max = X->Check.idim_max;

/*
//zero clear
  #pragma omp parallel for default(none) private(j) \
  shared(i_max,tmp_v0)
  for(j=0;j<i_max;j++){
    tmp_v0[j+1]=0.0;
  }
*/

  vecJ = cd_2d_allocate( 3, 3);
  matJ = cd_2d_allocate(4, 4);
  matJ2 = cd_2d_allocate(4, 4);
  matB = cd_2d_allocate(2,2);
  matJL = cd_1d_allocate(64*64);
  matI = cd_1d_allocate(64*64);

  //defmodelBoost(X->Boost.W0, X->Boost.R0, X->Boost.num_pivot, X->Boost.ishift_nspin, X->Boost.list_6spin_star, X->Boost.list_6spin_pair, 1, X->Boost.arrayJ, X->Boost.vecB);
  
  for(iloop=0; iloop < X->Boost.R0; iloop++){


    for(j=iloop*X->Boost.num_pivot; j < (iloop+1)*X->Boost.num_pivot; j++){
      
      num_J_star = (long unsigned int)X->Boost.list_6spin_star[j][0]; //(0,j) 
      ishift1    = (long unsigned int)X->Boost.list_6spin_star[j][1]; //(1,j) 
      ishift2    = (long unsigned int)X->Boost.list_6spin_star[j][2]; //(2,j) 
      ishift3    = (long unsigned int)X->Boost.list_6spin_star[j][3]; //(3,j)
      ishift4    = (long unsigned int)X->Boost.list_6spin_star[j][4]; //(4,j)
      ishift5    = (long unsigned int)X->Boost.list_6spin_star[j][5]; //(5,j)
      pivot_flag = (long unsigned int)X->Boost.list_6spin_star[j][6]; //(6,j)
      //pow1 = (int)pow(2.0,ishift1);
      //pow2 = (int)pow(2.0,ishift1+ishift2);
      //pow3 = (int)pow(2.0,ishift1+ishift2+ishift3);
      pow4 = (int)pow(2.0,ishift1+ishift2+ishift3+ishift4);
      pow5 = (int)pow(2.0,ishift1+ishift2+ishift3+ishift4+ishift5);
      //pow11= (int)pow(2.0,ishift1+1);
      //pow21= (int)pow(2.0,ishift1+ishift2+1);
      //pow31= (int)pow(2.0,ishift1+ishift2+ishift3+1);
      pow41= (int)pow(2.0,ishift1+ishift2+ishift3+ishift4+1);
      pow51= (int)pow(2.0,ishift1+ishift2+ishift3+ishift4+ishift5+1);

      for(k=0; k < (64*64); k++){
        matJL[k] = 0.0 + 0.0*I;
        matI[k]  = 0.0 + 0.0*I;
      }
      for(k=0; k < 64; k++){
        matI[k+64*k] = 1.0;
      }

      for(ell=0; ell < num_J_star; ell++){
        mi   = (long unsigned int)X->Boost.list_6spin_pair[j][0][ell]; //(1,ell,j)
        mj   = (long unsigned int)X->Boost.list_6spin_pair[j][1][ell]; //(2,ell,j)
        mri  = (long unsigned int)X->Boost.list_6spin_pair[j][2][ell]; //(3,ell,j)
        mrj  = (long unsigned int)X->Boost.list_6spin_pair[j][3][ell]; //(4,ell,j)
        mrk  = (long unsigned int)X->Boost.list_6spin_pair[j][4][ell]; //(5,ell,j)
        mrl  = (long unsigned int)X->Boost.list_6spin_pair[j][5][ell]; //(6,ell,j)
        indj = X->Boost.list_6spin_pair[j][6][ell]; //(7,ell,j)
        for(i1 = 0; i1 < 3; i1++){
          for(i2 = 0; i2 < 3; i2++){
            vecJ[i1][i2] = X->Boost.arrayJ[(indj-1)][i1][i2];
          }
        } 
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
        matJ[1][2] =(-1.0)*vecJ[0][2]-I*vecJ[1][2];
        //matJSS(2,4)=dcmplx(-1.0d0,0.0d0)*vecJ(3,1)-dcmplx(0.0d0,1.0d0)*vecJ(3,2)
        matJ[1][3] =(-1.0)*vecJ[2][0]-I*vecJ[2][1];
        //matJSS(3,1)= vecJ(3,1)+dcmplx(0.0d0,1.0d0)*vecJ(3,2)
        matJ[2][0] = vecJ[2][0]+I*vecJ[2][1];
        //matJSS(3,2)=dcmplx(-1.0d0,0.0d0)*vecJ(1,3)+dcmplx(0.0d0,1.0d0)*vecJ(2,3)
        matJ[2][1] =(-1.0)*vecJ[0][2]+I*vecJ[1][2];
        //matJSS(3,3)=dcmplx(-1.0d0,0.0d0)*vecJ(3,3)
        matJ[2][2] =(-1.0)*vecJ[2][2];
        //matJSS(3,4)= vecJ(1,1)+vecJ(2,2)+dcmplx(0.0d0,1.0d0)*vecJ(1,2)-dcmplx(0.0d0,1.0d0)*vecJ(2,1)
        matJ[2][3] = vecJ[0][0]+vecJ[1][1]+I*vecJ[0][1]-I*vecJ[1][0];
        //matJSS(4,1)= vecJ(1,3)+dcmplx(0.0d0,1.0d0)*vecJ(2,3)
        matJ[3][0] = vecJ[0][2]+I*vecJ[1][2];
        //matJSS(4,2)=dcmplx(-1.0d0,0.0d0)*vecJ(3,1)+dcmplx(0.0d0,1.0d0)*vecJ(3,2)
        matJ[3][1] =(-1.0)*vecJ[2][0]+I*vecJ[2][1];
        //matJSS(4,3)= vecJ(1,1)+vecJ(2,2)-dcmplx(0.0d0,1.0d0)*vecJ(1,2)+dcmplx(0.0d0,1.0d0)*vecJ(2,1)
        matJ[3][2] = vecJ[0][0]+vecJ[1][1]-I*vecJ[0][1]+I*vecJ[1][0];
        //matJSS(4,4)=dcmplx(-1.0d0,0.0d0)*vecJ(3,3)
        matJ[3][3] =(-1.0)*vecJ[2][2];
        
        matJ2[3][3] = matJ[0][0]; 
        matJ2[3][0] = matJ[0][1]; 
        matJ2[3][1] = matJ[0][2]; 
        matJ2[3][2] = matJ[0][3]; 
        matJ2[0][3] = matJ[1][0]; 
        matJ2[0][0] = matJ[1][1]; 
        matJ2[0][1] = matJ[1][2]; 
        matJ2[0][2] = matJ[1][3]; 
        matJ2[1][3] = matJ[2][0]; 
        matJ2[1][0] = matJ[2][1];
        matJ2[1][1] = matJ[2][2]; 
        matJ2[1][2] = matJ[2][3]; 
        matJ2[2][3] = matJ[3][0]; 
        matJ2[2][0] = matJ[3][1]; 
        matJ2[2][1] = matJ[3][2]; 
        matJ2[2][2] = matJ[3][3]; 

        for(ellri=0; ellri<2; ellri++){
        for(ellrj=0; ellrj<2; ellrj++){
        for(ellrk=0; ellrk<2; ellrk++){
        for(ellrl=0; ellrl<2; ellrl++){
          for(elli1=0; elli1<2; elli1++){
          for(ellj1=0; ellj1<2; ellj1++){
          for(elli2=0; elli2<2; elli2++){
          for(ellj2=0; ellj2<2; ellj2++){
            
            iSSL1 = elli1*(int)pow(2,mi) + ellj1*(int)pow(2,mj) + ellri*(int)pow(2,mri) + ellrj*(int)pow(2,mrj) + ellrk*(int)pow(2,mrk) + ellrl*(int)pow(2,mrl);
            iSSL2 = elli2*(int)pow(2,mi) + ellj2*(int)pow(2,mj) + ellri*(int)pow(2,mri) + ellrj*(int)pow(2,mrj) + ellrk*(int)pow(2,mrk) + ellrl*(int)pow(2,mrl);
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

      /* external magnetic field B */
      if(pivot_flag==1){
        matB[0][0] = + X->Boost.vecB[2]; // -BM
        matB[1][1] = - X->Boost.vecB[2]; // -BM
        //matB[0][1] = - X->Boost.vecB[0] + I*X->Boost.vecB[1]; // -BM
        //matB[1][0] = - X->Boost.vecB[0] - I*X->Boost.vecB[1]; // -BM
        matB[0][1] = - X->Boost.vecB[0] - I*X->Boost.vecB[1]; // -BM
        matB[1][0] = - X->Boost.vecB[0] + I*X->Boost.vecB[1]; // -BM
        for(ellri=0; ellri<2; ellri++){
        for(ellrj=0; ellrj<2; ellrj++){
        for(ellrk=0; ellrk<2; ellrk++){
        for(ellrl=0; ellrl<2; ellrl++){
        for(ellj1=0; ellj1<2; ellj1++){
          for(elli1=0; elli1<2; elli1++){
          for(elli2=0; elli2<2; elli2++){
            for(ellj2=0; ellj2<X->Boost.ishift_nspin; ellj2++){
              iSSL1 = elli1*(int)pow(2,ellj2) + ellj1*(int)pow(2,((ellj2+1)%6)) + ellri*(int)pow(2,((ellj2+2)%6)) + ellrj*(int)pow(2,((ellj2+3)%6)) + ellrk*(int)pow(2,((ellj2+4)%6)) + ellrl*(int)pow(2,((ellj2+5)%6));
              iSSL2 = elli2*(int)pow(2,ellj2) + ellj1*(int)pow(2,((ellj2+1)%6)) + ellri*(int)pow(2,((ellj2+2)%6)) + ellrj*(int)pow(2,((ellj2+3)%6)) + ellrk*(int)pow(2,((ellj2+4)%6)) + ellrl*(int)pow(2,((ellj2+5)%6));
              matJL[iSSL1+64*iSSL2] += matB[elli1][elli2];
            }
          } 
          } 
        }
        }
        }
        }
        }
      }
      /* external magnetic field B */
    
      iomp=i_max/(int)pow(2.0,ishift1+ishift2+ishift3+ishift4+ishift5+2);

      #pragma omp parallel default(none) private(arrayx,arrayz,arrayw,ell4,ell5,ell6,m0,Ipart1,TRANSA,TRANSB,M,N,K,LDA,LDB,LDC,ALPHA,BETA) \
      shared(matJL,matI,iomp,i_max,myrank,ishift1,ishift2,ishift3,ishift4,ishift5,pow4,pow5,pow41,pow51,tmp_v0,tmp_v1,tmp_v3)
      {

        arrayx = cd_1d_allocate(64*((int)pow(2.0,ishift4+ishift5-1)));
        arrayz = cd_1d_allocate(64*((int)pow(2.0,ishift4+ishift5-1)));
        arrayw = cd_1d_allocate(64*((int)pow(2.0,ishift4+ishift5-1)));

#pragma omp for
        for(ell6 = 0; ell6 < iomp; ell6++){
          Ipart1=pow51*2*ell6;
          for(ell5 = 0; ell5 < (int)pow(2.0, ishift5-1); ell5++){
            for(ell4 = 0; ell4 < (int)pow(2.0, ishift4-1); ell4++){
              for(m0 = 0; m0 < 16; m0++){        
                arrayz[(0 + m0 +64*(ell4+ell5*(int)pow(2.0,ishift4-1)))] = tmp_v1[(1 + m0+16*ell4          +pow41*ell5+Ipart1)];
                arrayz[(16+ m0 +64*(ell4+ell5*(int)pow(2.0,ishift4-1)))] = tmp_v1[(1 + m0+16*ell4+pow4     +pow41*ell5+Ipart1)];
                arrayz[(32+ m0 +64*(ell4+ell5*(int)pow(2.0,ishift4-1)))] = tmp_v1[(1 + m0+16*ell4+pow5     +pow41*ell5+Ipart1)];
                arrayz[(48+ m0 +64*(ell4+ell5*(int)pow(2.0,ishift4-1)))] = tmp_v1[(1 + m0+16*ell4+pow4+pow5+pow41*ell5+Ipart1)];
                tmp_v3[(1 + m0+16*ell4          +pow41*ell5+Ipart1)]=tmp_v1[(1 + m0+16*ell4          +pow41*ell5+Ipart1)];
                tmp_v3[(1 + m0+16*ell4+pow4     +pow41*ell5+Ipart1)]=tmp_v1[(1 + m0+16*ell4+pow4     +pow41*ell5+Ipart1)];
                tmp_v3[(1 + m0+16*ell4+pow5     +pow41*ell5+Ipart1)]=tmp_v1[(1 + m0+16*ell4+pow5     +pow41*ell5+Ipart1)];
                tmp_v3[(1 + m0+16*ell4+pow4+pow5+pow41*ell5+Ipart1)]=tmp_v1[(1 + m0+16*ell4+pow4+pow5+pow41*ell5+Ipart1)];
                arrayx[(0 + m0 +64*(ell4+ell5*(int)pow(2.0,ishift4-1)))] = tmp_v0[(1 + m0+16*ell4          +pow41*ell5+Ipart1)];
                arrayx[(16+ m0 +64*(ell4+ell5*(int)pow(2.0,ishift4-1)))] = tmp_v0[(1 + m0+16*ell4+pow4     +pow41*ell5+Ipart1)];
                arrayx[(32+ m0 +64*(ell4+ell5*(int)pow(2.0,ishift4-1)))] = tmp_v0[(1 + m0+16*ell4+pow5     +pow41*ell5+Ipart1)];
                arrayx[(48+ m0 +64*(ell4+ell5*(int)pow(2.0,ishift4-1)))] = tmp_v0[(1 + m0+16*ell4+pow4+pow5+pow41*ell5+Ipart1)];
              } 
            }
          }
          
          
          for(ell5 = 0; ell5 < (int)pow(2.0, ishift5-1); ell5++){
            for(ell4 = 0; ell4 < (int)pow(2.0, ishift4-1); ell4++){
              for(m0 = 0; m0 < 16; m0++){
                arrayz[(0 + m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))] = tmp_v1[(1 + m0+16*ell4          +pow41*ell5+pow51+Ipart1)];
                arrayz[(16+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))] = tmp_v1[(1 + m0+16*ell4+pow4     +pow41*ell5+pow51+Ipart1)];
                arrayz[(32+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))] = tmp_v1[(1 + m0+16*ell4+pow5     +pow41*ell5+pow51+Ipart1)];
                arrayz[(48+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))] = tmp_v1[(1 + m0+16*ell4+pow4+pow5+pow41*ell5+pow51+Ipart1)];
                tmp_v3[(1 + m0+16*ell4          +pow41*ell5+pow51+Ipart1)] = tmp_v1[(1 + m0+16*ell4          +pow41*ell5+pow51+Ipart1)];
                tmp_v3[(1 + m0+16*ell4+pow4     +pow41*ell5+pow51+Ipart1)] = tmp_v1[(1 + m0+16*ell4+pow4     +pow41*ell5+pow51+Ipart1)];
                tmp_v3[(1 + m0+16*ell4+pow5     +pow41*ell5+pow51+Ipart1)] = tmp_v1[(1 + m0+16*ell4+pow5     +pow41*ell5+pow51+Ipart1)];
                tmp_v3[(1 + m0+16*ell4+pow4+pow5+pow41*ell5+pow51+Ipart1)] = tmp_v1[(1 + m0+16*ell4+pow4+pow5+pow41*ell5+pow51+Ipart1)];
                arrayx[(0 + m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))] = tmp_v0[(1 + m0+16*ell4          +pow41*ell5+pow51+Ipart1)];
                arrayx[(16+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))] = tmp_v0[(1 + m0+16*ell4+pow4     +pow41*ell5+pow51+Ipart1)];
                arrayx[(32+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))] = tmp_v0[(1 + m0+16*ell4+pow5     +pow41*ell5+pow51+Ipart1)];
                arrayx[(48+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))] = tmp_v0[(1 + m0+16*ell4+pow4+pow5+pow41*ell5+pow51+Ipart1)];
              }
            
            }
          } 
          
          TRANSA = 'N';
          TRANSB = 'N';
          M = 64;
          N = (int)pow(2.0, ishift4+ishift5-1);
          K = 64;
          ALPHA = 1.0;
          LDA = 64;
          LDB = 64;
          BETA = 1.0;
          LDC = 64;
          
                zgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,matJL,&LDA,arrayz,&LDB,&BETA,arrayx,&LDC);
                //zgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,matI,&LDA,arrayz,&LDB,&BETA,arrayx,&LDC);
/*          
          for(ell5=0;ell5<(64*N);ell5++){
            arrayw[ell5]=0.0;
          }
          for(ell5=0;ell5<64;ell5++){
            for(ell4=0;ell4<64;ell4++){
              for(m0=0;m0<N;m0++){
                arrayw[(ell5+64*m0)] += matJL[(ell5+64*ell4)]*arrayz[(ell4+64*m0)];
              }
            }
          }
          for(ell5=0;ell5<64*N;ell5++){
            arrayx[ell5] += arrayw[ell5];
          }
*/          
        


          for(ell5 = 0; ell5 < (int)pow(2.0,ishift5-1); ell5++){
          for(ell4 = 0; ell4 < (int)pow(2.0,ishift4-1); ell4++){
            for(m0 = 0; m0 < 16; m0++){
              tmp_v1[(1 + m0+16*ell4          +pow41*ell5+Ipart1)]       = arrayx[(0 + m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)))];
              tmp_v1[(1 + m0+16*ell4+pow4     +pow41*ell5+Ipart1)]       = arrayx[(16+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)))];
              tmp_v1[(1 + m0+16*ell4+pow5     +pow41*ell5+Ipart1)]       = arrayx[(32+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)))];
              tmp_v1[(1 + m0+16*ell4+pow4+pow5+pow41*ell5+Ipart1)]       = arrayx[(48+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)))];
            }
          }
          }
          for(ell5 = 0; ell5 < (int)pow(2.0,ishift5-1); ell5++){
          for(ell4 = 0; ell4 < (int)pow(2.0,ishift4-1); ell4++){
            for(m0 = 0; m0 < 16; m0++){
              tmp_v1[(1 + m0+16*ell4          +pow41*ell5+pow51+Ipart1)] = arrayx[(0 + m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))];
              tmp_v1[(1 + m0+16*ell4+pow4     +pow41*ell5+pow51+Ipart1)] = arrayx[(16+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))];
              tmp_v1[(1 + m0+16*ell4+pow5     +pow41*ell5+pow51+Ipart1)] = arrayx[(32+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))];
              tmp_v1[(1 + m0+16*ell4+pow4+pow5+pow41*ell5+pow51+Ipart1)] = arrayx[(48+ m0+64*(ell4+ell5*(int)pow(2.0,ishift4-1)+(int)pow(2.0,ishift4+ishift5-2)))];
            }
          }
          }

        }/* omp parallel for */
        free_cd_1d_allocate(arrayz);
        free_cd_1d_allocate(arrayx);
        free_cd_1d_allocate(arrayw);
      }/* omp parallel */

      if(pivot_flag==1){
        iomp=i_max/(int)pow(2.0,X->Boost.ishift_nspin);
        #pragma omp parallel for default(none) private(ell4,ell5,ell6,m0,Ipart1,TRANSA,TRANSB,M,N,K,LDA,LDB,LDC,ALPHA,BETA) \
        firstprivate(iomp) shared(i_max,ishift1,ishift2,ishift3,ishift4,ishift5,pow4,pow5,pow41,pow51,X,tmp_v0,tmp_v1)
        for(ell5 = 0; ell5 < iomp; ell5++ ){
          for(ell4 = 0; ell4 < (int)pow(2.0,X->Boost.ishift_nspin); ell4++){
            tmp_v0[(1 + ell5+(i_max/(int)pow(2.0,X->Boost.ishift_nspin))*ell4)] = tmp_v1[(1 + ell4+((int)pow(2.0,X->Boost.ishift_nspin))*ell5)];
          } 
        }
        iomp=i_max/(int)pow(2.0,X->Boost.ishift_nspin);
        #pragma omp parallel for default(none) private(ell4,ell5) \
        firstprivate(iomp) shared(i_max,X,tmp_v1,tmp_v3)
        for(ell5 = 0; ell5 < iomp; ell5++ ){
          for(ell4 = 0; ell4 < (int)pow(2.0,X->Boost.ishift_nspin); ell4++){
            tmp_v1[(1 + ell5+(i_max/(int)pow(2.0,X->Boost.ishift_nspin))*ell4)] = tmp_v3[(1 + ell4+((int)pow(2.0,X->Boost.ishift_nspin))*ell5)];
          } 
        }
      }
      else{ 
        #pragma omp parallel for default(none) private(ell4) \
        shared(i_max,tmp_v0,tmp_v1,tmp_v3)
        for(ell4 = 0; ell4 < i_max; ell4++ ){
          tmp_v0[1 + ell4] = tmp_v1[1 + ell4];
          tmp_v1[1 + ell4] = tmp_v3[1 + ell4];
        }
      }/* if pivot_flag */

    }/* loop for j */

    /*
    ierr = MPI_Alltoall(&tmp_v1[1],(int)(i_max/nproc),MPI_DOUBLE_COMPLEX,&tmp_v3[1],(int)(i_max/nproc),MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD);
    ierr = MPI_Alltoall(&tmp_v0[1],(int)(i_max/nproc),MPI_DOUBLE_COMPLEX,&tmp_v2[1],(int)(i_max/nproc),MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD);
     */
    MPI_Alltoall(&tmp_v1[1],(int)(i_max/nproc),MPI_DOUBLE_COMPLEX,&tmp_v3[1],(int)(i_max/nproc),MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD);
    MPI_Alltoall(&tmp_v0[1],(int)(i_max/nproc),MPI_DOUBLE_COMPLEX,&tmp_v2[1],(int)(i_max/nproc),MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD);


    iomp=(int)pow(2.0,X->Boost.W0)/nproc;
    #pragma omp parallel for default(none) private(ell4,ell5,ell6) \
    firstprivate(iomp) shared(i_max,X,nproc,tmp_v0,tmp_v1,tmp_v2,tmp_v3)
    //for(ell4 = 0; ell4 < (int)pow(2.0,X->Boost.W0)/nproc; ell4++ ){
    for(ell4 = 0; ell4 < iomp; ell4++ ){
      for(ell5 = 0; ell5 < nproc; ell5++ ){
        for(ell6 = 0; ell6 < (int)(i_max/(int)pow(2.0,X->Boost.W0)); ell6++ ){
          tmp_v1[(1 + ell6+ell5*i_max/(int)pow(2.0,X->Boost.W0)+ell4*i_max/((int)pow(2.0,X->Boost.W0)/nproc))] = tmp_v3[(1 + ell6+ell4*i_max/(int)pow(2.0,X->Boost.W0)+ell5*i_max/nproc)];
          tmp_v0[(1 + ell6+ell5*i_max/(int)pow(2.0,X->Boost.W0)+ell4*i_max/((int)pow(2.0,X->Boost.W0)/nproc))] = tmp_v2[(1 + ell6+ell4*i_max/(int)pow(2.0,X->Boost.W0)+ell5*i_max/nproc)];
        }
      }   
    }


  }/* loop for iloop */

/*
  dam_pr= child_general_int_spin_MPIBoost
    (
     matJ, X, tmp_v0, tmp_v1);
  
  X->Large.prdct += dam_pr;
*/
//  c_free1(arrayz, (int)pow(2.0, 16));
//  c_free1(arrayx, (int)pow(2.0, 16));
//  c_free1(arrayw, (int)pow(2.0, 16));

  free_cd_2d_allocate(vecJ);
  free_cd_2d_allocate(matJ);
  free_cd_2d_allocate(matJ2);
  free_cd_2d_allocate(matB);
  free_cd_1d_allocate(matJL);
  free_cd_1d_allocate(matI);
#endif
  
}/*void general_int_spin_MPIBoost*/

