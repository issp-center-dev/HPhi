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

#include <bitcalc.h>
#include "common/setmemory.h"
#include "FileIO.h"
#include "sz.h"
#include "wrapperMPI.h"
#include "xsetmem.h"

/**
 * @file   sz.c
 * 
 * @brief  Generating Hilbert spaces
 * 
 * @version 0.2
 * @details 
 *
 * @version 0.1
 *
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 * 
 */

int omp_sz_Kondo_hacker(
        long unsigned int ib,
        long unsigned int ihfbit,
        struct BindStruct *X,
        long unsigned int *list_1_,
        long unsigned int *list_2_1_,
        long unsigned int *list_2_2_,
        long unsigned int *list_jb_
);

/** 
 * 
 * @brief generating Hilbert space
 * 
 * @param[inout] X 
 * @param[out] list_1_   list_1[icnt] = i (index of full Hilbert space) : icnt = index in the restricted Hilbert space 
 * @param[out] list_2_1_ icnt=list_2_1[]+list_2_2[] 
 * @param[out] list_2_2_ 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int sz
(
 struct BindStruct *X,
 long unsigned int *list_1_,
 long unsigned int *list_2_1_,
 long unsigned int *list_2_2_
 )
{
  FILE *fp,*fp_err;
  char sdt[D_FileNameMax],sdt_err[D_FileNameMax];
  long unsigned int *HilbertNumToSz;
  long unsigned int i,icnt; 
  long unsigned int ib,jb,ib_start,ib_end, sdim_div, sdim_rest;
    
  long unsigned int j;
  long unsigned int div;
  long unsigned int num_up,num_down;
  long unsigned int irght,ilft,ihfbit;
  long unsigned int *jbthread;

  //*[s] for omp parall
  int mythread;
  unsigned int  all_up,all_down,tmp_res,num_threads;
  long unsigned int tmp_1,tmp_2,tmp_3;
  long int **comb, **comb2;
  //*[e] for omp parall

  // [s] for Kondo
  unsigned int N_all_up, N_all_down;
  unsigned int all_loc;
  long unsigned int num_loc, div_down;
  unsigned int num_loc_up;
  int icheck_loc;
  int ihfSpinDown=0;
  // [e] for Kondo
    
  long unsigned int i_max=0;
  double idim=0.0;
  long unsigned int div_up;

  // [s] for general spin
  long int *list_2_1_Sz = NULL;
  long int *list_2_2_Sz = NULL;
  if(X->Def.iFlgGeneralSpin==TRUE){
    list_2_1_Sz = li_1d_allocate(X->Check.sdim+2);
    list_2_2_Sz = li_1d_allocate((X->Def.Tpow[X->Def.Nsite-1]*X->Def.SiteToBit[X->Def.Nsite-1]/X->Check.sdim)+2);
    for(j=0; j<X->Check.sdim+2;j++){
      list_2_1_Sz[j]=0;
      }
    for(j=0; j< (X->Def.Tpow[X->Def.Nsite-1]*X->Def.SiteToBit[X->Def.Nsite-1]/X->Check.sdim)+2; j++){
      list_2_2_Sz[j]=0;
    }
  }
  // [e] for general spin

  long unsigned int *list_jb;
  list_jb = lui_1d_allocate(X->Large.SizeOflistjb);
  for(i=0; i<X->Large.SizeOflistjb; i++){
    list_jb[i]=0;
  }

//hacker
  int hacker;
  long unsigned int tmp_i,tmp_j,tmp_pow,max_tmp_i;
  long unsigned int ia,ja;
  long unsigned int ibpatn=0;
//hacker

  int iSpnup, iMinup,iAllup;
  unsigned int N2=0;
  unsigned int N=0;
  fprintf(stdoutMPI, "%s", cProStartCalcSz);
  TimeKeeper(X, cFileNameSzTimeKeep, cInitalSz, "w");
  TimeKeeper(X, cFileNameTimeKeep, cInitalSz, "a");

  if(X->Check.idim_max!=0){
  switch(X->Def.iCalcModel){
  case HubbardGC:
  case HubbardNConserved:
  case Hubbard:
    N2=2*X->Def.Nsite;
    idim = pow(2.0,N2);
    break;
  case KondoGC:
  case Kondo:
    N2  = 2*X->Def.Nsite;
    N  =  X->Def.Nsite;
    idim = pow(2.0,N2);
    for(j=0;j<N;j++){
      fprintf(stdoutMPI, cStateLocSpin,j,X->Def.LocSpn[j]);
    }
    break;
  case SpinGC:
  case Spin:
    N=X->Def.Nsite;
    if(X->Def.iFlgGeneralSpin==FALSE){
      idim = pow(2.0, N);
    }
    else{
      idim=1;
      for(j=0; j<N; j++){
        idim *= X->Def.SiteToBit[j];
      }
    }
    break;
  }
  comb = li_2d_allocate(X->Def.Nsite+1,X->Def.Nsite+1);
  i_max=X->Check.idim_max;
  
  switch(X->Def.iCalcModel){
  case HubbardNConserved:
  case Hubbard:
  case KondoGC:
  case Kondo:
  case Spin:
    if(X->Def.iFlgGeneralSpin==FALSE){
      if(GetSplitBitByModel(X->Def.Nsite, X->Def.iCalcModel, &irght, &ilft, &ihfbit)!=0){
        exitMPI(-1);
      }
      X->Large.irght=irght;
      X->Large.ilft=ilft;
      X->Large.ihfbit=ihfbit;
      //fprintf(stdoutMPI, "idim=%lf irght=%ld ilft=%ld ihfbit=%ld \n",idim,irght,ilft,ihfbit);
    }
     else{
      ihfbit=X->Check.sdim;
      //fprintf(stdoutMPI, "idim=%lf ihfbit=%ld \n",idim, ihfbit);
    }
    break;
  default:
    break;
  }
  
  icnt=1;
  jb=0;

  if(X->Def.READ==1){
    if(Read_sz(X, irght, ilft, ihfbit, &i_max)!=0){
      exitMPI(-1);
    }
  }
  else{ 
    sprintf(sdt, cFileNameSzTimeKeep, X->Def.CDataFileHead);
#ifdef _OPENMP
    num_threads  = omp_get_max_threads();
#else
    num_threads=1;
#endif
    childfopenMPI(sdt,"a", &fp);
    fprintf(fp, "num_threads==%d\n",num_threads);
    fclose(fp);
    
    //*[s] omp parallel

    TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzStart, "a");
    TimeKeeper(X, cFileNameTimeKeep, cOMPSzStart, "a");
    switch(X->Def.iCalcModel){
    case HubbardGC:
      icnt = X->Def.Tpow[2*X->Def.Nsite-1]*2+0;/*Tpow[2*X->Def.Nsit]=1*/
      break;
      
    case SpinGC:
      if(X->Def.iFlgGeneralSpin==FALSE){
        icnt = X->Def.Tpow[X->Def.Nsite-1]*2+0;/*Tpow[X->Def.Nsit]=1*/
      }
      else{
        icnt = X->Def.Tpow[X->Def.Nsite-1]*X->Def.SiteToBit[X->Def.Nsite-1];
      }
      break;
      
    case KondoGC:
      // this part can not be parallelized
      jb = 0;
      num_loc=0;
      for(j=X->Def.Nsite/2; j< X->Def.Nsite ;j++){ // counting # of localized spins
        if(X->Def.LocSpn[j] != ITINERANT){ // //ITINERANT ==0 -> itinerant
          num_loc += 1;
        }
      }

      for(ib=0;ib<X->Check.sdim;ib++){
        list_jb[ib]=jb;
        i=ib*ihfbit;
        icheck_loc=1;
        for(j=(X->Def.Nsite+1)/2; j< X->Def.Nsite ;j++){
          div_up    = i & X->Def.Tpow[2*j];
          div_up    = div_up/X->Def.Tpow[2*j];
          div_down  = i & X->Def.Tpow[2*j+1];
          div_down  = div_down/X->Def.Tpow[2*j+1];
          if(X->Def.LocSpn[j] != ITINERANT){
            if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
              icheck_loc= icheck_loc;
            }
            else{
              icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
            }
          }
        }
        if(icheck_loc == 1){
          if(X->Def.Nsite%2==1 && X->Def.LocSpn[X->Def.Nsite/2] != ITINERANT){
            jb +=X->Def.Tpow[X->Def.Nsite-1-(X->Def.NLocSpn-num_loc)];
          }else{
            jb +=X->Def.Tpow[X->Def.Nsite-(X->Def.NLocSpn-num_loc)];
          }
        }
      }

      icnt = 0; 
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N2, X) shared(list_1_, list_2_1_, list_2_2_, list_jb)
      for(ib=0;ib<X->Check.sdim;ib++){
        icnt+=omp_sz_KondoGC(ib, ihfbit, X, list_1_, list_2_1_, list_2_2_, list_jb);
      }      
    break;

    case Hubbard:
      hacker = X->Def.read_hacker;
      if(hacker==0){
        // this part can not be parallelized
        jb = 0;
        for(ib=0;ib<X->Check.sdim;ib++){ // sdim = 2^(N/2)
          list_jb[ib] = jb;
          i           = ib*ihfbit;
          //[s] counting # of up and down electrons
          num_up      = 0;
          for(j=0;j<=N2-2;j+=2){ // even -> up spin
            div=i & X->Def.Tpow[j];
            div=div/X->Def.Tpow[j];
            num_up+=div;
          }
          num_down=0;
          for(j=1;j<=N2-1;j+=2){ // odd -> down spin
            div=i & X->Def.Tpow[j];
            div=div/X->Def.Tpow[j];
            num_down+=div;
          }
          //[e] counting # of up and down electrons
          tmp_res  = X->Def.Nsite%2; // even Ns-> 0, odd Ns -> 1
          all_up   = (X->Def.Nsite+tmp_res)/2;
          all_down = (X->Def.Nsite-tmp_res)/2;

          tmp_1 = Binomial(all_up,X->Def.Nup-num_up,comb,all_up);
          tmp_2 = Binomial(all_down,X->Def.Ndown-num_down,comb,all_down);
          jb   += tmp_1*tmp_2;
        }

        //#pragma omp barrier
        TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzMid, "a");
        TimeKeeper(X, cFileNameTimeKeep, cOMPSzMid, "a");

        icnt = 0;
        for(ib=0;ib<X->Check.sdim;ib++){
          icnt+=omp_sz(ib,ihfbit, X, list_1_, list_2_1_, list_2_2_, list_jb);
        }
        break;
      }else if(hacker==1){
        jbthread = lui_1d_allocate(nthreads);
        #pragma omp parallel default(none) \
        shared(X,list_jb,ihfbit,N2,nthreads,jbthread) \
        private(ib,i,j,num_up,num_down,div,tmp_res,tmp_1,tmp_2,jb,all_up,all_down, \
                comb2,mythread,sdim_div,sdim_rest,ib_start,ib_end)
        {
          jb = 0;
#ifdef _OPENMP
          mythread = omp_get_thread_num();
#else
          mythread = 0;
#endif
          comb2 = li_2d_allocate(X->Def.Nsite+1,X->Def.Nsite+1);
          //
          // explict loop decomposition is nessesary to fix the asignment to each thread
          //
          sdim_div = X->Check.sdim / nthreads;
          sdim_rest = X->Check.sdim % nthreads;
          if(mythread < sdim_rest){
            ib_start = sdim_div*mythread + mythread;
            ib_end = ib_start + sdim_div + 1;
          }
          else{
            ib_start = sdim_div*mythread + sdim_rest;
            ib_end = ib_start + sdim_div;
          }
          //
          for(ib=ib_start;ib<ib_end;ib++){
            list_jb[ib]=jb;

            i=ib*ihfbit;
            num_up=0;
            for(j=0;j<=N2-2;j+=2){
              div=i & X->Def.Tpow[j];
              div=div/X->Def.Tpow[j];
              num_up+=div;
            }
            num_down=0;
            for(j=1;j<=N2-1;j+=2){
              div=i & X->Def.Tpow[j];
              div=div/X->Def.Tpow[j];
              num_down+=div;
            }

            tmp_res  = X->Def.Nsite%2; // even Ns-> 0, odd Ns -> 1
            all_up   = (X->Def.Nsite+tmp_res)/2;
            all_down = (X->Def.Nsite-tmp_res)/2;

            tmp_1 = Binomial(all_up,X->Def.Nup-num_up,comb2,all_up);
            tmp_2 = Binomial(all_down,X->Def.Ndown-num_down,comb2,all_down);
            jb   += tmp_1*tmp_2;
          }
          free_li_2d_allocate(comb2);
          if(mythread != nthreads-1) jbthread[mythread+1] = jb;
          #pragma omp barrier
          #pragma omp single
          {
            jbthread[0] = 0;
            for(j=1;j<nthreads;j++){
              jbthread[j] += jbthread[j-1];
            }
          }
          for(ib=ib_start;ib<ib_end;ib++){
            list_jb[ib] += jbthread[mythread];
          }
        }//omp parallel
        free_lui_1d_allocate(jbthread);

        //#pragma omp barrier
        TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzMid, "a");
        TimeKeeper(X, cFileNameTimeKeep, cOMPSzMid, "a");

        icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, X) shared(list_1_, list_2_1_, list_2_2_, list_jb)
        for(ib=0;ib<X->Check.sdim;ib++){
          icnt+=omp_sz_hacker(ib,ihfbit,X,list_1_, list_2_1_, list_2_2_, list_jb);
        }
        break;
      }
      else{
        fprintf(stderr, "Error: CalcHS in ModPara file must be 0 or 1 for Hubbard model.");
        return -1;
      }
      
    case HubbardNConserved:
      hacker = X->Def.read_hacker;
      if(hacker==0){
        // this part can not be parallelized
        jb = 0;
        iSpnup=0;
        iMinup=0;
        iAllup=X->Def.Ne;
        if(X->Def.Ne > X->Def.Nsite){
          iMinup = X->Def.Ne-X->Def.Nsite;
          iAllup = X->Def.Nsite;
        }
        for(ib=0;ib<X->Check.sdim;ib++){
          list_jb[ib]=jb;
          i=ib*ihfbit;
          num_up=0;
          for(j=0;j<=N2-2;j+=2){
            div=i & X->Def.Tpow[j];
            div=div/X->Def.Tpow[j];
            num_up+=div;
          }
          num_down=0;
          for(j=1;j<=N2-1;j+=2){
            div=i & X->Def.Tpow[j];
            div=div/X->Def.Tpow[j];
            num_down+=div;
          }
          tmp_res  = X->Def.Nsite%2; // even Ns-> 0, odd Ns -> 1
          all_up   = (X->Def.Nsite+tmp_res)/2;
          all_down = (X->Def.Nsite-tmp_res)/2;

          for(iSpnup=iMinup; iSpnup<= iAllup; iSpnup++){
            tmp_1 = Binomial(all_up, iSpnup-num_up,comb,all_up);
            tmp_2 = Binomial(all_down, X->Def.Ne-iSpnup-num_down,comb,all_down);
            jb   += tmp_1*tmp_2;
          }
        }
        //#pragma omp barrier
        TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzMid, "a");
        TimeKeeper(X, cFileNameTimeKeep, cOMPSzMid, "a");

        icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N2, X) shared(list_1_, list_2_1_, list_2_2_, list_jb) 
        for(ib=0;ib<X->Check.sdim;ib++){
          icnt+=omp_sz(ib,ihfbit, X,list_1_, list_2_1_, list_2_2_, list_jb);
        }
        break;
      }
      else if(hacker==1){
        iMinup=0;
        iAllup=X->Def.Ne;
        if(X->Def.Ne > X->Def.Nsite){
          iMinup = X->Def.Ne-X->Def.Nsite;
          iAllup = X->Def.Nsite;
        }
        jbthread = lui_1d_allocate(nthreads);
        #pragma omp parallel default(none) \
        shared(X,iMinup,iAllup,list_jb,ihfbit,N2,nthreads,jbthread) \
        private(ib,i,j,num_up,num_down,div,tmp_res,tmp_1,tmp_2,iSpnup,jb,all_up,all_down,comb2, \
                mythread,sdim_rest,sdim_div,ib_start,ib_end)
        {
          jb = 0;
          iSpnup=0;
#ifdef _OPENMP
          mythread = omp_get_thread_num();
#else
          mythread = 0;
#endif
          comb2 = li_2d_allocate(X->Def.Nsite+1,X->Def.Nsite+1);
          //
          // explict loop decomposition is nessesary to fix the asignment to each thread
          //
          sdim_div = X->Check.sdim / nthreads;
          sdim_rest = X->Check.sdim % nthreads;
          if(mythread < sdim_rest){
            ib_start = sdim_div*mythread + mythread;
            ib_end = ib_start + sdim_div + 1;
          }
          else{
            ib_start = sdim_div*mythread + sdim_rest;
            ib_end = ib_start + sdim_div;
          }
          //
          for(ib=ib_start;ib<ib_end;ib++){
            list_jb[ib]=jb;
            i=ib*ihfbit;
            num_up=0;
            for(j=0;j<=N2-2;j+=2){
              div=i & X->Def.Tpow[j];
              div=div/X->Def.Tpow[j];
              num_up+=div;
            }
            num_down=0;
            for(j=1;j<=N2-1;j+=2){
              div=i & X->Def.Tpow[j];
              div=div/X->Def.Tpow[j];
              num_down+=div;
            }
            tmp_res  = X->Def.Nsite%2; // even Ns-> 0, odd Ns -> 1
            all_up   = (X->Def.Nsite+tmp_res)/2;
            all_down = (X->Def.Nsite-tmp_res)/2;

            for(iSpnup=iMinup; iSpnup<= iAllup; iSpnup++){
              tmp_1 = Binomial(all_up, iSpnup-num_up,comb2,all_up);
              tmp_2 = Binomial(all_down, X->Def.Ne-iSpnup-num_down,comb2,all_down);
              jb   += tmp_1*tmp_2;
            }
          }
          free_li_2d_allocate(comb2);
          if(mythread != nthreads-1) jbthread[mythread+1] = jb;
          #pragma omp barrier
          #pragma omp single
          {
            jbthread[0] = 0;
            for(j=1;j<nthreads;j++){
              jbthread[j] += jbthread[j-1];
            }
          }
          for(ib=ib_start;ib<ib_end;ib++){
            list_jb[ib] += jbthread[mythread];
          }
        }//omp parallel
        free_lui_1d_allocate(jbthread);
        //#pragma omp barrier
        TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzMid, "a");
        TimeKeeper(X, cFileNameTimeKeep, cOMPSzMid, "a");

        icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N2, X) shared(list_1_, list_2_1_, list_2_2_, list_jb) 
        for(ib=0;ib<X->Check.sdim;ib++){
          icnt+=omp_sz_hacker(ib,ihfbit, X,list_1_, list_2_1_, list_2_2_, list_jb);
        }

        break;
      }
      else{
        fprintf(stderr, "Error: CalcHS in ModPara file must be 0 or 1 for Hubbard model.");
        return -1;
      }
      
    case Kondo:
      // this part can not be parallelized
      N_all_up   = X->Def.Nup;
      N_all_down = X->Def.Ndown;
      fprintf(stdoutMPI, cStateNupNdown, N_all_up,N_all_down);

      jb = 0;
      num_loc=0;
      for(j=X->Def.Nsite/2; j< X->Def.Nsite ;j++){// counting localized # of spins
        if(X->Def.LocSpn[j] != ITINERANT){
          num_loc += 1;
        }
      }

      for(ib=0;ib<X->Check.sdim;ib++){ //sdim = 2^(N/2)
        list_jb[ib] = jb;
        i           = ib*ihfbit; // ihfbit=pow(2,((Nsite+1)/2))
        num_up      = 0;
        num_down    = 0;
        icheck_loc  = 1;

        for(j=X->Def.Nsite/2; j< X->Def.Nsite ;j++){
          div_up    = i & X->Def.Tpow[2*j];
          div_up    = div_up/X->Def.Tpow[2*j];
          div_down  = i & X->Def.Tpow[2*j+1];
          div_down  = div_down/X->Def.Tpow[2*j+1];
          if(X->Def.LocSpn[j] == ITINERANT){
            num_up   += div_up;        
            num_down += div_down;  
          }else{    
            num_up   += div_up;     
            num_down += div_down;
            if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){ // odd site
              icheck_loc= icheck_loc;
              ihfSpinDown=div_down;
              if(div_down ==0){
                num_up += 1;
              }
            }else{
              icheck_loc   = icheck_loc*(div_up^div_down);// exclude empty or doubly occupied site
            }
          }
        }

        if(icheck_loc == 1){ // itinerant of local spins without holon or doublon
          tmp_res  = X->Def.Nsite%2; // even Ns-> 0, odd Ns -> 1
          all_loc =  X->Def.NLocSpn-num_loc; // # of local spins
          all_up   = (X->Def.Nsite+tmp_res)/2-all_loc;
          all_down = (X->Def.Nsite-tmp_res)/2-all_loc;
          if(X->Def.Nsite%2==1 && X->Def.LocSpn[X->Def.Nsite/2] != ITINERANT){
            all_up   = (X->Def.Nsite)/2-all_loc;
            all_down = (X->Def.Nsite)/2-all_loc;
          }

          for(num_loc_up=0; num_loc_up <= all_loc; num_loc_up++){
            tmp_1 = Binomial(all_loc, num_loc_up, comb, all_loc);
            if( X->Def.Nsite%2==1 && X->Def.LocSpn[X->Def.Nsite/2] != ITINERANT){
              if(ihfSpinDown !=0){
                tmp_2 = Binomial(all_up, X->Def.Nup-num_up-num_loc_up, comb, all_up);
                tmp_3 = Binomial(all_down, X->Def.Ndown-num_down-(all_loc-num_loc_up), comb, all_down);
              }
              else{
                tmp_2 = Binomial(all_up, X->Def.Nup-num_up-num_loc_up, comb, all_up);
                tmp_3 = Binomial(all_down, X->Def.Ndown-num_down-(all_loc-num_loc_up), comb, all_down);
              }
            }
            else{
              tmp_2 = Binomial(all_up, X->Def.Nup-num_up-num_loc_up, comb, all_up);
              tmp_3 = Binomial(all_down, X->Def.Ndown-num_down-(all_loc-num_loc_up), comb, all_down);
            }
            jb   += tmp_1*tmp_2*tmp_3;
          }
        }

      }
      //#pragma omp barrier
      TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzMid, "a");
      TimeKeeper(X, cFileNameTimeKeep, cOMPSzMid, "a");

      hacker = X->Def.read_hacker;
      if(hacker==0){
        icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N2, X) shared(list_1_, list_2_1_, list_2_2_, list_jb)
        for(ib=0;ib<X->Check.sdim;ib++){
          icnt+=omp_sz_Kondo(ib,ihfbit, X, list_1_, list_2_1_, list_2_2_, list_jb);
        }
      }else if(hacker==1){
        icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N2, X) shared(list_1_, list_2_1_, list_2_2_, list_jb)
        for(ib=0;ib<X->Check.sdim;ib++){
          icnt+=omp_sz_Kondo_hacker(ib,ihfbit, X, list_1_, list_2_1_, list_2_2_, list_jb);
        } 
      }
      break;

    case Spin:
      // this part can not be parallelized
      if(X->Def.iFlgGeneralSpin==FALSE){
        hacker = X->Def.read_hacker;
        //printf(" rank=%d:Ne=%ld ihfbit=%ld sdim=%ld\n", myrank,X->Def.Ne,ihfbit,X->Check.sdim);
        // using hacker's delight only + no open mp 
        if(hacker        ==  -1){
          icnt    = 1;
          tmp_pow = 1;
          tmp_i   = 0;
          jb      = 0;
          ja      = 0;
          while(tmp_pow < X->Def.Tpow[X->Def.Ne]){
            tmp_i   += tmp_pow;
            tmp_pow  = tmp_pow*2;
          }
          //printf("DEBUG: %ld %ld %ld %ld\n",tmp_i,X->Check.sdim,X->Def.Tpow[X->Def.Ne],X->Def.Nsite);
          if(X->Def.Nsite%2==0){
            max_tmp_i = X->Check.sdim*X->Check.sdim;
          }else{
            max_tmp_i = X->Check.sdim*X->Check.sdim*2-1;
          }  
          while(tmp_i<max_tmp_i){
            list_1_[icnt]=tmp_i;
           
            ia= tmp_i & irght;
            ib= tmp_i & ilft;
            ib= ib/ihfbit; 
            if(ib==ibpatn){
              ja=ja+1;
            }else{
              ibpatn = ib;
              ja     = 1;
              jb     = icnt-1;
            }
            
            list_2_1_[ia] = ja+1;
            list_2_2_[ib] = jb+1;
            tmp_j = snoob(tmp_i);
            tmp_i =        tmp_j;
            icnt        +=  1;
          }
          icnt = icnt-1;
          // old version + hacker's delight
        }else if(hacker  ==  1){
          jb = 0;
          for(ib=0;ib<X->Check.sdim;ib++){
            list_jb[ib]=jb;
            i=ib*ihfbit;
            num_up=0;
            for(j=0;j<N; j++){
              div_up = i & X->Def.Tpow[j];
              div_up = div_up/X->Def.Tpow[j];
              num_up +=div_up;
            }
            all_up   = (X->Def.Nsite+1)/2;
            tmp_1 = Binomial(all_up,X->Def.Ne-num_up,comb,all_up);
            jb   += tmp_1;
          }
          //#pragma omp barrier
          TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzMid, "a");
          TimeKeeper(X, cFileNameTimeKeep, cOMPSzMid, "a");

          icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N, X, list_1_, list_2_1_, list_2_2_, list_jb)
          for(ib=0;ib<X->Check.sdim;ib++){
            icnt+=omp_sz_spin_hacker(ib,ihfbit,N,X, list_1_, list_2_1_, list_2_2_, list_jb);
          }
          //printf(" rank=%d ib=%ld:Ne=%d icnt=%ld :idim_max=%ld N=%d\n", myrank,ib,X->Def.Ne,icnt,X->Check.idim_max,N);
          // old version
        }else if(hacker  ==  0){
          jb = 0;
          for(ib=0;ib<X->Check.sdim;ib++){
            list_jb[ib]=jb;
            i=ib*ihfbit;
            num_up=0;
            for(j=0;j<N; j++){
              div_up = i & X->Def.Tpow[j];
              div_up = div_up/X->Def.Tpow[j];
              num_up +=div_up;
            }
            all_up   = (X->Def.Nsite+1)/2;
            tmp_1 = Binomial(all_up,X->Def.Ne-num_up,comb,all_up);
            jb   += tmp_1;
          }
          //#pragma omp barrier
          TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzMid, "a");
          TimeKeeper(X, cFileNameTimeKeep, cOMPSzMid, "a");

          icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ihfbit, N, X) shared(list_1_, list_2_1_, list_2_2_, list_jb)
          for(ib=0;ib<X->Check.sdim;ib++){
            icnt+=omp_sz_spin(ib,ihfbit,N,X,list_1_, list_2_1_, list_2_2_, list_jb);
          }
        }
        else{
          fprintf(stderr, "Error: CalcHS in ModPara file must be -1 or 0 or 1 for Spin model.");
          return -1;
        }
      }else{
        unsigned int Max2Sz=0;
        unsigned int irghtsite=1;
        long unsigned int itmpSize=1;
        int i2Sz=0;
        for(j=0; j<X->Def.Nsite; j++){
          itmpSize *= X->Def.SiteToBit[j];
          if(itmpSize==ihfbit){
            break;
          }
          irghtsite++;
        }
        for(j=0; j<X->Def.Nsite; j++){
          Max2Sz += X->Def.LocSpn[j];
        }

        HilbertNumToSz = lui_1d_allocate(2*Max2Sz+1);
        for(ib=0; ib<2*Max2Sz+1; ib++){
          HilbertNumToSz[ib]=0;
        }

        for(ib =0; ib<ihfbit; ib++){
          i2Sz=0;
          for(j=1; j<= irghtsite; j++){
            i2Sz += GetLocal2Sz(j,ib, X->Def.SiteToBit, X->Def.Tpow);
          }
          list_2_1_Sz[ib]=i2Sz;
          HilbertNumToSz[i2Sz+Max2Sz]++;
        }
        jb = 0;
        long unsigned int ilftdim=(X->Def.Tpow[X->Def.Nsite-1]*X->Def.SiteToBit[X->Def.Nsite-1])/ihfbit;
        for(ib=0;ib<ilftdim;ib++){
          list_jb[ib]=jb;
          i2Sz=0;
          for(j=1;j<=(N-irghtsite); j++){
            i2Sz += GetLocal2Sz(j+irghtsite,ib*ihfbit, X->Def.SiteToBit, X->Def.Tpow);
          }
          list_2_2_Sz[ib]=i2Sz;
          if((X->Def.Total2Sz- i2Sz +(int)Max2Sz)>=0 && (X->Def.Total2Sz- i2Sz) <= (int)Max2Sz){
            jb += HilbertNumToSz[X->Def.Total2Sz- i2Sz +Max2Sz];
          }
        }

        TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzMid, "a");
        TimeKeeper(X, cFileNameTimeKeep, cOMPSzMid, "a");

        icnt = 0;
#pragma omp parallel for default(none) reduction(+:icnt) private(ib) firstprivate(ilftdim, ihfbit,  X)  shared(list_1_, list_2_1_, list_2_2_, list_2_1_Sz, list_2_2_Sz,list_jb)
        for(ib=0;ib<ilftdim; ib++){
          icnt+=omp_sz_GeneralSpin(ib,ihfbit,X, list_1_, list_2_1_, list_2_2_, list_2_1_Sz, list_2_2_Sz,list_jb);
        }

        free_lui_1d_allocate(HilbertNumToSz);
      }
      
      break;
    default:
      exitMPI(-1);
       
    }    
    i_max=icnt;
    //fprintf(stdoutMPI, "Debug: Xicnt=%ld \n",icnt);
    TimeKeeper(X, cFileNameSzTimeKeep, cOMPSzFinish, "a");
    TimeKeeper(X, cFileNameTimeKeep, cOMPSzFinish, "a");

  }

  if(X->Def.iFlgCalcSpec == CALCSPEC_NOT){
    if(X->Def.iCalcModel==HubbardNConserved){
      X->Def.iCalcModel=Hubbard;
    }
  }
  
  //Error message
  //i_max=i_max+1;
  if(i_max!=X->Check.idim_max){
    fprintf(stderr, "%s", cErrSz);
    fprintf(stderr, cErrSz_ShowDim, i_max, X->Check.idim_max);
    strcpy(sdt_err,cFileNameErrorSz);
    if(childfopenMPI(sdt_err,"a",&fp_err)!=0){
      exitMPI(-1);
    }
    fprintf(fp_err, "%s",cErrSz_OutFile);
    fclose(fp_err);
    exitMPI(-1);
  }

  free_li_2d_allocate(comb);
  }
  fprintf(stdoutMPI, "%s", cProEndCalcSz);

    free(list_jb);
    if(X->Def.iFlgGeneralSpin==TRUE){
        free(list_2_1_Sz);
        free(list_2_2_Sz);
    }
  return 0;    
}

/** 
 * 
 * @file   sz.c
 * @brief  calculating binomial coefficients
 * 
 * @param[in] n      n for @f$_nC_k = \frac{n!}{(n-k)!k!}@f$
 * @param[in] k      k for @f$_nC_k = \frac{n!}{(n-k)!k!}@f$
 * @param[out] comb   binomial coefficients @f$_nC_k@f$
 * @param[in] Nsite  # of sites
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
long int Binomial(int n,int k,long int **comb,int Nsite){
  // nCk, Nsite=max(n)
  int tmp_i,tmp_j;

  if(n==0 && k==0){
    return 1;
  } 
  else if(n<0 || k<0 || n<k){
    return 0;
  }
  
  for(tmp_i=0;tmp_i<=Nsite;tmp_i++){
    for(tmp_j=0;tmp_j<=Nsite;tmp_j++){
      comb[tmp_i][tmp_j] = 0;
    }
  }

  comb[0][0] = 1;
  comb[1][0] = 1;
  comb[1][1] = 1;
  for(tmp_i=2;tmp_i<=n;tmp_i++){
    for(tmp_j=0;tmp_j<=tmp_i;tmp_j++){
      if(tmp_j==0){
        comb[tmp_i][tmp_j] = 1;
      }else if(tmp_j==tmp_i){
        comb[tmp_i][tmp_j] = 1;
      }else{
        comb[tmp_i][tmp_j] = comb[tmp_i-1][tmp_j-1]+comb[tmp_i-1][tmp_j];
      }
    }
  }
  return comb[n][k];
}

/** 
 * @brief calculating restricted Hilbert space for Hubbard systems
 * 
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[in] X
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia) 
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int omp_sz(
                 long unsigned int ib,    //!<[in]
                 long unsigned int ihfbit, //!<[in]
                 struct BindStruct *X,     //!<[in]
                 long unsigned int *list_1_, //!<[out]
                 long unsigned int *list_2_1_,//!<[out]
                 long unsigned int *list_2_2_,//!<[out]
                 long unsigned int *list_jb_ //!<[in]
                 )
{
  long unsigned int i,j; 
  long unsigned int ia,ja,jb;
  long unsigned int div_down, div_up;
  long unsigned int num_up,num_down;
  long unsigned int tmp_num_up,tmp_num_down;
    
  jb = list_jb_[ib];
  i  = ib*ihfbit;
    
  num_up   = 0;
  num_down = 0;
  for(j=0;j< X->Def.Nsite ;j++){
    div_up    = i & X->Def.Tpow[2*j];
    div_up    = div_up/X->Def.Tpow[2*j];
    div_down  = i & X->Def.Tpow[2*j+1];
    div_down  = div_down/X->Def.Tpow[2*j+1];
    num_up += div_up;
    num_down += div_down;
  }
  
  ja=1;
  tmp_num_up   = num_up;
  tmp_num_down = num_down;

  if(X->Def.iCalcModel==Hubbard){
    for(ia=0;ia<X->Check.sdim;ia++){
      i=ia;
      num_up =  tmp_num_up;
      num_down =  tmp_num_down;
      
      for(j=0;j<X->Def.Nsite;j++){
        div_up    = i & X->Def.Tpow[2*j];
        div_up    = div_up/X->Def.Tpow[2*j];
        div_down  = i & X->Def.Tpow[2*j+1];
        div_down  = div_down/X->Def.Tpow[2*j+1];
        num_up += div_up;
        num_down += div_down;
      }
      if(num_up == X->Def.Nup && num_down == X->Def.Ndown){
        list_1_[ja+jb]=ia+ib*ihfbit;
        list_2_1_[ia]=ja+1;
        list_2_2_[ib]=jb+1;
        ja+=1;
      } 
    }
  }
  else if(X->Def.iCalcModel==HubbardNConserved){
    for(ia=0;ia<X->Check.sdim;ia++){
      i=ia;
      num_up =  tmp_num_up;
      num_down =  tmp_num_down;
      
      for(j=0;j<X->Def.Nsite;j++){
        div_up    = i & X->Def.Tpow[2*j];
        div_up    = div_up/X->Def.Tpow[2*j];
        div_down  = i & X->Def.Tpow[2*j+1];
        div_down  = div_down/X->Def.Tpow[2*j+1];
        num_up += div_up;
        num_down += div_down;
      }
      if( (num_up+num_down) == X->Def.Ne){
        list_1_[ja+jb]=ia+ib*ihfbit;
        list_2_1_[ia]=ja+1;
        list_2_2_[ib]=jb+1;
        ja+=1;
      } 
    }  
  }
  ja=ja-1;    
  return ja; 
}
/** 
 * @brief efficient version of calculating restricted Hilbert space for Hubbard systems  using snoob
 * details of snoob is found in S.H. Warren, Hacker’s Delight, second ed., Addison-Wesley, ISBN: 0321842685, 2012.
 *
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[in] X
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia) 
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 */
int omp_sz_hacker(long unsigned int ib,
                        long unsigned int ihfbit,
                        struct BindStruct *X,
                        long unsigned int *list_1_,
                        long unsigned int *list_2_1_,
                        long unsigned int *list_2_2_,
                        long unsigned int *list_jb_
                        )
{
  long unsigned int i,j; 
  long unsigned int ia,ja,jb;
  long unsigned int div_down, div_up;
  long unsigned int num_up,num_down;
  long unsigned int tmp_num_up,tmp_num_down;
    
  jb = list_jb_[ib];
  i  = ib*ihfbit;
    
  num_up   = 0;
  num_down = 0;
  for(j=0;j< X->Def.Nsite ;j++){
    div_up    = i & X->Def.Tpow[2*j];
    div_up    = div_up/X->Def.Tpow[2*j];
    div_down  = i & X->Def.Tpow[2*j+1];
    div_down  = div_down/X->Def.Tpow[2*j+1];
    num_up += div_up;
    num_down += div_down;
  }
  
  ja=1;
  tmp_num_up   = num_up;
  tmp_num_down = num_down;

  if(X->Def.iCalcModel==Hubbard){
    if(tmp_num_up <= X->Def.Nup && tmp_num_down <= X->Def.Ndown){ //do not exceed Nup and Ndown
      ia = X->Def.Tpow[X->Def.Nup+X->Def.Ndown-tmp_num_up-tmp_num_down]-1;
      if(ia < X->Check.sdim){
        num_up   =  tmp_num_up;
        num_down =  tmp_num_down;
        for(j=0;j<X->Def.Nsite;j++){
          div_up    = ia & X->Def.Tpow[2*j];
          div_up    = div_up/X->Def.Tpow[2*j];
          div_down  = ia & X->Def.Tpow[2*j+1];
          div_down  = div_down/X->Def.Tpow[2*j+1];
          num_up   += div_up;
          num_down += div_down;
        }
        if(num_up == X->Def.Nup && num_down == X->Def.Ndown){
          list_1_[ja+jb]=ia+ib*ihfbit;
          list_2_1_[ia]=ja+1;
          list_2_2_[ib]=jb+1;
          ja+=1;
        }
        if(ia!=0){
          ia = snoob(ia);
          while(ia < X->Check.sdim){
            num_up   =  tmp_num_up;
            num_down =  tmp_num_down;
            for(j=0;j<X->Def.Nsite;j++){
              div_up    = ia & X->Def.Tpow[2*j];
              div_up    = div_up/X->Def.Tpow[2*j];
              div_down  = ia & X->Def.Tpow[2*j+1];
              div_down  = div_down/X->Def.Tpow[2*j+1];
              num_up   += div_up;
              num_down += div_down;
            }
            if(num_up == X->Def.Nup && num_down == X->Def.Ndown){
              list_1_[ja+jb]=ia+ib*ihfbit;
              list_2_1_[ia]=ja+1;
              list_2_2_[ib]=jb+1;
              ja+=1;
            }
            ia = snoob(ia);
          }
        } 
      } 
    }
  }
  else if(X->Def.iCalcModel==HubbardNConserved){
    if(tmp_num_up+tmp_num_down <= X->Def.Ne){ //do not exceed Ne
      ia = X->Def.Tpow[X->Def.Ne-tmp_num_up-tmp_num_down]-1;
      if(ia < X->Check.sdim){
        list_1_[ja+jb]=ia+ib*ihfbit;
        list_2_1_[ia]=ja+1;
        list_2_2_[ib]=jb+1;
        ja+=1;
        if(ia!=0){
          ia = snoob(ia);
          while(ia < X->Check.sdim){
            list_1_[ja+jb]=ia+ib*ihfbit;
            list_2_1_[ia]=ja+1;
            list_2_2_[ib]=jb+1;
            ja+=1;
            ia = snoob(ia);
          }
        } 
      }  
    }
  }
  ja=ja-1;    
  return ja; 
}

/** 
 * @brief calculating restricted Hilbert space for Kondo systems
 *
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[in] X
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia) 
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 */
int omp_sz_Kondo(
                       long unsigned int ib,        //[in]
                       long unsigned int ihfbit,    //[in]
                       struct BindStruct *X,        //[in]
                       long unsigned int *list_1_,  //[out]
                       long unsigned int *list_2_1_,//[out]
                       long unsigned int *list_2_2_,//[out]
                       long unsigned int *list_jb_  //[in]
                       )
{
  long unsigned int i,j; 
  long unsigned int ia,ja,jb;
  long unsigned int div_down, div_up;
  long unsigned int num_up,num_down;
  long unsigned int tmp_num_up,tmp_num_down;
  int icheck_loc;
    
  jb = list_jb_[ib];
  i  = ib*ihfbit;
    
  num_up   = 0;
  num_down = 0;
  icheck_loc=1;
  for(j=X->Def.Nsite/2; j< X->Def.Nsite ;j++){
    div_up    = i & X->Def.Tpow[2*j];
    div_up    = div_up/X->Def.Tpow[2*j];
    div_down  = i & X->Def.Tpow[2*j+1];
    div_down  = div_down/X->Def.Tpow[2*j+1];

    if(X->Def.LocSpn[j] == ITINERANT){
      num_up   += div_up;        
      num_down += div_down;  
    }else{    
      num_up   += div_up;        
      num_down += div_down;
      if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
        icheck_loc= icheck_loc;
      }
      else{
        icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
      }
    }
  }
  
  ja=1;
  tmp_num_up   = num_up;
  tmp_num_down = num_down;
  if(icheck_loc ==1){
    for(ia=0;ia<X->Check.sdim;ia++){
      i=ia;
      num_up =  tmp_num_up;
      num_down =  tmp_num_down;
      icheck_loc=1;
      for(j=0;j<(X->Def.Nsite+1)/2;j++){
        div_up    = i & X->Def.Tpow[2*j];
        div_up    = div_up/X->Def.Tpow[2*j];
        div_down  = i & X->Def.Tpow[2*j+1];
        div_down  = div_down/X->Def.Tpow[2*j+1];

        if(X->Def.LocSpn[j] ==  ITINERANT){
          num_up   += div_up;        
          num_down += div_down;  
        }else{    
          num_up   += div_up;        
          num_down += div_down;  
          if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
            icheck_loc= icheck_loc;
          }
          else{
            icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
          }
        }
      }
      
      if(icheck_loc == 1 && X->Def.LocSpn[X->Def.Nsite/2] != ITINERANT && X->Def.Nsite%2==1){
        div_up    = ia & X->Def.Tpow[X->Def.Nsite-1];
        div_up    = div_up/X->Def.Tpow[X->Def.Nsite-1];
        div_down  = (ib*ihfbit) & X->Def.Tpow[X->Def.Nsite];
        div_down  = div_down/X->Def.Tpow[X->Def.Nsite];
        icheck_loc= icheck_loc*(div_up^div_down);
      }
      
      if(num_up == X->Def.Nup && num_down == X->Def.Ndown && icheck_loc==1){
        list_1_[ja+jb]=ia+ib*ihfbit;
        /*
        list_2_1_[ia]=ja;
        list_2_2_[ib]=jb;
         */
        list_2_1_[ia]=ja+1;
        list_2_2_[ib]=jb+1;
        //printf("DEBUG: rank=%d, list_1[%d]=%d, list_2_1_[%d]=%d, list_2_2_[%d]=%d\n", myrank, ja+jb, list_1_[ja+jb], ia, list_2_1[ia], ib, list_2_2[ib]);
        ja+=1;
      }
    }
  }
  ja=ja-1;    
  return ja; 
}
/** 
 * @brief calculating restricted Hilbert space for Kondo-GC systems
 *
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[in] X
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia) 
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 */
int omp_sz_Kondo_hacker(
                       long unsigned int ib,
                       long unsigned int ihfbit,
                       struct BindStruct *X,
                       long unsigned int *list_1_,
                       long unsigned int *list_2_1_,
                       long unsigned int *list_2_2_,
                       long unsigned int *list_jb_
                       )
{
  long unsigned int i,j; 
  long unsigned int ia,ja,jb;
  long unsigned int div_down, div_up;
  long unsigned int num_up,num_down;
  long unsigned int tmp_num_up,tmp_num_down;
  int icheck_loc;
    
  jb = list_jb_[ib];
  i  = ib*ihfbit;
    
  num_up   = 0;
  num_down = 0;
  icheck_loc=1;
  for(j=X->Def.Nsite/2; j< X->Def.Nsite ;j++){
    div_up    = i & X->Def.Tpow[2*j];
    div_up    = div_up/X->Def.Tpow[2*j];
    div_down  = i & X->Def.Tpow[2*j+1];
    div_down  = div_down/X->Def.Tpow[2*j+1];

    if(X->Def.LocSpn[j] == ITINERANT){
      num_up   += div_up;        
      num_down += div_down;  
    }else{    
      num_up   += div_up;        
      num_down += div_down;
      if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
        icheck_loc= icheck_loc;
      }
      else{
        icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubly occupied site
      }
    }
  }
//[s] get ja  
  ja           = 1;
  tmp_num_up   = num_up;
  tmp_num_down = num_down;
  if(icheck_loc ==1){
    //for(ia=0;ia<X->Check.sdim;ia++){
    ia = X->Def.Tpow[X->Def.Nup+X->Def.Ndown-tmp_num_up-tmp_num_down]-1;
    //ia = 1;
    //if(ia < X->Check.sdim && ia!=0){
    //ia = snoob(ia);
    while(ia < X->Check.sdim && ia!=0){
    // for(ia=0;ia<X->Check.sdim;ia++){
        //[s] proceed ja
        i        = ia;
        num_up   =  tmp_num_up;
        num_down =  tmp_num_down;
        icheck_loc=1;
        for(j=0;j<(X->Def.Nsite+1)/2;j++){
          div_up    = i & X->Def.Tpow[2*j];
          div_up    = div_up/X->Def.Tpow[2*j];
          div_down  = i & X->Def.Tpow[2*j+1];
          div_down  = div_down/X->Def.Tpow[2*j+1];

          if(X->Def.LocSpn[j] ==  ITINERANT){
            num_up   += div_up;        
            num_down += div_down;  
          }else{    
            num_up   += div_up;        
            num_down += div_down;  
            if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
              icheck_loc= icheck_loc;
            }
            else{
              icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
            }
          }
        }
        
        if(icheck_loc == 1 && X->Def.LocSpn[X->Def.Nsite/2] != ITINERANT && X->Def.Nsite%2==1){
          div_up    = ia & X->Def.Tpow[X->Def.Nsite-1];
          div_up    = div_up/X->Def.Tpow[X->Def.Nsite-1];
          div_down  = (ib*ihfbit) & X->Def.Tpow[X->Def.Nsite];
          div_down  = div_down/X->Def.Tpow[X->Def.Nsite];
          icheck_loc= icheck_loc*(div_up^div_down);
        }
        
        if(num_up == X->Def.Nup && num_down == X->Def.Ndown && icheck_loc==1){
          //printf("ia=%ud ja=%ud \n",ia,ja);
          list_1_[ja+jb]=ia+ib*ihfbit;
          list_2_1_[ia]=ja+1;
          list_2_2_[ib]=jb+1;
          ja+=1;
        }
        ia = snoob(ia);
        //[e] proceed ja
        //ia+=1;
      //}
    }
  }
//[e] get ja
  ja=ja-1;    
  return ja; 
}


/** 
 * 
 * 
 * @param ib 
 * @param ihfbit 
 * @param N2 
 * @param X 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int omp_sz_KondoGC(
                         long unsigned int ib,  //!<[in]
                         long unsigned int ihfbit,//!<[in]
                         struct BindStruct *X,    //!<[in]
                         long unsigned int *list_1_, //!<[out]
                         long unsigned int *list_2_1_,//!<[out]
                         long unsigned int *list_2_2_,//!<[out]
                         long unsigned int *list_jb_//!<[in]
                         )
{
  long unsigned int i,j; 
  long unsigned int ia,ja,jb;
  long unsigned int div_down, div_up;
  int icheck_loc;
    
  jb = list_jb_[ib];
  i  = ib*ihfbit;
  icheck_loc=1;
  for(j=X->Def.Nsite/2; j< X->Def.Nsite ;j++){
    div_up    = i & X->Def.Tpow[2*j];
    div_up    = div_up/X->Def.Tpow[2*j];
    div_down  = i & X->Def.Tpow[2*j+1];
    div_down  = div_down/X->Def.Tpow[2*j+1];
    if(X->Def.LocSpn[j] !=  ITINERANT){
      if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
        icheck_loc= icheck_loc;
      }
      else{
        icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
      }
    }
  }

  ja=1;
  if(icheck_loc ==1){
    for(ia=0;ia<X->Check.sdim;ia++){
      i=ia;
      icheck_loc =1;
      for(j=0;j<(X->Def.Nsite+1)/2;j++){
        div_up    = i & X->Def.Tpow[2*j];
        div_up    = div_up/X->Def.Tpow[2*j];
        div_down  = i & X->Def.Tpow[2*j+1];
        div_down  = div_down/X->Def.Tpow[2*j+1];
        if(X->Def.LocSpn[j] !=  ITINERANT){
          if(X->Def.Nsite%2==1 && j==(X->Def.Nsite/2)){
            icheck_loc= icheck_loc;
          }
          else{
            icheck_loc   = icheck_loc*(div_up^div_down);// exclude doubllly ocupited site
          }
        }
      }

      if(icheck_loc == 1 && X->Def.LocSpn[X->Def.Nsite/2] != ITINERANT && X->Def.Nsite%2==1){
        div_up    = ia & X->Def.Tpow[X->Def.Nsite-1];
        div_up    = div_up/X->Def.Tpow[X->Def.Nsite-1];
        div_down  = (ib*ihfbit) & X->Def.Tpow[X->Def.Nsite];
        div_down  = div_down/X->Def.Tpow[X->Def.Nsite];
        icheck_loc= icheck_loc*(div_up^div_down);
      }
      
      if(icheck_loc==1){
        list_1_[ja+jb]=ia+ib*ihfbit;
        list_2_1_[ia]=ja+1;
        list_2_2_[ib]=jb+1;
        ja+=1;
      }
    }
  }
  ja=ja-1;
    
  return ja; 
}

/** 
 * @brief calculating restricted Hilbert space for spin-1/2 systems
 *
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[in] X          
 * @param[in] N ???
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia)
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 */
int omp_sz_spin(
                      long unsigned int ib, 
                      long unsigned int ihfbit,
                      unsigned int N,
                      struct BindStruct *X,
                      long unsigned int *list_1_,
                      long unsigned int *list_2_1_,
                      long unsigned int *list_2_2_,
                      long unsigned int *list_jb_
                      )
{
  long unsigned int i,j,div; 
  long unsigned int ia,ja,jb;
  long unsigned int num_up;
  unsigned int tmp_num_up;
  
  jb = list_jb_[ib];
  i  = ib*ihfbit;
  num_up=0;
  for(j=0;j<N;j++){
    div=i & X->Def.Tpow[j];
    div=div/X->Def.Tpow[j];
    num_up+=div;
  }
  ja=1;
  tmp_num_up   = num_up;
  
  for(ia=0;ia<ihfbit;ia++){
    i=ia;
    num_up =  tmp_num_up;
    for(j=0;j<N;j++){
      div=i & X->Def.Tpow[j];
      div=div/X->Def.Tpow[j];
      num_up+=div;
    }

    if(num_up == X->Def.Ne){
      list_1_[ja+jb]=ia+ib*ihfbit;
      list_2_1_[ia]=ja+1;
      list_2_2_[ib]=jb+1;
      ja+=1;
    } 
  }
  ja=ja-1;
  return ja; 
}


/** 
 * @brief efficient version of calculating restricted Hilbert space for spin-1/2 systems 
 * details of snoob is found in S.H. Warren, Hacker’s Delight, second ed., Addison-Wesley, ISBN: 0321842685, 2012.
 *
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[in] X          
 * @param[in] N ???
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia)
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 */
int omp_sz_spin_hacker(
                             long unsigned int ib, 
                             long unsigned int ihfbit,
                             unsigned int N,
                             struct BindStruct *X,
                             long unsigned int *list_1_,
                             long unsigned int *list_2_1_,
                             long unsigned int *list_2_2_,
                             long unsigned int *list_jb_
                             )
{
  long unsigned int i,j,div; 
  long unsigned int ia,ja,jb;
  long unsigned int num_up;
  unsigned int tmp_num_up;
  
  jb = list_jb_[ib];
  i  = ib*ihfbit;
  num_up=0;
  for(j=0;j<N;j++){
    div=i & X->Def.Tpow[j];
    div=div/X->Def.Tpow[j];
    num_up+=div;
  }
  ja=1;
  tmp_num_up   = num_up;
  
  // using hacker's delight
  if(tmp_num_up<=X->Def.Ne && (X->Def.Ne-tmp_num_up)<= X->Def.Nsite-1){ // do not exceed Ne
    ia = X->Def.Tpow[X->Def.Ne-tmp_num_up]-1;
    if(ia<ihfbit ){          // do not exceed Ne
      list_1_[ja+jb] = ia+ib*ihfbit;
      list_2_1_[ia]  = ja+1;
      list_2_2_[ib]  = jb+1;
      ja           += 1;

      if(ia!=0){
        ia = snoob(ia);
        while(ia < ihfbit){
          //fprintf(stdoutMPI, " X: ia= %ld ia=%ld \n", ia,ia);
          list_1_[ja+jb]    = ia+ib*ihfbit;
          list_2_1_[ia]     = ja+1;
          list_2_2_[ib]     = jb+1;
          ja+=1;
          ia = snoob(ia);
        }
      }
    }
  }
  ja=ja-1;
  return ja; 
}

/** 
 * @brief calculating restricted Hilbert space for general spin systems (S>1/2)
 *
 * @param[in] ib   upper half bit of i    
 * @param[in] ihfbit 2^(Ns/2) 
 * @param[in] X          
 * @param[out] list_1_    list_1_[icnt] = i : i is divided into ia and ib (i=ib*ihfbit+ia) 
 * @param[out] list_2_1_  list_2_1_[ib] = jb  
 * @param[out] list_2_2_  list_2_2_[ia] = ja  : icnt=jb+ja
 * @param[out] list_2_1_Sz_  
 * @param[out] list_2_2_Sz_  
 * @param[in] list_jb_   list_jb_[ib]  = jb  
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 */
int omp_sz_GeneralSpin(
                             long unsigned int ib, 
                             long unsigned int ihfbit,
                             struct BindStruct *X,
                             long unsigned int *list_1_,
                             long unsigned int *list_2_1_,
                             long unsigned int *list_2_2_,
                             long int *list_2_1_Sz_,
                             long int *list_2_2_Sz_,
                             long unsigned int *list_jb_
                             )
{
  long unsigned int ia,ja,jb;  
  int list_2_2_Sz_ib=0;
  int tmp_2Sz=0;
  jb = list_jb_[ib];
  list_2_2_Sz_ib =list_2_2_Sz_[ib];
  ja=1;
  for(ia=0;ia<ihfbit;ia++){
    tmp_2Sz=list_2_1_Sz_[ia]+list_2_2_Sz_ib;
    if(tmp_2Sz == X->Def.Total2Sz){
      list_1_[ja+jb]=ia+ib*ihfbit;
      list_2_1_[ia]=ja+1;
      list_2_2_[ib]=jb+1;
      ja+=1;
    } 
  }
  ja=ja-1;
  return ja; 
}

/** 
 * @brief reading the list of the restricted Hilbert space
 * 
 * @param[in] X 
 * @param[in] irght 
 * @param[in] ilft 
 * @param[in] ihfbit 
 * @param[in] i_max 
 * 
 * @return 
 * @author Takahiro Misawa (The University of Tokyo)
 * @author Kazuyoshi Yoshimi (The University of Tokyo)
 */
int Read_sz
(
 struct BindStruct *X,
 const long unsigned int irght,
 const long unsigned int ilft,
 const long unsigned int ihfbit,
 long unsigned int *i_max
 )
{
  FILE *fp,*fp_err;
  char sdt[D_FileNameMax];
  char buf[D_FileNameMax];
    
  long unsigned int icnt=0; 
  long unsigned int ia,ib;
  long unsigned int ja=0;
  long unsigned int jb=0;
  long unsigned int ibpatn=0;
  long unsigned int dam; 

  TimeKeeper(X,cFileNameSzTimeKeep,cReadSzStart, "a");
  TimeKeeper(X,cFileNameTimeKeep,cReadSzStart, "a");

  switch(X->Def.iCalcModel){
  case Hubbard:
  case HubbardGC:
  case Spin:
  case SpinGC:
    sprintf(sdt,cFileNameListModel, X->Def.Nsite, X->Def.Nup, X->Def.Ndown);
    break;
  case Kondo:
    sprintf(sdt,"ListForKondo_Ns%d_Ncond%d.dat",X->Def.Nsite,X->Def.Ne);
    break;
  }
  if(childfopenMPI(sdt,"r", &fp)!=0){
    exitMPI(-1);
  }  

  if(fp == NULL){
    if(childfopenMPI(cFileNameErrorSz,"a",&fp_err)!=0){
      exitMPI(-1);
    }
    fprintf(fp_err, "%s", cErrSz_NoFile);
    fprintf(stderr, "%s", cErrSz_NoFile);
    fprintf(fp_err, cErrSz_NoFile_Show,sdt);
    fprintf(stderr, cErrSz_NoFile_Show, sdt);
    fclose(fp_err);
  }else{
    while(NULL != fgetsMPI(buf,sizeof(buf),fp)){  
      dam=atol(buf);  
      list_1[icnt]=dam;
            
      ia= dam & irght;
      ib= dam & ilft;
      ib=ib/ihfbit; 
            
      if(ib==ibpatn){
        ja=ja+1;
      }else{
        ibpatn=ib;
        ja=1;
        jb=icnt-1;
      }
            
      list_2_1[ia]=ja;
      list_2_2[ib]=jb;
      icnt+=1;
                
    }
    fclose(fp);
    *i_max=icnt-1;
  }

  TimeKeeper(X, cFileNameSzTimeKeep, cReadSzEnd, "a");
  TimeKeeper(X, cFileNameTimeKeep, cReadSzEnd, "a");

  return 0;
}
