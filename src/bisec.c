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

#include "stdlib.h"
#include "math.h"
#include "bisec.h"

/** 
 * @brief obtain the eigen energies by using bisection method
 * This method is based the TITPACK ver 2.
 * 
 * @param alpha coefficient of tridiagonal matrix in the Lanczos method
 * @param beta coefficient of tridiagonal matrix in the Lanczos method 
 * @param ndim dimension of tridiagonal matrix in the Lanczos method
 * @param E  eigen energies
 * @param ne number of excited states  
 * @param eps precision for bisection method
 * @author Takahiro Misawa (The University of Tokyo)
 */
void bisec(double *alpha, //!<[in]
           double *beta,  //!<[in]
           int ndim,  //!<[in]
           double *E, //!<[out]
           int ne,  //!<[in]
           double eps //!<[in]
){ 

    long int i,j,k;
    int numneg,ipass,i_lim;
    double range,*b2;
    double epsabs;
    double a,b,c,g;


    b2=(double *)malloc((ndim+1)*sizeof(double));
/* initial bound*/
      
      
    range=fabs(*(alpha+1))+fabs(*(beta+1));

    for(k=2;k<=ndim-1;k++){
        if(range<fabs(*(beta+k-1)+fabs(*(alpha+k))+fabs(*(beta+k)))){
            range=fabs(*(beta+k-1)+fabs(*(alpha+k))+fabs(*(beta+k)));
        }
    }
      
    if(range<fabs(*(beta+ndim-1))+fabs(*(alpha+ndim))){
        range=fabs(*(beta+ndim-1))+fabs(*(alpha+ndim));
    }
      
    range=-range;
    b2[1]=0.0;
    
    for(i=2;i<=ndim;i++){
        b2[i]=*(beta+i-1)*(*(beta+i-1));    
    }       
      
    epsabs=fabs(range)*eps;
      
    for(i=1;i<=ne;i++){
        *(E+i)=-range;
        b=range;
    }
    

/* bisection method*/
    for(k=1;k<=ne;k++){
        a=*(E+k);
        
        for(j=1;j<=100;j++){  
                c=(a+b)/2.0;
            if(fabs(a-b)<epsabs){
                break;
            }
            
            numneg=0;
            g=1.0;
            ipass=0;
            for(i=1;i<=ndim;i++){
                if(ipass==0){
                    g=c-*(alpha+i)-b2[i]/g;
                }else if(ipass==1){
                    ipass=2;
                }else{
                    g=c-*(alpha+i);
                    ipass=0;
                }

                if(ipass==0){
                    if(g<=0.0){
                        numneg=numneg+1;
                    }
                    
                    if(fabs(g)<=fabs(b2[i]*epsabs*eps)){
                        ipass=1;
                    }
                }
            }  
                
            numneg=ndim-numneg;
            if(numneg<k){
                b=c;
            }else{
                a=c;
                if(numneg>ne){
                    i_lim=ne;
                }else{
                    i_lim=numneg;
                }
                for(i=k;i<=i_lim;i++){
                    *(E+i)=c;
                }
            }
        }
    }
    free(b2);
}
