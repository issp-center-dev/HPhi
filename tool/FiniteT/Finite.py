import numpy as np
import os
import copy
import math
import cmath

with open("output/Eigenvalue.dat") as f:
     data      = f.read()
     data      = data.split("\n")
     max_eigen  = len(data)-1
     Energy     = np.zeros([max_eigen],dtype=np.float64)
     print("max_eigen= ", max_eigen)
     #[s] count not empty elements
     cnt = 0
     for i in range(0,len(data)):
         if data[i]: # if data[i] is not empty
             tmp        = data[i].split()
             Energy[cnt] = float(tmp[1])
             cnt       += 1

f        = open("FullDiag.dat", 'wt')
e_temp   = -4
int_temp =  1
temp     = -1
while temp < 100:
    if int_temp == 19:
        int_temp  = 1
        e_temp   += 1
    temp     = (int_temp/2.0+0.5)*pow(10,e_temp)
    print(int_temp,e_temp,temp)
    beta     = 1.0/temp
    int_temp += 1
    Z        = 0.0
    all_E    = 0.0
    all_E2   = 0.0
    for cnt in range(0,max_eigen):
        norm_ene  = Energy[cnt]-Energy[0]
        all_E    += Energy[cnt]*math.exp(-beta*norm_ene)
        all_E2   += (Energy[cnt]**2)*math.exp(-beta*norm_ene)
        Z        += math.exp(-beta*norm_ene)
    tmp_e = all_E/Z
    tmp_C = beta*beta*(all_E2/Z-(all_E/Z)**2)
    f.write(" {0:16f} ".format(temp) \
    +" {0:.16f}   ".format(tmp_e)         \
    +" {0:.16f}   ".format(tmp_C)           \
    +"\n")
    #print(beta,tmp)
f.close()
