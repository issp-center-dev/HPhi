import numpy as np
import os
import sys
import copy
import math
import cmath

def main():
    #[s] read modpara
    list_mod      =['Nsite','Lanczos_max','NumAve','ExpecInterval'] 
    dict_mod      = func_mod("modpara.def",list_mod)    
    max_set       = int(dict_mod['NumAve']) 
    max_eigen     = int(dict_mod['Lanczos_max'])
    #[e] read modpara

    #
    Ref_ave_Temp     = np.zeros([max_eigen],dtype=np.float64)
    Ref_err_Temp     = np.zeros([max_eigen],dtype=np.float64)
    Ref_ave_Ene      = np.zeros([max_eigen],dtype=np.float64)
    Ref_err_Ene      = np.zeros([max_eigen],dtype=np.float64)
    # 
    ave_Temp     = np.zeros([max_eigen],dtype=np.float64)
    err_Temp     = np.zeros([max_eigen],dtype=np.float64)
    # 
    InvTemp      = np.zeros([max_set,max_eigen],dtype=np.float64)
    ave_InvTemp  = np.zeros([max_eigen],dtype=np.float64)
    err_InvTemp  = np.zeros([max_eigen],dtype=np.float64)
    Ene          = np.zeros([max_set,max_eigen],dtype=np.float64)
    ave_Ene      = np.zeros([max_eigen],dtype=np.float64)
    err_Ene      = np.zeros([max_eigen],dtype=np.float64)
    Ene2         = np.zeros([max_set,max_eigen],dtype=np.float64)
    ave_Ene2     = np.zeros([max_eigen],dtype=np.float64)
    err_Ene2     = np.zeros([max_eigen],dtype=np.float64)
    #
    Spc          = np.zeros([max_set,max_eigen],dtype=np.float64)
    ave_Spc      = np.zeros([max_eigen],dtype=np.float64)
    err_Spc      = np.zeros([max_eigen],dtype=np.float64)
    Ent          = np.zeros([max_set,max_eigen-1],dtype=np.float64)
    ave_Ent      = np.zeros([max_eigen-1],dtype=np.float64)
    err_Ent      = np.zeros([max_eigen-1],dtype=np.float64)
    #

    for cnt_set in range(0,max_set):
        with open("output/SS_rand%d.dat" % (cnt_set)) as f:
            data      = f.read()
            data      = data.split("\n")
            #[s] count not empty elements
            for i in range(1,len(data)):
                tmp_i              = i-1
                tmp                = data[i].split()
                #print(tmp)
                if len(tmp)>1: # if data[i] is not empty
                    InvTemp[cnt_set][tmp_i]   = float(tmp[0])
                    Ene[cnt_set][tmp_i]       = float(tmp[1])
                    Ene2[cnt_set][tmp_i]      = float(tmp[2])
                    Spc[cnt_set][tmp_i]       = (float(tmp[2])-float(tmp[1])**2)*float(tmp[0])*float(tmp[0])

    ave_InvTemp = np.mean(InvTemp,axis=0)
    err_InvTemp = np.std(InvTemp,axis=0,ddof=1)
    ave_Ene  = np.mean(Ene,axis=0)
    err_Ene  = np.std(Ene,axis=0,ddof=1)
    ave_Spc  = np.mean(Spc,axis=0)
    err_Spc  = np.std(Spc,axis=0,ddof=1)

    for cnt_set in range(0,max_set):
        #with open("ave_tpq_%d.dat" % (cnt_set), 'w') as f:
        tmp_Ent = 0.0
        for i in range(0,max_eigen-1):
            tmp_Ent         +=  Spc[cnt_set][i]*(1.0-InvTemp[cnt_set][i]/InvTemp[cnt_set][i+1])
            Ent[cnt_set][i]  =  tmp_Ent

    ave_Ent  = np.mean(Ent,axis=0)
    err_Ent  = np.std(Ent,axis=0,ddof=1)

    with open("ave_tpq.dat", 'w') as f:
        for cnt in range(1,max_eigen-1):
            temp          =   1.0/ave_InvTemp[cnt]
            temp_err      =   err_InvTemp[cnt]/ave_InvTemp[cnt]**2 
            ave_Temp[cnt] =   temp      
            err_Temp[cnt] =   temp_err      
            print(" %.16f  " % (temp), end="", file=f)#1
            print(" %.16f  " % (temp_err), end="", file=f)#2
            print(" %.16f  " % (ave_Ene[cnt]), end="", file=f) #3
            print(" %.16f  " % (err_Ene[cnt]), end="", file=f) #4
            print(" %.16f  " % (ave_Spc[cnt]), end="", file=f) #5
            print(" %.16f  " % (err_Spc[cnt]), end="", file=f) #6
            print(" %.16f  " % (ave_Ent[cnt]), end="", file=f) #7
            print(" %.16f  " % (err_Ent[cnt]), end="", file=f) #8
            print(" ", file=f)

    with open("reference.dat") as f:
        data      = f.read()
        data      = data.split("\n")
        #[s] count not empty elements
        for i in range(0,len(data)):
            tmp_i              = i+1
            tmp                = data[i].split()
            if len(tmp)>1: # if data[i] is not empty
                Ref_ave_Temp[tmp_i]   = float(tmp[0])
                Ref_err_Temp[tmp_i]   = float(tmp[1])
                Ref_ave_Ene[tmp_i]    = float(tmp[2])
                Ref_err_Ene[tmp_i]    = float(tmp[3])

    result = 0
    for cnt in range(1,max_eigen-2):
        #print(cnt,Ref_ave_Temp[cnt],ave_Temp[cnt],Ref_err_Temp[cnt],Ref_ave_Ene[cnt],ave_Ene[cnt],Ref_err_Ene[cnt]);
        diff_temp = abs(Ref_ave_Temp[cnt]-ave_Temp[cnt])
        diff_ene  = abs(Ref_ave_Ene[cnt] -ave_Ene[cnt])
        if diff_temp > max(2*Ref_err_Temp[cnt],1e-8) :
            result = -1
            print("fatatl error in temp ")
            print(cnt,Ref_ave_Temp[cnt],ave_Temp[cnt],Ref_err_Temp[cnt],Ref_ave_Ene[cnt],ave_Ene[cnt],Ref_err_Ene[cnt]);
        if diff_ene  > max(2*Ref_err_Ene[cnt],1e-8) :
            result = -1
            print("fatatl error in ene ")
            print(cnt,Ref_ave_Temp[cnt],ave_Temp[cnt],Ref_err_Temp[cnt],Ref_ave_Ene[cnt],ave_Ene[cnt],Ref_err_Ene[cnt]);

    sys.exit(result)

def func_mod(in_file,list_param):
    with open(in_file) as f:
        data = f.read()
    data = data.split("\n")
    dict_param = {}
    for i in data:
        tmp = i.split()
        if len(tmp)>0:
            for i in list_param:
                if i == tmp[0]:
                    dict_param[i] = tmp[1]
    return dict_param

if __name__ == "__main__":
    main()
