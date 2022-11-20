import numpy as np
import os
import copy
import sys
import math
import cmath
import random
import Norm_TPQ

def main():
    #[s] set param
    num_sample  = int(sys.argv[1]) 
    max_BS      = int(sys.argv[2]) 
    print(num_sample,max_BS)
    dir_Norm    = "dir_Norm"
    os.makedirs(dir_Norm, exist_ok=True)
    #[e] set param

    #[s] Calc Norm
    log_Z,phys_Z,phys_Ene,phys_Ene2,phys_Sz,phys_Sz2,phys_InvTemp = Norm_TPQ.mainBasic(num_sample,dir_Norm)
    #[e] Calc Norm
    #[s] Basic BS ene,spc,ent
    BS_Basic(max_BS,log_Z,phys_Z,phys_Ene,phys_Ene2,phys_Sz,phys_Sz2,phys_InvTemp)
    #[e] Basic BS ene,spc,ent
    #[s] Ohter phys. quantities
    dir_Phys       = "dir_Phys"
    os.makedirs(dir_Phys, exist_ok=True)
    #header         = "Sq"
    #BS_Phys(max_BS,phys_Z,dir_Norm,dir_Phys,header)
    #[e] Ohter phys. quantities
 

def BS_Basic(max_BS,log_Z,phys_Z,phys_Ene,phys_Ene2,phys_Sz,phys_Sz2,phys_InvTemp):
    num_sample  = phys_Z.shape[0]
    max_eigen   = len(log_Z) 
    max_set     = num_sample 
    BS_sample   = num_sample
    #[s] Boot strap [BS]
    BS_Z         = np.zeros([max_BS,max_eigen],dtype=np.float64)
    BS_Ene       = np.zeros([max_BS,max_eigen],dtype=np.float64)
    BS_Spc       = np.zeros([max_BS,max_eigen],dtype=np.float64)
    BS_Sz        = np.zeros([max_BS,max_eigen],dtype=np.float64)
    BS_Chi       = np.zeros([max_BS,max_eigen],dtype=np.float64)
    BS_Ent       = np.zeros([max_BS,max_eigen],dtype=np.float64)
    for set_BS in range(0,max_BS):
        print(set_BS)
        for cnt in range(0,max_eigen):
            tmp_z    = 0.0
            tmp_ene  = 0.0
            tmp_ene2 = 0.0
            tmp_Sz   = 0.0
            tmp_Sz2  = 0.0
            for cnt_BS in range(0,BS_sample):
                rand_set  = random.randrange(max_set)
                #print(set_BS,cnt_BS,rand_set)
                tmp_z    += phys_Z[rand_set][cnt]
                tmp_ene  += phys_Ene[rand_set][cnt]
                tmp_ene2 += phys_Ene2[rand_set][cnt]
                tmp_Sz   += phys_Sz[rand_set][cnt]
                tmp_Sz2  += phys_Sz2[rand_set][cnt]
            BS_Z[set_BS][cnt]    = tmp_z/BS_sample
            BS_Ene[set_BS][cnt]  = tmp_ene/tmp_z
            BS_Spc[set_BS][cnt]  = tmp_ene2/tmp_z-(tmp_ene/tmp_z)**2
            BS_Sz[set_BS][cnt]   = tmp_Sz/tmp_z
            BS_Chi[set_BS][cnt]  = tmp_Sz2/tmp_z-(tmp_Sz/tmp_z)**2
            BS_Ent[set_BS][cnt]  = math.log(tmp_z/BS_sample)+tmp_ene/tmp_z*phys_InvTemp[0][cnt]+log_Z[cnt]

    ave_BS_Z    = np.mean(BS_Z,axis=0)
    err_BS_Z    = np.std(BS_Z,axis=0,ddof=1)
    ave_BS_Ene  = np.mean(BS_Ene,axis=0)
    err_BS_Ene  = np.std(BS_Ene,axis=0,ddof=1)
    ave_BS_Spc  = np.mean(BS_Spc,axis=0)
    err_BS_Spc  = np.std(BS_Spc,axis=0,ddof=1)
    ave_BS_Sz   = np.mean(BS_Sz,axis=0)
    err_BS_Sz   = np.std(BS_Sz,axis=0,ddof=1)
    ave_BS_Chi  = np.mean(BS_Chi,axis=0)
    err_BS_Chi  = np.std(BS_Chi,axis=0,ddof=1)
    ave_BS_Ent  = np.mean(BS_Ent,axis=0)
    err_BS_Ent  = np.std(BS_Ent,axis=0,ddof=1)

    with open("BS_MaxBS%d.dat" % (max_BS), 'w') as f:
        for cnt in range(0,max_eigen):
            beta    = phys_InvTemp[0][cnt]
            ave_Ene = ave_BS_Ene[cnt]
            err_Ene = err_BS_Ene[cnt]
            ave_Spc = beta**2*(ave_BS_Spc[cnt])
            err_Spc = beta**2*err_BS_Spc[cnt]
            ave_Sz  = ave_BS_Sz[cnt]
            err_Sz  = err_BS_Sz[cnt]
            ave_Chi = beta*(ave_BS_Chi[cnt])
            err_Chi = beta*err_BS_Chi[cnt]
            ave_Ent = ave_BS_Ent[cnt]
            err_Ent = err_BS_Ent[cnt]
            #print(err_Ene,err_Spc)
            if cnt ==0:  
                print(" # inf ", end="", file=f)#1
            else:
                print(" %.16f  " % (1.0/phys_InvTemp[0][cnt]), end="", file=f)#1
            print(" %.16f  " % (ave_Ene), end="", file=f) #2
            print(" %.16f  " % (err_Ene), end="", file=f) #3
            print(" %.16f  " % (ave_Spc), end="", file=f) #4
            print(" %.16f  " % (err_Spc), end="", file=f) #5
            print(" %.16f  " % (ave_Ent), end="", file=f) #6
            print(" %.16f  " % (err_Ent), end="", file=f) #7
            print(" %.16f  " % (ave_Sz), end="", file=f) #8
            print(" %.16f  " % (err_Sz), end="", file=f) #9
            print(" %.16f  " % (ave_Chi), end="", file=f) #10
            print(" %.16f  " % (err_Chi), end="", file=f) #11
            print(" %.16f  " % (ave_BS_Z[cnt]), end="", file=f) #12
            print(" %.16f  " % (err_BS_Z[cnt]), end="", file=f) #13
            print(" %d     " % (cnt), end="", file=f) #14
            print(" ", file=f)
    #[e] Boot strap [BS]

def BS_Phys(max_BS,phys_Z,dir_Norm,dir_Phys,header):
    num_sample  = phys_Z.shape[0]
    max_set     = num_sample 
    BS_sample   = num_sample
    #[s] BS for ohter physical quantities
    dir_Phys       = "dir_Phys"
    header         = "Sq"
    Norm_Phys      = Norm_TPQ.mainPhys(num_sample,dir_Norm,dir_Phys,header)
    max_eigen_phys = Norm_Phys.shape[1]
    num_phys       = Norm_Phys.shape[2]
    print("num_phys=",num_phys,max_eigen_phys)
    #[s]BS
    BS_Phys      = np.zeros([max_BS,max_eigen_phys,num_phys],dtype=np.float64)
    for set_BS in range(0,max_BS):
        print("set_BS",set_BS)
        for cnt in range(0,max_eigen_phys):
            tmp_phys   = np.zeros([num_phys],dtype=np.float64)
            tmp_z      = 0.0
            for cnt_BS in range(0,BS_sample):
                rand_set  = random.randrange(max_set)
                true_cnt  = int(Norm_Phys[rand_set][cnt][0])
                tmp_z    += phys_Z[rand_set][true_cnt]
                for cnt_phys in range(0,num_phys):
                    tmp_phys[cnt_phys]  += Norm_Phys[rand_set][cnt][cnt_phys]
            #print(tmp_z)
            for cnt_phys in range(0,num_phys):
                BS_Phys[set_BS][cnt][cnt_phys]  = tmp_phys[cnt_phys]/tmp_z
    #[e]BS
    ave_BS_Phys = np.mean(BS_Phys,axis=0)
    err_BS_Phys = np.std(BS_Phys,axis=0,ddof=1)

    with open("%s_MaxBS%d.dat" % (header,max_BS), 'w') as f:
        for cnt in range(0,max_eigen_phys):
            if cnt == 0:
                print(" #  inf " , end="", file=f)#1
            else:
                temp    = 1.0/Norm_Phys[0][cnt][1]
                print(" %.16f  " % (temp), end="", file=f)#1
            for cnt_phys in range(2,num_phys):
                print(" %.16f  " % (ave_BS_Phys[cnt][cnt_phys]), end="", file=f) #2n+1
                print(" %.16f  " % (err_BS_Phys[cnt][cnt_phys]), end="", file=f) #2n+2 
            print(" ", file=f)

if __name__ == "__main__":
    main()
