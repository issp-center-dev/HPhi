import numpy as np
import os
import copy
import math
import cmath
import sys

def mainBasic(num_sample,dir_Norm):
    #num_sample = int(sys.argv[1])
    #dir_Norm = "dir_Norm_2"
    #tpq_type = sys.argv[3]
    tpq_type = "cTPQ" # fixed for cTPQ
    Norm,InvTemp,Ene,Ene2,Spc,Sz,Sz2,chi = read_file(num_sample)
    num_step    =  Norm.shape[1]

    ave_InvTemp = np.mean(InvTemp,axis=0)
    #err_InvTemp = np.std(InvTemp,axis=0,ddof=1)
    ave_Ene    = np.mean(Ene,axis=0)
    err_Ene    = np.std(Ene,axis=0,ddof=1)
    ave_Spc    = np.mean(Spc,axis=0)
    err_Spc    = np.std(Spc,axis=0,ddof=1)
    ave_Sz     = np.mean(Sz,axis=0)
    err_Sz     = np.std(Sz,axis=0,ddof=1)
    ave_chi    = np.mean(chi,axis=0)
    err_chi    = np.std(chi,axis=0,ddof=1)

    with open("%s/SimpleAverage.dat" % (dir_Norm),'w') as f:
         for k in range(0,num_step):
            if k==0:
                print("# inf %12.f %.12f %.12f %.12f %.12f %.12f %.12f %.12f " % (ave_Ene[k],err_Ene[k],ave_Spc[k],err_Spc[k],ave_Sz[k],err_Sz[k],ave_chi[k],err_chi[k]),file=f)
            else:
                temp = 1.0/ave_InvTemp[k]
                print(" %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f " % (temp,ave_Ene[k],err_Ene[k],ave_Spc[k],err_Spc[k],ave_Sz[k],err_Sz[k],ave_chi[k],err_chi[k]),file=f)

    phys_Z,phys_Ene,phys_Ene2,phys_Sz,phys_Sz2,phys_InvTemp  = CalcBasic(Norm,Ene,Ene2,Sz,Sz2,InvTemp,tpq_type)

    log_Z  = np.zeros([num_step],dtype=np.float64)
    with open("%s/Norm.dat" % (dir_Norm),'w') as f:
         tmp_log_Z =  -math.log(Norm[0][0])
         for k in range(0,num_step):
             tmp_log_Z += math.log(Norm[0][k])
             log_Z[k]   = tmp_log_Z
             print("%d %.16f %.16f " % (k,log_Z[k],InvTemp[0][k]),file=f)
 
    for cnt in range(num_sample):
        with open("%s/TPQ_%d.dat" % (dir_Norm,cnt),'w') as f2:
            for k in range(0,num_step):
                 print("%d %.12f %.12f %.12f %.12f %.12f %.12f" % (k,phys_InvTemp[cnt][k],phys_Z[cnt][k],phys_Ene[cnt][k],phys_Ene2[cnt][k],phys_Sz[cnt][k],phys_Sz2[cnt][k]),file=f2)
   
    return log_Z,phys_Z,phys_Ene,phys_Ene2,phys_Sz,phys_Sz2,phys_InvTemp

def mainPhys(num_sample,dir_Norm,dir_Phys,header):
    IPL_InvTemp,IPL_Z =  read_Norm(num_sample,dir_Norm)
    Phys              =  read_phys(num_sample,dir_Phys,header)
    tpq_type          =  "cTPQ"
    num_step          =  Phys.shape[1]
    num_phys          =  Phys.shape[2]

    print("num_phys",num_phys)
    Norm_Phys         = CalcPhys(IPL_Z,IPL_InvTemp,Phys,tpq_type)
    #Note: Norm_Phys[0] = cnt, Norm_Phys[1] = inve_temp, Norm_Phys[2-] physical quantities
    for cnt in range(num_sample):
        print(cnt)
        with open("%s/Norm_%s_set%d.dat" % (dir_Phys,header,cnt),'w') as f2:
            for k in range(0,num_step):
                true_cnt  = int(Norm_Phys[cnt][k][0])
                print("%d     "% (true_cnt),end="",file=f2)
                print("%.12f  "% (IPL_Z[cnt][true_cnt]),end="",file=f2)
                for cnt_phys in range(1,num_phys):
                    print(" %12.12f " % (Norm_Phys[cnt][k][cnt_phys]),end="",file=f2)
                print(" ",file=f2)
    return Norm_Phys


def read_file(num_sample):
    in_file = "output/Norm_rand0.dat"

    with open("%s" % (in_file)) as f:
        data      = f.read()
        data      = data.split("\n")
    num_step    = int(len(data)-2)
    Norm        = np.zeros([num_sample,num_step],dtype=np.float64)
    InvTemp     = np.zeros([num_sample,num_step],dtype=np.float64)
    Ene         = np.zeros([num_sample,num_step],dtype=np.float64)
    Ene2        = np.zeros([num_sample,num_step],dtype=np.float64)
    Spc         = np.zeros([num_sample,num_step],dtype=np.float64)
    Sz          = np.zeros([num_sample,num_step],dtype=np.float64)
    Sz2         = np.zeros([num_sample,num_step],dtype=np.float64)
    chi         = np.zeros([num_sample,num_step],dtype=np.float64)

    for cnt_samp in range(num_sample):
        in_file = "output/Norm_rand%d.dat" % (cnt_samp)
        print(in_file)

        with open("%s" % (in_file)) as f:
            data      = f.read()
            data      = data.split("\n")
        for i in range(1,num_step+1):
            cnt       = i-1
            tmp       = data[i].split()
            Norm[cnt_samp][cnt] = float(tmp[1])**2
            #print(cnt,Norm[cnt])

        in_file = "output/SS_rand%d.dat" %(cnt_samp)
        print(in_file)

        with open("%s" % (in_file)) as f:
            data      = f.read()
            data      = data.split("\n")
            #print(len(data))
        for i in range(1,num_step+1):
            cnt        = i-1
            tmp        = data[i].split()
            InvTemp[cnt_samp][cnt]  = float(tmp[0])
            Ene[cnt_samp][cnt]   = float(tmp[1])
            Ene2[cnt_samp][cnt]  = float(tmp[2])
            Spc[cnt_samp][cnt]   = float(tmp[0])*float(tmp[0])*(float(tmp[2])-float(tmp[1])**2)

        in_file = "output/Flct_rand%d.dat" %(cnt_samp)
        print(in_file)

        with open("%s" % (in_file)) as f:
            data      = f.read()
            data      = data.split("\n")
            #print(len(data))
        for i in range(1,num_step+1):
            cnt        = i-1
            tmp        = data[i].split()
            Sz[cnt_samp][cnt]   = float(tmp[5])
            Sz2[cnt_samp][cnt]   = float(tmp[6])
            chi[cnt_samp][cnt]  = float(tmp[0])*(float(tmp[6])-float(tmp[5])**2)

    return  Norm,InvTemp,Ene,Ene2,Spc,Sz,Sz2,chi

def CalcBasic(Norm,Ene,Ene2,Sz,Sz2,InvTemp,tpq_type):
    num_sample  = Norm.shape[0]
    num_step    = Norm.shape[1]
    print(num_step)
    phys_InvTemp   = np.zeros([num_sample,num_step],dtype=np.float64)
    phys_Z      = np.zeros([num_sample,num_step],dtype=np.float64)
    phys_Ene    = np.zeros([num_sample,num_step],dtype=np.float64)
    phys_Ene2   = np.zeros([num_sample,num_step],dtype=np.float64)
    phys_Sz     = np.zeros([num_sample,num_step],dtype=np.float64)
    phys_Sz2    = np.zeros([num_sample,num_step],dtype=np.float64)
    #NB: cnt_samp=0, Z is always 1
    for k in range(0,num_step):
        phys_Z[0][k]       = 1.0       
        phys_Ene[0][k]     = Ene[0][k]      
        phys_Ene2[0][k]    = Ene2[0][k]
        phys_Sz[0][k]      = Sz[0][k]      
        phys_Sz2[0][k]     = Sz2[0][k]
        phys_InvTemp[0][k] = InvTemp[0][k]
 
    for cnt_samp in range(1,num_sample):
        #[s] k=0 for each sample, the 1st norm should be 1
        k = 0
        tot_Z      = 1
        tot_Ene    = Ene[cnt_samp][k]
        tot_Ene2   = Ene2[cnt_samp][k]
        tot_Sz     = Sz[cnt_samp][k]
        tot_Sz2    = Sz2[cnt_samp][k]
        phys_Z[cnt_samp][k]    = tot_Z        
        phys_Ene[cnt_samp][k]  = tot_Ene        
        phys_Ene2[cnt_samp][k] = tot_Ene2 
        phys_Sz[cnt_samp][k]  = tot_Sz
        phys_Sz2[cnt_samp][k] = tot_Sz2
        phys_InvTemp[cnt_samp][k] =  InvTemp[cnt_samp][k]
        #[e] k=0
 
        for k in range(1,num_step):
            if tpq_type == "mTPQ": #  this part will be modified, now not used
                if InvTemp[0][k] > InvTemp[cnt_samp][k]: 
                    ext_k  = k+1
                    if InvTemp[0][k] > InvTemp[cnt_samp][k+1]:
                        ext_k  = k+2
                        #print("A fatal error")
                elif InvTemp[0][k] < InvTemp[cnt_samp][k]:
                    ext_k  = k-1
                    if InvTemp[0][k] < InvTemp[cnt_samp][k-1]:
                        ext_k  = k-2
                        #print("B fatal error k=%d cnt_sampl=%d: %f %f %f" % (k,cnt_samp,InvTemp[0][k],InvTemp[cnt_samp][k],InvTemp[cnt_samp][k-1]))
                else:
                    print("fatal")
                ratio_beta = (InvTemp[0][k]-InvTemp[cnt_samp][k])/(InvTemp[cnt_samp][ext_k]-InvTemp[cnt_samp][k])
                IPL_Z      = (Norm[cnt_samp][ext_k]-Norm[cnt_samp][k])*ratio_beta+Norm[cnt_samp][k]
                IPL_Ene    = (Ene[cnt_samp][ext_k]-Ene[cnt_samp][k])*ratio_beta+Ene[cnt_samp][k]
                IPL_Ene2   = (Ene2[cnt_samp][ext_k]-Ene2[cnt_samp][k])*ratio_beta+Ene2[cnt_samp][k]
                IPL_Sz     = (Sz[cnt_samp][ext_k]-Sz[cnt_samp][k])*ratio_beta+Sz[cnt_samp][k]
                IPL_Sz2    = (Sz2[cnt_samp][ext_k]-Sz2[cnt_samp][k])*ratio_beta+Sz2[cnt_samp][k]
                IPL_InvTemp   = InvTemp[0][k]
            elif tpq_type == "cTPQ":
                IPL_Z       =  Norm[cnt_samp][k]
                IPL_Ene     =  Ene[cnt_samp][k]
                IPL_Ene2    =  Ene2[cnt_samp][k]
                IPL_Sz      =  Sz[cnt_samp][k]
                IPL_Sz2     =  Sz2[cnt_samp][k]
                IPL_InvTemp =  InvTemp[cnt_samp][k]
            tot_Z      = tot_Z*IPL_Z/Norm[0][k] # normalized by Z @ samp=0
            tot_Ene    = IPL_Ene*tot_Z
            tot_Ene2   = IPL_Ene2*tot_Z
            tot_Sz     = IPL_Sz*tot_Z
            tot_Sz2    = IPL_Sz2*tot_Z
            phys_Z[cnt_samp][k]    = tot_Z        
            phys_Ene[cnt_samp][k]  = tot_Ene        
            phys_Ene2[cnt_samp][k] = tot_Ene2 
            phys_Sz[cnt_samp][k]  = tot_Sz      
            phys_Sz2[cnt_samp][k] = tot_Sz2
            phys_InvTemp[cnt_samp][k] = IPL_InvTemp
              
    return phys_Z,phys_Ene,phys_Ene2,phys_Sz,phys_Sz2,phys_InvTemp
    
def CalcPhys(IPL_Z,IPL_InvTemp,Phys,tpq_type):
    num_sample  = Phys.shape[0]
    num_step    = Phys.shape[1]
    num_phys    = Phys.shape[2]
    print(num_sample,num_step,num_phys)
    Norm_Phys   = np.zeros([num_sample,num_step,num_phys],dtype=np.float64)
    for k in range(0,num_step):
        for cnt_phys in range(0,num_phys):
            Norm_Phys[0][k][cnt_phys]  = Phys[0][k][cnt_phys]
 
    for cnt_samp in range(1,num_sample):
        for k in range(1,num_step):
            true_cnt                  = int(Phys[cnt_samp][k][0])
            Norm_Phys[cnt_samp][k][0] = true_cnt
            Norm_Phys[cnt_samp][k][1] = IPL_InvTemp[cnt_samp][true_cnt]
            temp      = Phys[cnt_samp][k][0]
            if tpq_type=="mTPQ":
                IPL_temp = 1.0/IPL_InvTemp[cnt_samp][true_cnt]
                if IPL_temp > temp: 
                    ext_k  = k-1
                elif IPL_temp < temp:
                    ext_k  = k+1
                else:
                    print("fatal")
                ratio_T   = (IPL_temp-Phys[cnt_samp][k][1])/(Phys[cnt_samp][ext_k][1]-Phys[cnt_samp][k][1])
                for cnt_phys in range(2,num_phys):
                    #Norm_Phys[cnt_samp][k][cnt_phys]=IPL_Z[cnt_samp][true_cnt]*Phys[cnt_samp][k][cnt_phys]
                    DPhys = Phys[cnt_samp][ext_k][cnt_phys]-Phys[cnt_samp][k][cnt_phys]
                    Norm_Phys[cnt_samp][k][cnt_phys]=IPL_Z[cnt_samp][true_cnt]*(DPhys*ratio_T+Phys[cnt_samp][k][cnt_phys])
            elif tpq_type=="cTPQ":
                for cnt_phys in range(2,num_phys):
                    Norm_Phys[cnt_samp][k][cnt_phys]=IPL_Z[cnt_samp][true_cnt]*Phys[cnt_samp][k][cnt_phys]
    return Norm_Phys

def read_phys(num_sample,dir_Phys,header):
    in_file = "%s/%s_set0.dat"  % (dir_Phys,header)
    print(in_file)

    with open("%s" % (in_file)) as f:
        data      = f.read()
        data      = data.split("\n")
    num_step    = int(len(data)-1)
    tmp         = data[1].split()
    num_phys    = len(tmp)
    print(num_step)
    print(num_phys)
    Phys        = np.zeros([num_sample,num_step,num_phys],dtype=np.float64)

    for cnt_samp in range(num_sample):
        in_file = "%s/%s_set%d.dat" % (dir_Phys,header,cnt_samp)
        print(in_file)
        with open("%s" % (in_file)) as f:
            data      = f.read()
            data      = data.split("\n")
            #print(len(data))
        for i in range(0,num_step):
            cnt        = i
            tmp        = data[i].split()
            for cnt_phys in range(len(tmp)):
                Phys[cnt_samp][cnt][cnt_phys]  = float(tmp[cnt_phys])
    return  Phys

def read_Norm(num_sample,dir_Norm):
    in_file = "%s/TPQ_0.dat"%(dir_Norm)
    with open("%s" % (in_file)) as f:
        data      = f.read()
        data      = data.split("\n")
    num_step    = int(len(data)-2)
    print("num_step=",num_step)
    IPL_Z       = np.zeros([num_sample,num_step],dtype=np.float64)
    IPL_InvTemp = np.zeros([num_sample,num_step],dtype=np.float64)

    for cnt_samp in range(num_sample):
        in_file = "%s/TPQ_%d.dat" % (dir_Norm,cnt_samp)
        print(in_file)
        with open("%s" % (in_file)) as f:
            data      = f.read()
            data      = data.split("\n")
        for i in range(1,num_step+1):
            cnt       = i-1
            tmp       = data[i].split()
            IPL_InvTemp[cnt_samp][cnt] = float(tmp[1])
            IPL_Z[cnt_samp][cnt]       = float(tmp[2])

    return  IPL_InvTemp,IPL_Z



if __name__ == "__main__":
    main()
