import numpy as np
import os
import copy
import math
import cmath
import read    #using read.py#
import hphi_io #using hphi_io.py#

#cnt_name   = "input"
#output_dir = "dir_"+"{}".format(cnt_name)
output_dir = "./"
#os.makedirs(output_dir,exist_ok=True)
tmp_sdt  = "pair.txt"
num      = read.func_count(tmp_sdt)
#[s] set input.
list_input = ['Ns','exct']
dict_input = read.func_input(list_input)
max_site = int(dict_input["Ns"])
exct     = 2
All_N    = max_site
#[e] set input.
#[s] interaction
Ising     = np.zeros([max_site,max_site],dtype=np.float64)
Exchange  = np.zeros([max_site,max_site],dtype=np.float64)
PairLift  = np.zeros([max_site,max_site],dtype=np.float64)
InterAll  = np.zeros([max_site,2,max_site,2,max_site,2,max_site,2],dtype=np.complex128)
#[e] interaction
#[s] read pair.txt
siteI   = np.zeros([num],dtype=np.int64)
siteJ   = np.zeros([num],dtype=np.int64)
intT1   = np.zeros([num],dtype=np.compat.unicode)
intT2   = np.zeros([num],dtype=np.compat.unicode)
para    = np.zeros([num],dtype=np.float64)
read.func_readpair(tmp_sdt,siteI,siteJ,intT1,intT2,para)
#[e] read pair.txt

for cnt in range(num):
    #[s] for diagonal part
    if intT1[cnt] == 'x' and intT2[cnt] == 'x':
        PairLift[siteI[cnt]][siteJ[cnt]] += 0.25*para[cnt]
        Exchange[siteI[cnt]][siteJ[cnt]] += 0.25*para[cnt]
    if intT1[cnt] == 'y' and intT2[cnt] == 'y':
        PairLift[siteI[cnt]][siteJ[cnt]] += -0.25*para[cnt]
        Exchange[siteI[cnt]][siteJ[cnt]] += 0.25*para[cnt]
    if intT1[cnt] == 'z' and intT2[cnt] == 'z':
        Ising[siteI[cnt]][siteJ[cnt]] += para[cnt]
    #[e] for diagonal part
    # S^{+}=P, S^{-}=M
    # Sx=(1/2)[  P(I) + M(I)]
    # Sy=(i/2)[ -P(I) + M(I)]
    # Sz=(1/2)[  U(I) - D(I)]
    #[s] for non-diagonal part
    if (intT1[cnt] == 'x' and intT2[cnt] == 'y') or (intT1[cnt] == 'y' and intT2[cnt] == 'x'):
        I = siteI[cnt]
        J = siteJ[cnt]
        # x(I)y(J)  =  (i/4) [P(I) + M(I)]*[-P(J) + M(J)]
        # x(J)y(I)  =  (i/4) [P(J) + M(J)]*[-P(I) + M(I)]      
        # x(I)y(J) + x(J)y(I) = (i/2)*[-P(I)P(J)+M(I)M(J)]
        InterAll[I][0][I][1][J][0][J][1] += complex(0,-0.5*para[cnt]) # -(i/2)*P(I)P(J)
        InterAll[J][1][J][0][I][1][I][0] += complex(0,0.5*para[cnt])  #  (i/2)*M(J)M(I)
    if (intT1[cnt] == 'x' and intT2[cnt] == 'z') or (intT1[cnt] == 'z' and intT2[cnt] == 'x'):
        I = siteI[cnt]
        J = siteJ[cnt]
        # x(I)z(J)   =  (1/4) [P(I) + M(I)] * [U(J)-D(J)]
        # x(J)z(I)   =  (1/4) [P(J) + M(J)] * [U(I)-D(I)] 
        InterAll[I][0][I][1][J][0][J][0] +=  0.25*para[cnt] #  (1/4) * P(I)U(J) 
        InterAll[J][0][J][0][I][1][I][0] +=  0.25*para[cnt] #  (1/4) * U(J)M(I) 

        InterAll[I][0][I][1][J][1][J][1] += -0.25*para[cnt] # -(1/4) * P(I)D(J) 
        InterAll[J][1][J][1][I][1][I][0] += -0.25*para[cnt] # -(1/4) * D(J)M(I) 

        InterAll[J][0][J][1][I][0][I][0] +=  0.25*para[cnt] #  (1/4) * P(J)U(I) 
        InterAll[I][0][I][0][J][1][J][0] +=  0.25*para[cnt] #  (1/4) * U(I)M(J) 

        InterAll[J][1][J][0][I][1][I][1] += -0.25*para[cnt] # -(1/4) * M(J)D(I) 
        InterAll[I][1][I][1][J][0][J][1] += -0.25*para[cnt] # -(1/4) * D(I)P(J) 
    if (intT1[cnt] == 'y' and intT2[cnt] == 'z') or (intT1[cnt] == 'z' and intT2[cnt] == 'y'):
        I = siteI[cnt]
        J = siteJ[cnt]
        # y(I)z(J)         =  (i/4) [-P(I) + M(I)] * [U(J)-D(J)]
        # y(J)z(I)         =  (i/4) [-P(J) + M(J)] * [U(I)-D(I)]
        InterAll[I][0][I][1][J][0][J][0] += complex(0,-0.25*para[cnt]) # -(i/4)* P(I)U(J)
        InterAll[J][0][J][0][I][1][I][0] += complex(0,0.25*para[cnt])  #  (i/4)* U(J)M(I)

        InterAll[I][0][I][1][J][1][J][1] += complex(0,0.25*para[cnt])  #  (i/4)* P(I)D(J)
        InterAll[J][1][J][1][I][1][I][0] += complex(0,-0.25*para[cnt]) # -(i/4)* D(J)M(I)

        InterAll[J][0][J][1][I][0][I][0] += complex(0,-0.25*para[cnt]) # -(i/4)* P(J)U(I)
        InterAll[I][0][I][0][J][1][J][0] += complex(0,0.25*para[cnt])  #  (i/4)* U(I)M(J)

        InterAll[J][1][J][0][I][1][I][1] += complex(0,-0.25*para[cnt]) # -(i/4)* M(J)D(I) 
        InterAll[I][1][I][1][J][0][J][1] += complex(0,0.25*para[cnt])  #  (i/4)* D(I)P(J)
    #if intT1[cnt] == 'y' and intT2[cnt] == 'x':
    #    print('should be xy')
    #if intT1[cnt] == 'z' and intT2[cnt] == 'x':
    #    print('should be xz')
    #if intT1[cnt] == 'z' and intT2[cnt] == 'y':
     #   print('should be yz')
    #[e] for non-diagonal part
#[s] io hamiltonians
hphi_io.func_io("{}".format(output_dir)+"Ising.def",Ising,"two")
hphi_io.func_io("{}".format(output_dir)+"Exchange.def",Exchange,"two")
hphi_io.func_io("{}".format(output_dir)+"PairLift.def",PairLift,"two")
hphi_io.func_io_all("{}".format(output_dir)+"InterAll.def",max_site,InterAll)
#[e] io hamiltonians
 
f        = open("{}".format(output_dir)+"calcmod_cg.def", 'wt')
f.write("  #CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, 4:Time-evolution"+"\n")
f.write("  #CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC"+"\n")
f.write("  #Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart"+"\n")
f.write("  #CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, 4:Restart, 5:Restart&Save"+"\n")
f.write("  CalcType   3"+"\n")
f.write("  CalcModel   4"+"\n")
f.write("  ReStart   0"+"\n")
f.write("  CalcSpec   0"+"\n")
f.write("  CalcEigenVec   0"+"\n")
f.write("  InitialVecType   0"+"\n")
f.write("  InputEigenVec   0"+"\n")
f.write("  OutputEigenVec   0"+"\n")
f.close()

f        = open("{}".format(output_dir)+"calcmod_tpq.def", 'wt')
f.write("  #CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, 4:Time-evolution"+"\n")
f.write("  #CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC"+"\n")
f.write("  #Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart"+"\n")
f.write("  #CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, 4:Restart, 5:Restart&Save"+"\n")
f.write("  CalcType   1"+"\n")
f.write("  CalcModel   4"+"\n")
f.write("  ReStart   0"+"\n")
f.write("  CalcSpec   0"+"\n")
f.write("  CalcEigenVec   0"+"\n")
f.write("  InitialVecType   0"+"\n")
f.write("  InputEigenVec   0"+"\n")
f.write("  OutputEigenVec   0"+"\n")
f.close()


f        = open("{}".format(output_dir)+"namelist_cg.def", 'wt')
f.write("  ModPara       modpara.def"+"\n")
f.write("  CalcMod       calcmod_cg.def"+"\n")
f.write("  LocSpin       locspn.def"+"\n")
f.write("  Ising         Ising.def"+"\n")
f.write("  Exchange      Exchange.def"+"\n")
f.write("  Pairlift      PairLift.def"+"\n")
f.write("  InterAll      InterAll.def"+"\n")
f.write("  OneBodyG      greenone.def"+"\n")
f.write("  TwoBodyG      greentwo.def"+"\n")
f.close()

f        = open("{}".format(output_dir)+"namelist_tpq.def", 'wt')
f.write("  ModPara       modpara.def"+"\n")
f.write("  CalcMod       calcmod_tpq.def"+"\n")
f.write("  LocSpin       locspn.def"+"\n")
f.write("  Ising         Ising.def"+"\n")
f.write("  Exchange      Exchange.def"+"\n")
f.write("  Pairlift      PairLift.def"+"\n")
f.write("  InterAll      InterAll.def"+"\n")
f.write("  TwoBodyG      greentwo.def"+"\n")
f.close()
 

f        = open("./{}/".format(output_dir)+"modpara.def", 'wt')
f.write("--------------------  "+"\n")
f.write("Model_Parameters   0  "+"\n")
f.write("--------------------  "+"\n")
f.write("HPhi_Cal_Parameters  "+"\n")
f.write("--------------------  "+"\n")
f.write("CDataFileHead  zvo  "+"\n")
f.write("CParaFileHead  zqp  "+"\n")
f.write("--------------------  "+"\n")
f.write("Nsite          {}".format(max_site)+"\n")
f.write("Ncond          {}".format(max_site)+"\n")
f.write("Lanczos_max    2000   "+"\n")
f.write("initial_iv     -1    "+"\n")
f.write("exct           {}".format(exct)+"\n")
f.write("LanczosEps     14   "+"\n")
f.write("LanczosTarget  2   "+"\n")
f.write("LargeValue     30  "+"\n")
f.write("NumAve         5   "+"\n")
f.write("ExpecInterval  20  "+"\n")
f.close()

num_loc  = max_site
f        = open("{}".format(output_dir)+"locspn.def", 'wt')
f.write("==================="+"\n")
f.write("loc "+"{0:8d}".format(num_loc)+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
for all_i in range(0,All_N):
     f.write(" {0:8d} ".format(all_i)+" 1 " \
     +"\n")
f.close()

num_green  = 0
f        = open("{}".format(output_dir)+"greenone.def", 'wt')
f.write("==================="+"\n")
f.write("loc "+"{0:8d}".format(num_green)+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
f.close()

num_green  = 6*All_N*All_N
f        = open("{}".format(output_dir)+"greentwo.def", 'wt')
f.write("==================="+"\n")
f.write("loc "+"{0:8d}".format(num_green)+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
for all_i in range(0,All_N):
    for all_j in range(0,All_N):
        f.write(" {0:8d} ".format(all_i)+" 0 " \
           +" {0:8d} ".format(all_i)+" 0 "     \
           +" {0:8d} ".format(all_j)+" 0 "     \
           +" {0:8d}   ".format(all_j)+" 0 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 0 " \
           +" {0:8d} ".format(all_i)+" 0 "     \
           +" {0:8d} ".format(all_j)+" 1 "     \
           +" {0:8d}   ".format(all_j)+" 1 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 1 " \
           +" {0:8d} ".format(all_i)+" 1 "     \
           +" {0:8d} ".format(all_j)+" 0 "     \
           +" {0:8d}   ".format(all_j)+" 0 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 1 " \
           +" {0:8d} ".format(all_i)+" 1 "     \
           +" {0:8d} ".format(all_j)+" 1 "     \
           +" {0:8d}   ".format(all_j)+" 1 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 0 " \
           +" {0:8d} ".format(all_i)+" 1 "     \
           +" {0:8d} ".format(all_j)+" 1 "     \
           +" {0:8d}   ".format(all_j)+" 0 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 1 " \
           +" {0:8d} ".format(all_i)+" 0 "     \
           +" {0:8d} ".format(all_j)+" 0 "     \
           +" {0:8d}   ".format(all_j)+" 1 "   \
           +"\n")
f.close()
