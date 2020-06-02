import numpy as np
import os
import copy
import math
import cmath
import read    #using read.py#

#[s] read input.txt
list_lat =['Lx','Ly','Lz','orb_num'] # list for lattice parametes
dict_lat = read.func_input(list_lat)     # read input.txt
print(list_lat)
print(dict_lat)
#[e] read input.txt

#[s] dict_lat -> parameters 
Lx       = int(dict_lat['Lx'])
Ly       = int(dict_lat['Ly'])
Lz       = int(dict_lat['Lz'])
orb_num  = int(dict_lat['orb_num'])
#[e] dict_lat -> parameters 
print('Lx = ',Lx)
print('Ly = ',Ly)
print('Lz = ',Lz)
print('orb_nun = ',orb_num)

All_N = Lx*Ly*Lz*orb_num

#[s] initialize
list_org   = [Lx,Ly,Lz,orb_num]
#[e] initialize
    
num_green  = 6*All_N*All_N
f        = open("greentwo.def", 'wt')
f.write("==================="+"\n")
f.write("loc "+"{0:8d}".format(num_green)+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
#[s] z and orb
list_trans = [0,0,0,0] # for z
#[e] z and orb
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
        f.write(" {0:8d} ".format(all_i)+" 0 " \
           +" {0:8d} ".format(all_i)+" 1 "     \
           +" {0:8d} ".format(all_j)+" 1 "     \
           +" {0:8d}   ".format(all_j)+" 0 "   \
           +"\n")
f.close()
