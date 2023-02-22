import numpy as np
import math
import cmath

def func_site(all_i,list_org):
    Lx      = list_org[0]
    Ly      = list_org[1]
    Lz      = list_org[2]
    orb_num = list_org[3]

    orb  = all_i%orb_num
    site = (all_i-orb)/orb_num
    x    = site%Lx
    tmp  = (site-x)/Lx
    y    = tmp%Ly
    z    = (tmp-y)/Ly 
    list_site = [int(x),int(y),int(z),int(orb)]
    return list_site

def func_input(list_lat):
    with open('input.txt') as f:
        data = f.read()
    data = data.split("\n")
    #print('check input',data)
    #print(len(data))
    dict_lat = {}
  
    for i in data:
        tmp = i.split()
        if len(tmp)>0:
            for i in list_lat:
                if i == tmp[0]:
                    dict_lat[i] = tmp[1]
    return dict_lat

#[s] read input.txt
list_lat =['Lx','Ly','Lz','orb_num'] # list for lattice parametes
dict_lat = func_input(list_lat)     # read input.txt
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
print('orb_num = ',orb_num)

All_N = Lx*Ly*Lz*orb_num

#[s] initialize
list_org   = [Lx,Ly,Lz,orb_num]
#list_trans = [0,0,0,0]  # x,y,z,orb
#list_site  = [0,0,0,0]  # x,y,z,orb
#[e] initialize

#[s] allocate
Sq     = np.zeros([Lx+1,Ly+1],dtype=np.float64)
Sq_err = np.zeros([Lx+1,Ly+1],dtype=np.float64)
#[e] allocate
max_num = 1
for num_cg in range(0,max_num):
    file_name = "output/zvo_cisajscktalt_eigen"+"{0:1d}".format(num_cg)+".dat"

    with open(file_name) as f:
        data      = f.read()
        data      = data.split("\n")
        print(len(data))
    #[s] count not empty elements
    cnt = 0
    for i in range(0,len(data)):
        if data[i]: # if data[i] is not empty
            cnt += 1
    #print(cnt)
    cnt_max = cnt
    #[e] count not empty elements
   
    with open("Sq_eigen%d.dat" % (num_cg), 'w') as f:
        for kx in range(0,Lx+1):
            for ky in range(0,Ly+1):
                tmp_Sq  = 0.0
                tmp_Sz  = 0.0
                for cnt in range(0,cnt_max):
                    tmp             = data[cnt].split()  
                    all_i           = int(tmp[0]) 
                    list_site       = func_site(all_i,list_org)
                    i_x             = list_site[0]
                    i_y             = list_site[1]
                    #
                    all_j           = int(tmp[4]) 
                    list_site       = func_site(all_j,list_org)
                    j_x             = list_site[0]
                    j_y             = list_site[1]
    
                    theta           = 2*math.pi*kx*(i_x-j_x)/Lx+2*math.pi*ky*(i_y-j_y)/Ly
                    if cnt%6   == 0 or cnt%6 == 3:
                       tmp_Sq += 0.25*float(tmp[8])*math.cos(theta)
                       tmp_Sz += 0.25*float(tmp[8])*math.cos(theta)
                    elif cnt%6 == 1 or cnt%6 == 2:
                       tmp_Sq += -0.25*float(tmp[8])*math.cos(theta)
                       tmp_Sz += -0.25*float(tmp[8])*math.cos(theta)
                    else:
                       tmp_Sq += 0.5*float(tmp[8])*math.cos(theta)
                tmp_Sq  = tmp_Sq/(Lx*Ly)
                tmp_Sz  = tmp_Sz/(Lx*Ly)
                print(kx,ky,tmp_Sq,tmp_Sz)
                print("%d %d %12.8f %12.8f  " % (kx,ky,tmp_Sq,tmp_Sz), file=f)
            print(" " , file=f)
            print(" ")
