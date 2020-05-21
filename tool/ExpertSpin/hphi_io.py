from __future__ import print_function
import numpy as np
import math
import cmath
import itertools

def func_io_all(file_name,max_site,param):
    #[s] interaction 
    print("  output interaction = ",file_name)
    #[e] interaction 
    num_param = len(np.nonzero(param)[0])
    f        = open(file_name, 'wt')
    f.write("==================="+"\n")
    f.write("num "+"{0:8d}".format(num_param)+"\n")
    f.write("==================="+"\n")
    f.write("==================="+"\n")
    f.write("==================="+"\n")
    cnt = 0
    for all_i in range(0,max_site):
        for all_j in range(0,max_site):
            for spn_0,spn_1,spn_2,spn_3 in itertools.product([0,1],repeat=4):
                #print(spn_0,spn_1,spn_2,spn_3)
                tmp_param  = param[all_i][spn_0][all_i][spn_1][all_j][spn_2][all_j][spn_3]
                if abs(tmp_param) != 0:
                    cnt       += 1
                    #print(all_i,spn_0,all_i,spn_1,all_j,spn_2,all_j,spn_3,tmp_param.real,tmp_param.imag)
                    f.write(" {0:8d} ".format(all_i) \
                    +" {0:8d}   ".format(spn_0)     \
                    +" {0:8d}   ".format(all_i)     \
                    +" {0:8d}   ".format(spn_1)     \
                    +" {0:8d}   ".format(all_j)     \
                    +" {0:8d}   ".format(spn_2)     \
                    +" {0:8d}   ".format(all_j)     \
                    +" {0:8d}   ".format(spn_3)     \
                    +" {0:8f}   ".format(tmp_param.real) \
                    +" {0:8f}   ".format(tmp_param.imag) \
                    +"\n")
    f.close()

def func_io(file_name,param,out_type):
    #[s] interaction 
    print("  output interaction = ",file_name)
    #[e] interaction 
    num_param = len(np.nonzero(param)[0])
    f        = open(file_name, 'wt')
    f.write("==================="+"\n")
    f.write("num "+"{0:8d}".format(num_param)+"\n")
    f.write("==================="+"\n")
    f.write("==================="+"\n")
    f.write("==================="+"\n")
    cnt = 0
    for all_i in np.nonzero(param)[0]:
      all_j      = np.nonzero(param)[1][cnt]
      tmp_param  = param[all_i][all_j]
      cnt       += 1
      #print(all_i,all_j)
      if out_type == "two":
          f.write(" {0:8d} ".format(all_i) \
          +" {0:8d}   ".format(all_j)     \
          +" {0:8f}   ".format(tmp_param) \
          +"\n")
    f.close()

def func_mag(file_name,max_site,mag_h):
    print(file_name)
    num_param = 4*max_site
    f        = open(file_name, 'wt')
    f.write("==================="+"\n")
    f.write("num "+"{0:8d}".format(num_param)+"\n")
    f.write("==================="+"\n")
    f.write("==================="+"\n")
    f.write("==================="+"\n")
    for cnt in range(0,max_site):
      all_i     = cnt
      # 
      f.write(" {0:8d} ".format(all_i) \
      +" {0:8d}   ".format(0)          \
      +" {0:8d}   ".format(all_i)      \
      +" {0:8d}   ".format(1)          \
      +" {0:8f}   ".format(-0.5*mag_h)          \
      +" {0:8f}   ".format(0.5*mag_h)          \
      +"\n")
      # 
      f.write(" {0:8d} ".format(all_i) \
      +" {0:8d}   ".format(1)          \
      +" {0:8d}   ".format(all_i)      \
      +" {0:8d}   ".format(0)          \
      +" {0:8f}   ".format(-0.5*mag_h)          \
      +" {0:8f}   ".format(-0.5*mag_h)          \
      +"\n")
      #
      f.write(" {0:8d} ".format(all_i) \
      +" {0:8d}   ".format(0)          \
      +" {0:8d}   ".format(all_i)      \
      +" {0:8d}   ".format(0)          \
      +" {0:8f}   ".format(-0.5*mag_h)          \
      +" {0:8f}   ".format(0.0)          \
      +"\n")
      #
      f.write(" {0:8d} ".format(all_i) \
      +" {0:8d}   ".format(1)          \
      +" {0:8d}   ".format(all_i)      \
      +" {0:8d}   ".format(1)          \
      +" {0:8f}   ".format(0.5*mag_h)          \
      +" {0:8f}   ".format(0.0)          \
      +"\n")
      # 
    f.close()
