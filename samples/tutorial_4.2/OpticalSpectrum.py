import subprocess
import numpy as np
import os
import sys
import read
import lattice

class CalcOpt:
    def __init__(self,exct,NOmega,OmegaMin,OmegaMax,OmegaIm,sign,header, output_dir="./output"):
        self.exct       = exct
        self.NOmega     = NOmega
        self.OmegaMin   = OmegaMin
        self.OmegaMax   = OmegaMax
        self.OmegaIm    = OmegaIm
        self.sign       = sign
        self.header     = header
        self.output_dir = output_dir

    def Make_Opt_Input(self):
        #[s] words to be changed
        ex_mod = []
        ex_mod.append('OmegaOrg')
        ex_mod.append('OmegaMin')
        ex_mod.append('OmegaMax')
        ex_mod.append('NOmega')
        ex_mod.append('OmegaIm')
        #[e] words to be changed
        #[s] omega P & M
        OmegaMax_P = self.OmegaMax
        OmegaMin_P = self.OmegaMin
        OmegaMax_M = self.OmegaMin+(self.OmegaMax-self.OmegaMin)/NOmega
        OmegaMin_M = -1.0*self.OmegaMax+(self.OmegaMax-self.OmegaMin)/NOmega
        OmegaIm_P  = self.OmegaIm
        OmegaIm_M  = -1.0*self.OmegaIm
        #[e] omega P & M
        with open("modpara.def") as f:
            lines = f.readlines()
            with open("Smodpara.def", "w") as fex:
                for line in lines:
                    words = line.split()
                    if len(words) == 0:
                        continue
                    for tmp_word in ex_mod:
                        if words[0] == tmp_word:
                            print(tmp_word)
                            continue
                    fex.write(line)
                fex.write("OmegaOrg       {}\n".format(self.energy_list[self.exct]))
                fex.write("NOmega         {}\n".format(self.NOmega))
                if self.sign > 0:
                    fex.write("OmegaMin       {}\n".format(OmegaMin_P))
                    fex.write("OmegaMax       {}\n".format(OmegaMax_P))
                    fex.write("OmegaIm        {}\n".format(OmegaIm_P))
                else:
                    fex.write("OmegaMin       {}\n".format(OmegaMin_M))
                    fex.write("OmegaMax       {}\n".format(OmegaMax_M))
                    fex.write("OmegaIm        {}\n".format(OmegaIm_M))

        with open("calcmod.def") as f:
            lines = f.readlines()
            with open("Scalcmod.def", "w") as fex:
                for line in lines:
                    words = line.split()
                    if words[0] == "CalcSpec" or words[0] == "InputEigenVec":
                        continue
                    if words[0] == "OutputEigenVec" or words[0] == "ReStart":
                        continue
                    fex.write(line)
                fex.write("  ReStart         0\n")
                fex.write("  CalcSpec        3\n")
                fex.write("  InputEigenVec   1\n")
                fex.write("  OutputEigenVec  0\n")

        with open("namelist.def") as f:
            lines = f.readlines()
            with open("Snamelist.def", "w") as fex:
                for line in lines:
                    words = line.split()
                    if len(words) == 0:
                        continue
                    if words[0] == "CalcMod" or words[0] == "SpectrumVec":
                        continue
                    if words[0] == "ModPara" or words[0] == "PairExcitation":
                        continue
                    fex.write(line)
                fex.write("  ModPara           Smodpara.def\n")
                fex.write("  CalcMod           Scalcmod.def\n")
                fex.write("  PairExcitation    Current.def\n")
                fex.write("  SpectrumVec       {}_eigenvec_{}\n".format(self.header,self.exct))

    def get_energies(self):
        energy_list = []
        with open(os.path.join(output_dir, "{}_energy.dat".format(header))) as f:
            lines = f.readlines()
            for line in lines:
                words = line.split()
                if len(words) != 0 and words[0] == "Energy":
                    energy_list.append(float(words[1]))
        self.energy_list = energy_list
        self.ene_min = energy_list[0]
        self.ene_max = energy_list[len(energy_list) - 1]
        return energy_list


output_dir = "./output"
header     = "zvo"
exct       = 0
NOmega     = 2000
OmegaMin   = 0.0
OmegaMax   = 20.0
OmegaIm    = 0.1
sign       = float(sys.argv[1])
print('check sign',sign,type(sign))
#
calcopt=CalcOpt(exct,NOmega,OmegaMin,OmegaMax,OmegaIm,sign,header)
energy_list=calcopt.get_energies()
calcopt.Make_Opt_Input()
#
#[s] read input.txt
list_lat =['Lx','Ly','Lz','orb_num'] # list for lattice parameters
dict_lat = read.func_input(list_lat) # read input.txt
#[e] read input.txt
#[s] dict_lat -> parameters
Lx = int(dict_lat['Lx'])
Ly = int(dict_lat['Ly'])
Lz = int(dict_lat['Lz'])
orb_num  = int(dict_lat['orb_num'])
#[e] dict_lat -> parameters
print('Lx = ',Lx)
print('Ly = ',Ly)
print('Lz = ',Lz)
print('orb_num = ',orb_num)
All_N=Lx*Ly*Lz*orb_num
list_org = [Lx,Ly,Lz,orb_num]
list_trans = [0, 0, 0, 0]  # x,y,z,orb
with open("Current.def", 'w') as f:
    print("======optical conductivity ",  file=f)
    print("NCurrent 2                 ",  file=f)
    print("======optical conductivity ",  file=f)
    print("======optical conductivity ",  file=f)
    print("======optical conductivity ",  file=f)
    for leftright in range(2):
        print("%4d" % (4*All_N), file=f)
        for spin_i in range(0,2):
            for all_i in range(0,All_N):
                list_trans[0] = 1 # only +1 for x direction
                list_site = lattice.func_site(all_i, list_org)
                all_j     = lattice.func_strans(list_trans, list_site, list_org)
                print("%4d %4d %4d %4d 1 0 1 " % (all_i, spin_i, all_j, spin_i), file=f)
                print("%4d %4d %4d %4d 1 0 -1 " % (all_j, spin_i, all_i, spin_i), file=f)
