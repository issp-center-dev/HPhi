#title           :testTE.py
#description     :This is a test tool for time evolution method
#author          :Kazuyoshi Yoshimi
#date            :20171221
#version         :1.0
#python_version  :2.7.14  
#==============================================================================

# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import shutil
import subprocess
import argparse
#import scipy.misc as scm

def StandardIni(L=4, t=1, U=4, nelec=4, TwoSz=0, J=1, Model = "\"Fermion Hubbard\""):
    Info = {}
    Info["L"] = L
    Info["model"] = Model
    Info["method"] = "\"Lanczos\""
    Info["lattice"] = "\"chain\""
    if Model == "\"Fermion Hubbard\"":
        Info["t"] = t
        Info["U"] = U
        Info["nelec"] = nelec
    elif Model == "\"Spin\"":
        Info["L"] = 2*L
        Info["J"] = J
    elif Model == "\"Kondo\"":
        Info["t"] = t
        Info["U"] = U
        Info["J"] = t
        Info["nelec"] = nelec/2
    Info["EigenVecIO"] = "\"Out\""
    Info["2Sz"] = TwoSz
    return Info
    
def MakeStandard(Info, fileName):
    """
    :type Info: dictionary
    :type fileName: string 
    :rtype: None
    """

    f=open(fileName, "w")
    for _info in Info:
        strInfo=str(_info)+" = "+ str(Info[_info]) + "\n"
        f.write(strInfo)
    f.close()
    
def MakeInterAll(Info, InteractionType):
    """
    :type Info: dictionary
    :type InteractionType: string 
    :rtype: InterAllArray: list
    """

    L=Info["L"]
    InterAllArray=[]
    if(Info["model"]=="\"Fermion Hubbard\""):
        if(InteractionType == "Normal"):        
            for x in range(0, L):
                lattice_origin=x
                lattice_forward=(x+1)%L
                InterAllArray.append([lattice_origin, 0, lattice_forward, 0, lattice_forward, 1, lattice_origin, 1 ])
        elif(InteractionType == "Diagonal"):
            for x in range(0, L):
                lattice_origin=x
                lattice_forward=(x+1)%L
                InterAllArray.append([lattice_origin, 0, lattice_forward, 1, lattice_forward, 1, lattice_origin, 0 ])
    elif(Info["model"]=="\"Spin\""):
            for x in range(0, L):
                lattice_origin=x
                lattice_forward=(x+1)%L
                InterAllArray.append([lattice_origin, 0, lattice_origin, 1, lattice_forward, 1, lattice_forward, 0 ])
    elif(Info["model"]=="\"Kondo\""):  
            for x in range(0, L):
                lattice_origin=x
                lattice_forward=(x+1)%L
                InterAllArray.append([lattice_origin, 0, lattice_origin, 1, lattice_forward, 1, lattice_forward, 0 ])
    else:
        return False        
    return InterAllArray

def PrintInterAll(InterAllArray, filename):
    """
    :type InterAllArray: list
    :type filename: string
    :rtype: None
    """
    str_filename=filename
    f = open(str_filename,"w")
    f.write("========================\n")
    f.write("NInterAll "+ str(2*len(InterAllArray)) + "\n")
    f.write("========================\n")
    f.write("========zInterAll=======\n")
    f.write("========================\n")
    for _list in InterAllArray:
        f.write(str(_list[0])+" "+str(_list[1])+ " " + \
                str(_list[2]) + " " + str(_list[3]) + " "  + \
                str(_list[4]) +" "+str(_list[5])+ " " + \
                str(_list[6]) + " " + str(_list[7])+ " 1.0 0.0 \n")
        f.write(str(_list[6])+" "+str(_list[7])+ " " + \
                str(_list[4]) + " " + str(_list[5]) + " "  + \
                str(_list[2]) +" "+str(_list[3])+ " " + \
                str(_list[0]) + " " + str(_list[1])+ " 1.0 0.0 \n")
    f.close()

def PrintTE(InterAllArray, InfoTE, filename):
    """
    :type InterAllArray: list
    :type InfoTE: dictionay
    :type filename: string
    :rtype: None
    """
    #InfoTE = [NTimeStep, dt]
    dt = InfoTE["dt"]
    str_filename=filename
    f = open(str_filename,"w")
    f.write("========================\n")
    f.write("AllTimeStep "+ str(InfoTE["NTimeSteps"]) + "\n")
    f.write("========================\n")
    f.write("========TwoBody Time Evolution=======\n")
    f.write("========================\n")
    t=0.0
    for _n in range(InfoTE["NTimeSteps"]):
        f.write(str(t)+ " " + str(2*len(InterAllArray))+ "\n")
        for _list in InterAllArray:
            f.write(str(_list[0])+" "+str(_list[1])+ " " + \
                    str(_list[2]) + " " + str(_list[3]) + " "  + \
                    str(_list[4]) +" "+str(_list[5])+ " " + \
                    str(_list[6]) + " " + str(_list[7])+ " 1.0 0.0 \n")
            f.write(str(_list[6])+" "+str(_list[7])+ " " + \
                    str(_list[4]) + " " + str(_list[5]) + " "  + \
                    str(_list[2]) +" "+str(_list[3])+ " " + \
                    str(_list[0]) + " " + str(_list[1])+ " 1.0 0.0 \n")
        t += dt
    f.close()


def replaceStr(str_file_Name, str_origin, str_after):
    """
    :type str_file_Name: string
    :type str_origin: string
    :type str_after: string
    :rtype: None
    """

    rf=None
    wf=None
    try:
        rf=open(str_file_Name, "r")
        wf=open("temp_file", "w")
        for line in rf:
            if line.find(str_origin) != -1:
                line = line.replace(str_origin, str_after)
            wf.write(line)
    finally:
        rf.close()
        wf.close()

    if os.path.isfile(str_file_Name) and os.path.isfile("temp_file"):
        os.remove(str_file_Name)
        os.rename("temp_file", str_file_Name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='testTE.py',
        description="""test tool for time evolution method. \n
        1. Calculate ground state by Lanczos method (namelist1.def). \n
        2. Calculate physical quantities with zInterAll.dat 
        by time evolution method (namelist2.def, output:SS2.dat) \n
        3. Calculate physical quantities with tetwo.def 
        by time evolution method (namelist3.def, output:SS3.dat)
        4. Check whethere values of SS2.dat mathes those of SS3.dat (by manual). \n        
        """,
        add_help=True,
    )

    parser.add_argument('-p', '--path', action='store', dest='path',
                        nargs='?', default='./HPhi', type=str, choices=None,
                        help=('Path to HPhi.'),
                        metavar=None)
    parser.add_argument('-m', '--model', action='store', dest='model',
                        nargs='?', default='Fermion Hubbard', type=str, choices=None,
                        help=('model: Fermion Hubbard, Spin, Kondo'),
                        metavar=None)
    parser.add_argument('-t', '--type', action='store', dest='interaction',
                        nargs='?', default='Normal', type=str, choices=None,
                        help=('Interaction type: Normal, Diagonal'),
                        metavar=None)
    parser.add_argument('-mpi', '--mpi', action='store', dest='mpi',
                        nargs='?', default="", type=str, choices=None,
                        help=('MPI command'),
                        metavar=None)
    
    args = parser.parse_args()
    PathToHPhi=args.path
    Model = "\""+args.model+"\""
    InteractionType = args.interaction
    MPIRUN=str(args.mpi)
    InfoTE={"NTimeSteps":100, "dt":0.01}

    # Make standard input file

    sin=StandardIni(Model=Model)
    MakeStandard(sin, "stan.in")
    subprocess.call(PathToHPhi+" -sdry ./stan.in", shell=True)

    #[s] Make CalcMod file
    #For Lanczos
    replaceStr("calcmod.def", "CalcType   0", "CalcType   0")
    shutil.copyfile("calcmod.def", "calcmod1.def")
    #For TE
    replaceStr("calcmod.def", "CalcType   0", "CalcType   4")
    replaceStr("calcmod.def", "InputEigenVec   0", "InputEigenVec   1")
    replaceStr("calcmod.def", "OutputEigenVec   1", "OutputEigenVec   0")
    shutil.copyfile("calcmod.def", "calcmod2.def")
    #[e] Make CalcMod file

    #[s] Make Namelist file
    #For Lanczos
    replaceStr("namelist.def", "CalcMod  calcmod.def", "CalcMod  calcmod1.def")
    shutil.copyfile("namelist.def", "namelist1.def")
    #For TE
    replaceStr("namelist.def", "CalcMod  calcmod1.def", "CalcMod  calcmod2.def")
    replaceStr("namelist.def", "ModPara  modpara1.def", "ModPara  modpara2.def")
    f=open("namelist.def", "a")
    f.write("       InterAll  zInterAll.def\n")
    f.write("       TETwoBody  tetwo2.def\n")
    f.close()
    shutil.copyfile("namelist.def", "namelist2.def")
    #[e] Make Namelist file

    #[s] Make ModPara file
    #For Lanczos
    shutil.copyfile("modpara.def", "modpara1.def")
    #For TE
    replaceStr("modpara.def", "Lanczos_max    2000", "Lanczos_max    100")
    f=open("modpara", "a")
    f.write("       ExpandCoef     10\n")
    f.close()
    shutil.copyfile("modpara.def", "modpara2.def")
    #[e] Make ModPara file    

    #[s] Modify InterAll
    InterAllArray = MakeInterAll(sin, InteractionType)
    #[e] 
    
    PrintInterAll(InterAllArray, "zInterAll.def")
    #Lanczos
    subprocess.call(MPIRUN+" "+PathToHPhi+" -e ./namelist1.def", shell=True)
    #TE1 (not set components of tetwo.def)
    PrintInterAll([],  "zInterAll.def")    
    PrintTE(InterAllArray, InfoTE, "tetwo2.def")    
    subprocess.call(MPIRUN+" "+PathToHPhi+" -e ./namelist2.def", shell=True)
