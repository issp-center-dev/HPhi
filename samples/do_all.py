import os
import sys
import subprocess
import glob
import time

#Set path(command) to Hphi
if len(sys.argv) == 2:
    path_to_HPhi = os.path.abspath(sys.argv[1])
else:
    print("Error")
    print("Usage: python do_all.py path_to_HPhi.")
    print("path_to_HPhi: relative or absolute path to HPhi.")
    exit(1)

#get all dictioary of tutorials
dir_all = sorted(glob.glob("tutorial*"))
current_dir = os.getcwd()
for idx, dir_name in enumerate(dir_all):
    #if dir_name != "tutorial_4.3":
    #    continue
    print("#########Test {}/{}#############".format(idx+1, len(dir_all)))
    print("Check: {}".format(dir_name))
    dir_name = os.path.join(current_dir, dir_name) 
    os.chdir(dir_name)
    #Run script
    if os.path.exists("run.sh"):
        print("Start calculation.")
        start = time.time()
        with open("run.sh", "r") as fr:
            lines = fr.readlines()
            for line in lines[1:]:
                line = line.replace("$1", path_to_HPhi)
                cmd = line.split()
                with open("std.out", "w") as fw:
                    subprocess.call(cmd, stdout=fw)
        with open("check_file.txt", "r") as fr:
            lines = fr.readlines()
            bool_test = True
            for line in lines:
                path_to_file = os.path.join(dir_name, "output", line.strip())
                if path_to_file is False:
                    bool_test = False
                    print("Error: {} does not exist".format(line.strip()))
        elapsed_time = time.time() - start
        print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")
        if bool_test is True:
            print("  OK  : {}".format(dir_name))
    else:
        print("Error : run.sh does not exist in {}".format(dir_name))
