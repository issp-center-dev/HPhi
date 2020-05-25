import os
import subprocess
import glob

#get all dictioary of tutorials
dir_all = glob.glob("tutorial*")
path_to_HPhi = "/Users/k-yoshimi/program/HPhi_box/build/src/HPhi"
current_dir = os.getcwd()
for dir_name in dir_all:
    print("Check: {}".format(dir_name))
    dir_name = os.path.join(current_dir, dir_name) 
    os.chdir(dir_name)
    #Run script
    if os.path.exists("run.sh"):
        print("Start calculation.")
        with open("run.sh", "r") as fr:
            lines = fr.readlines()
            for line in lines[1:]:
                line = line.replace("$1", path_to_HPhi)
                cmd = line.split()
                subprocess.call(cmd, stdout=subprocess.PIPE)
        print("Check file.")
        with open("check_file.txt", "r") as fr:
            lines = fr.readlines()
            bool_test = True
            for line in lines:
                path_to_file = os.path.join(dir_name, "output", line.strip())
                if path_to_file is False:
                    bool_test = False
                    print("Error: {} does not exist".format(line.strip()))
        if bool_test is True:
            print("Success to do {}".format(dir_name))
        print("Finish calculation")
    else:
        print("Error: run.sh does not exist in {}".format(dir_name))
