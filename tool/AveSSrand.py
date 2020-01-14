import sys
import os.path
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    prog='AveSSrand.py',
    description='Tool for TPQ calculations.',
    epilog='end',
    add_help=True,
)

parser.add_argument('-n', action='store', dest='nsamples',
                    nargs='?', default=5, type=int, choices=None,
                    help=('Total number of samples.'),
                    metavar=None)

parser.add_argument('-o', action='store', dest='output',
                    nargs='?', default="ave_TPQ.dat", type=str, choices=None,
                    help=('Name of output file.'),
                    metavar=None)


args = parser.parse_args()
set_num=args.nsamples
np.zeros(set_num)

#Get sample size
str_ = "./output/SS_rand0.dat"
with open(str_, "r") as f:
    lines = f.readlines()
DataTmp = np.zeros( (len(lines)-2, set_num) )
DataEne = np.zeros( (len(lines)-2, set_num) )
DataC = np.zeros( (len(lines)-2, set_num) )

for num in range(set_num):
    str_ = "./output/SS_rand"+str(num)+".dat"
    if os.path.isfile(str_) is not True:
        print('The file {} does not exist.'.format(str_))
        sys.exit()
    with open(str_, "r") as f:
        lines = f.readlines()        
        for count, line in enumerate(lines[2:]):
            data = line.split()
            DataTmp[count][num]=1.0/float(data[0])
            DataEne[count][num]=float(data[1])
            DataC[count][num]=pow(float(data[0]),2)*( float(data[2])-pow(float(data[1]),2) )

with open(args.output, "w") as f:
    for count, datatmp in enumerate(DataTmp):
        f.write("{} {} {} {} {} {}\n".format(datatmp.mean(), datatmp.std(), \
                                           DataEne[count].mean(), DataEne[count].std(),\
                                           DataC[count].mean(), DataC[count].std()))
