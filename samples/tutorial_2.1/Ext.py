import numpy as np
import os
import copy
import math
import cmath
import sys

in_file = sys.argv[1]
print(in_file)

with open("%s" % (in_file)) as f:
    data      = f.read()
    data      = data.split("\n")
    print(len(data))
    i  = 1
    di = 1 
    with open("Ext_%s" % (in_file), 'w') as f:
        while i < len(data):
            i+=di
            if i < 20:
                di = 1
            elif i < 50:
	            di = 5
            elif i < 100:
	            di = 10
            elif i < 200:
	            di = 20
            elif i < 400:
	            di = 40
            elif i < 800:
	            di = 80
            elif i < 1600:
	            di = 160
            elif i < 3200:
	            di = 320
            elif i < 6400:
	            di = 640
            elif i < 12800:
	            di = 1280
            if i < len(data):
                 print(" %s  " % (data[i]), file=f)

    #for cnt in range(0,max_eigen-1):
        #print(" ", file=f)
