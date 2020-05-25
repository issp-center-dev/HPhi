# Sample program for calculating dynamical spin structure factors
# by using HPhi
# 
# Copyright (C) 2016 Takeo Kato (The Univeirsity of Tokyo)

# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 

# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#-------------------------------------------------------------

# [Note] We check that this python script works in HPhi ver.3.0.0.
# If this script does not work, please report it in HPhi-garally's issue.

import subprocess
import math
import sys
import os

#Set path(command) to Hphi
if len(sys.argv) == 2:
    path_to_HPhi = os.path.abspath(sys.argv[1])
else:
    print("Error")
    print("Usage: python do_all.py path_to_HPhi.")
    print("path_to_HPhi: relative or absolute path to HPhi.")
    exit(1)


L = 12
f = open('StdFace.def','w')
f.write('L = %i\n' % (L))
std_text = """model = "Spin"
method = "Lanczos"
lattice = "chain"
J = 1.0
2Sz = 0
2S = 1
"""
f.write(std_text)
f.close()
subprocess.call([path_to_HPhi,'-sdry','StdFace.def'])

f = open('calcmod.def','a')
f.write('OutputEigenVec 1\n')
f.close()

# first run (get eigenvalues, eigenvectors)
subprocess.call([path_to_HPhi,'-e','namelist.def'])
# get the ground-state energy
f = open('output/zvo_energy.dat','r')
energy = 0.0
for line in f:
    line_array = line.split('  ')
    if (line_array[0] == 'Energy') :
        energy = float(line_array[1])
f.close()

f = open('calcmod.def','a')
f.write('CalcSpec 1\n')
f.close()

f = open('modpara.def','a')
f.write('OmegaMin -10\n')
f.write('OmegaMax 0\n')
f.write('NOmega 100\n')
f.write('OmegaIm 0.1\n') 
f.close()

for i in range(1, L):
    f = open('pair.def','w')
    pair_text = """===============================
NCisAitCjtAjs      %i
===============================
====== PairExcitation =======
===============================
""" % (L*2)
    f.write(pair_text)
    for j in range(L):
        wr = math.cos(2.0*math.pi*float(i)*float(j)/float(L))
        wi = math.sin(2.0*math.pi*float(i)*float(j)/float(L))
        f.write('%i 0 %i 0 0 %f %f\n' % (j,j,wr,wi))
        f.write('%i 1 %i 1 0 %f %f\n' % (j,j,-wr,-wi))        
    f.close()
    print('Wavenumber %i\n' % (i))
    subprocess.call([path_to_HPhi,'-e','namelist.def'])
    subprocess.call(['cp','output/zvo_DynamicalGreen.dat','spectrum%i.dat' % (i)])
f = open('spectrum.dat','w')
for i in range(1, L):
    g = open('spectrum%i.dat' % (i % L),'r')
    for line in g:
        line_array = line.split(' ')
        f.write('%i %f %f\n' % (i,float(line_array[0])-energy,-float(line_array[3])))
    g.close()
    f.write('\n')
subprocess.call(["mv", "spectrum.dat", "spectrum.dat.bak"])
subprocess.call(["rm spectrum*.dat"], shell=True)
subprocess.call(["mv", "spectrum.dat.bak", "spectrum.dat"])
subprocess.call(["cp", "spectrum.dat", "output/spectrum.dat"])

f.close()
