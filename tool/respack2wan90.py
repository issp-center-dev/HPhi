#!/usr/bin/python
#
# DCore -- Integrated DMFT software for correlated electrons
# Copyright (C) 2017 The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
from __future__ import print_function
import sys
import numpy


def res2wan(name_in, name_out):
    #
    print("  Convert from \"{0}\" to \"{1}\"".format(name_in, name_out))
    #
    # Input
    #
    with open(name_in, 'r') as f:
        for ii in range(3):
            line = f.readline()  # Skip
            print("    "+line, end="")
        line = f.readlines()
    #
    # Parse from RESPACK format
    #
    temp1 = [[]]
    nr = 0
    for iline in range(len(line)):
        if line[iline] == "\n" or line[iline] == " \n":
            temp1.append([])
            nr += 1
        else:
            temp1[nr].append(line[iline].split())
    #
    print("        Number of R : ", nr)
    norb = int(numpy.sqrt(len(temp1[0])-1)+0.1)
    print("    Number of bands : ", norb)

    irvec = numpy.zeros((nr, 3), numpy.int_)
    hopping = numpy.zeros((nr, norb, norb), numpy.complex_)

    for ir in range(nr):
        for ii in range(3):
            irvec[ir, ii] = int(temp1[ir][0][ii])

        ii = 0
        for iorb in range(norb):
            for jorb in range(norb):
                ii += 1
                hopping[ir, int(int(temp1[ir][ii][0]))-1, int(int(temp1[ir][ii][1]))-1] = \
                    float(temp1[ir][ii][2]) + 1.0j * float(temp1[ir][ii][3])
    #
    # Output to wannier90 format
    #
    with open(name_out, 'w') as f:

        print("Converted from RESPACK", file=f)
        print(norb, file=f)
        print(nr, file=f)
        for ir in range(nr):
            print("    1", end="", file=f)
            if ir % 15 == 14:
                print("", file=f)
        if nr % 15 != 0:
            print("", file=f)
        for ir in range(nr):
            for iorb in range(norb):
                for jorb in range(norb):
                    print("%5d%5d%5d%5d%5d%12.6f%12.6f" %
                          (irvec[ir, 0], irvec[ir, 1], irvec[ir, 2], jorb+1, iorb+1,
                           hopping[ir, jorb, iorb].real, hopping[ir, jorb, iorb].imag), file=f)


def ref2geom(filename):
    #
    print("  Convert from \"./dir-wfn/dat.lattice\" and \"./dir-wan/dat.wan-center\" to \"{0}\"".format(filename))
    #
    # Geometry (for HPhi and mVMC)
    #
    avec = numpy.zeros((3, 3), numpy.float_)
    with open("./dir-wfn/dat.lattice", 'r') as fi:
        for ii in range(3):
            line = fi.readline()
            itemlist = line.split()
            for jj in range(3):
                avec[ii, jj] = float(itemlist[jj])
    #
    bvec = numpy.linalg.inv(avec)
    #
    #
    with open("./dir-wan/dat.wan-center", 'r') as fi:
        for ii in range(2):
            line = fi.readline()  # skip
            print("    " + line, end="")
        line = fi.readlines()
    nwan = len(line)
    centre = numpy.zeros((nwan, 3), numpy.float_)
    for iwan in range(nwan):
        itemlist = line[iwan].split()
        for ii in range(3):
            centre[iwan, ii] = float(itemlist[ii])
    centre = numpy.dot(centre, bvec)
    #
    # Bohr -> Angstrom
    #
    avec[:, :] *= 0.529177249
    #
    with open(filename, 'w') as fo:
        for ii in range(3):
            print("%f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=fo)
        print("%d" % nwan, file=fo)
        for iwan in range(nwan):
            print("%f %f %f" % (centre[iwan, 0], centre[iwan, 1], centre[iwan, 2]), file=fo)


def respack2wan90(seedname):

    res2wan("./dir-wan/dat.h_mat_r", seedname + "_hr.dat")
    res2wan("./dir-intW/dat.Wmat", seedname + "_ur.dat")
    res2wan("./dir-intJ/dat.Jmat", seedname + "_jr.dat")

    ref2geom(seedname + "_geom.dat")


if __name__ == '__main__':

    args = sys.argv
    seedname = "zvo"

    if len(args) == 1:
        seedname = "zvo"
    elif len(args) == 2:
        seedname = args[1]
    else:
        print("\nUsage:\n")
        print("  $ respack2wan90.py [seedname]\n")
        exit(-1)

    respack2wan90(seedname)
