#!/bin/sh
source ~/.bashrc
echo "start to calculate eigenvec"
./HPhi -e ./namelist_calceigenvec.def > output_eigenvec.out
echo "start to output information"
./HPhi -e ./namelist_outputinfo.def > outputinfo.out
echo "start to restart spectrum"
./HPhi -e ./namelist_restart.def > outputrestart.out
echo "start to restart2 spectrum"
./HPhi -e ./namelist_restart.def > outputrestart.out
