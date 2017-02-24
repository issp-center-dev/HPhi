#!/bin/bash
grep -A 3 "Lattice Vectors" ${1} | awk 'NR>1{print $2, $3, $4}'
grep -A 1000 " Final State" ${1} | awk '$1=="WF"&&$2=="centre"&&$3=="and"&&$4=="spread"{print $7, $8, $9}'|sed -e "s/,//g"
