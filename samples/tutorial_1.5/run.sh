#!/bin/sh
$1 -sdry stan1.in
python3 MakeGreen.py
$1 -e namelist.def
python3 CalcSq.py 
