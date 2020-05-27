#!/bin/sh
python3 MakeDef.py
cp dir_input/*.def .
$1 -e namelist.def
