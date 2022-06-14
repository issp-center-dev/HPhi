#!/bin/sh
$1 -s stan1.in
python3 Finite.py
$1 -s stan2.in
python3 AveSSrand.py
$1 -s stan3.in
python3 AveFlct.py
