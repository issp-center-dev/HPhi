#!/bin/sh
$1 -s stan1.in
$1 -s stan2a.in
python3 AveSSrand.py
