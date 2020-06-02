find ./* -name "*.def"          -exec rm {} \;
find ./* -name "*.dat"          -exec rm {} \;
find ./* -name "*.gp"           -exec rm {} \;
find ./* -name "std.out"        -exec rm {} \;
find ./* -name "tmp_*"          -exec rm {} \;
find ./* -name "DG_*"           -exec rm -r {} \;
find ./* -name "output"         -exec rm -r {} \;
find ./* -name "dir_input"      -exec rm -r {} \;
find ./* -name "__pycache__"    -exec rm -r {} \;
