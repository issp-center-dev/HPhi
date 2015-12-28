#!/usr/local/bin/perl
  $fname="calcmod.def";
  open(FILE,">$fname");
  printf FILE "CalcType $ARGV[0] \n";
  printf FILE "FlgFiniteTemperature 0 \n";
  printf FILE "CalcModel $ARGV[1]  \n";
  close(FILE);
