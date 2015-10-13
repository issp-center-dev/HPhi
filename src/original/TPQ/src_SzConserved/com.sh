#fccpx -Kfast,openmp -Nsrc,sta  EDmain.c -o ED
#fccpx -Dlapack -Kfast,openmp,SPARC64IXfx,parallel -Kmemalias,alias_const   EDmain.c -o ED -SSL2BLAMP
#icc -D"TPQ_"   -openmp -openmp-report2 -O3 EDmain.c -mkl=parallel -DDSFMT_MEXP=19937 -DHAVE_SSE2  dSFMT.c -o ED

icc  -openmp -openmp-report2 -O3 EDmain.c -DDSFMT_MEXP=19937 -DHAVE_SSE2  dSFMT.c -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread -o ED
#gcc -O3 -fopenmp -Wall EDmain.c -DDSFMT_MEXP=19937 -DHAVE_SSE2  dSFMT.c -o ED -llapack -lm
#gcc -g -fopenmp -Wall EDmain.c -DDSFMT_MEXP=19937 -DHAVE_SSE2  dSFMT.c -o ED -llapack -lm
#gcc -g -fopenmp EDmain.c -DDSFMT_MEXP=19937 -DHAVE_SSE2  dSFMT.c -o ED -llapack -lm 
