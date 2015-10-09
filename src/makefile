include make.sys

OBJS =	CG_EigenVector.c \
HPhiTrans.c \
HPhiMain.c \
FirstMultiply.c \
Lanczos_EigenValue.c \
Lanczos_EigenVector.c \
FileIO.c \
sz.c \
Multiply.c \
bisec.c \
check.c \
diagonalcalc.c \
dSFMT.o \
diagonalcalc.c \
expec_cisajscktaltdc.c \
expec_energy.c \
expec_totalspin.c \
FirstMultiply.c \
Lanczos_EigenValue.c \
Lanczos_EigenVector.c \
lapack_diag.c \
makeHam.c \
matrixlapack.c \
mfmemory.c \
mltply.c \
CalcByTPQ.c \
output.c \
output_list.c \
phys.c \
readdef.c \
sgn.c \
vec12.c \
xfreemem_large.c \
xsetmem.c \
log.c \

all:
	make libStdFace
	make HPhi

HPhi:$(OBJS) make.sys libStdFace.a
	$(CC) -g HPhiMain.c $(LAPACK_FLAGS) $(MTFLAGS) -I$(INCLUDE_DIR) dSFMT.o libStdFace.a -o $@ 

dSFMT.o:dSFMT.c make.sys
	$(CC) -c $< $(MTFLAGS) -I$(INCLUDE_DIR)

libStdFace:
	cd StdFace; make libStdFace.a

clean:
	cd StdFace; make clean
	rm -f HPhi *.o *.a
