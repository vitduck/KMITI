CC = icc
COMPFLAGS = -g -std=gnu99 -O2

MKL_ROOT = /scratch/optpar01/apps/build/oneAPI/mkl/2021.3.0

CXXFLAGS = $(COMPFLAGS) 

ifeq ($(REPORT), yes)
	CXXFLAGS+=-qopt-report=5
endif

.SUFFIXES: .o .c

all: matmul dgemm

matmul: 
	$(CC) $(CXXFLAGS) matmul.c -o matmul.x

dgemm: 
	$(CC) $(CXXFLAGS) -mkl=parallel -L$(MKL_ROOT)/lib/intel64 dgemm.c -o dgemm.x

hotspot: 
	advixe-cl -collect survey -project-dir ./vtune -- ./matmul.x

memory: 
	advixe-cl -collect map -project-dir ./vtune -- ./matmul.x
	
roofline:
	advixe-cl -collect tripcounts -flop -project-dir ./vtune -- ./matmul.x

open-gui:
	advixe-gui vtune/vtune.advixeproj >/dev/null 2>&1 &

clean: 
	rm -f matmul.x dgemm.x *.optrpt

clean-all: 
	rm -f matmul.x dgemm.x *.optrpt 
	rm -rf ./vtune
