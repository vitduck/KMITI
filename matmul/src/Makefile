CC = icc
CXXFLAGS = -g -std=c99

ifeq ($(REPORT), yes)
	CXXFLAGS+=-qopt-report=5
endif

ifeq ($(DOUBLE), yes) 
	CXXFLAGS+=-DDOUBLE
endif 

all: opt gemm

gemm:
	$(CC) $(CXXFLAGS) -O2 -qmkl=parallel -L$(MKLROOT)/lib/intel64 gemm.c -o gemm.x

no_loop-opt:
	$(CC) $(CXXFLAGS) -O1 matmul.c -o matmul_nlo.x

opt: 
	$(CC) $(CXXFLAGS) -O2 matmul.c -o matmul_opt.x

vec: 
	$(CC) $(CXXFLAGS) -O2 -xCORE-AVX512 -qopt-zmm-usage=high matmul.c -o matmul_vec.x

roofline:
	rm -rf ./roofline
	advisor -collect survey -project-dir ./roofline -- ./matmul.x 2000 2000 2000
	advisor -collect tripcounts -flop -project-dir ./roofline -- ./matmul.x 2000 2000 2000

gen9:
	advisor-python $(APM)/run_oa.py gen9 --config gen9_gt4 --collect basic --no-assume-dependencies -- ./matmul.x 2000 2000 2000

gen12: 
	rm -rf ./gen12
	advisor-python $(APM)/run_oa.py gen12 --config gen12_dg1 --collect full --no-assume-dependencies -- ./matmul.x 2000 2000 2000

clean: 
	rm -f *.x *.optrpt
