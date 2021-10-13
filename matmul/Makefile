CC = icc
COMPFLAGS = -g -std=c99 -O2 -xCORE-AVX512

CXXFLAGS = $(COMPFLAGS) 

ifeq ($(REPORT), yes)
	CXXFLAGS+=-qopt-report=5
endif

all: matmul dgemm

matmul:
	$(CC) $(CXXFLAGS) matmul.c -o matmul.x

dgemm:
	$(CC) $(CXXFLAGS) -mkl=parallel -L$(MKLROOT)/lib/intel64 dgemm.c -o dgemm.x

roofline:
	rm -rf ./roofline
	advisor -collect survey -project-dir ./roofline -- ./matmul.x 2000 2000 2000
	advisor -collect tripcounts -flop --enable-cache-simulation -project-dir ./roofline -- ./matmul.x 2000 2000 2000

dependencies:
	rm -rf ./dep
	advisor --collect=survey --project-dir=./dep -- ./matmul.x 2000 2000 2000
	advisor --collect=dependencies --project-dir=./dep --select=matmul.c:71,matmul.c:83 -- ./matmul.x 2000 2000 2000

offload:
	rm -rf ./offload
	advisor-python $(APM)/run_oa.py offload --config gen9_gt4 --collect basic -- ./matmul.x 2000 2000 2000

gen9:
	rm -rf ./gen9
	advisor-python $(APM)/run_oa.py gen9 --config gen9_gt4 --collect full --no-assume-dependencies -- ./matmul.x 2000 2000 2000

gen11: 
	rm -rf ./gen11
	advisor-python $(APM)/run_oa.py gen11 --config gen11_icl --collect full --no-assume-dependencies -- ./matmul.x 2000 2000 2000

gen12: 
	rm -rf ./gen12
	advisor-python $(APM)/run_oa.py gen12 --config gen12_dg1 --collect full --no-assume-dependencies -- ./matmul.x 2000 2000 2000

gen12_config:
	rm -rf ./gen12_config
	advisor-python $(APM)/run_oa.py gen12_config --config gen12_dg1 --config scalers.toml --collect full --no-assume-dependencies -- ./matmul.x 2000 2000 2000

open-gui:
	advisor-gui advisor/advisor.advixeproj >/dev/null 2>&1 &

clean: 
	rm -f *.x *.optrpt

clean-all: 
	rm -f *.x *.optrpt 
	rm -rf ./vtune
