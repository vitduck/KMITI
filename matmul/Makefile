CC = icc
COMPFLAGS = -g -std=c11

CXXFLAGS = $(COMPFLAGS) 

ifeq ($(REPORT), yes)
	CXXFLAGS+=-qopt-report=5
endif

all: sgemm O1 novec block O2 block_O2

sgemm:
	$(CC) $(CXXFLAGS) -mkl=parallel -L$(MKLROOT)/lib/intel64 sgemm.c -o sgemm.x

O1:
	$(CC) $(CXXFLAGS) -O1 matmul.c -o matmul_O1.x

novec:
	$(CC) $(CXXFLAGS) -O2 -no-vec matmul.c -o matmul_novec.x

nounroll: 
	$(CC) $(CXXFLAGS) -O2 -no-vec -unroll0 matmul.c -o matmul_nounroll.x

block: 
	$(CC) $(CXXFLAGS) -O2 -no-vec matmul_block.c -o matmul_block.x

O2:
	$(CC) $(CXXFLAGS) -O2 -xCORE-AVX512 -qopt-zmm-usage=high matmul.c -o matmul_O2.x

block_O2: 
	$(CC) $(CXXFLAGS) -O2 -xCORE-AVX512 -qopt-zmm-usage=high matmul_block.c -o matmul_block_O2.x

roofline_O1:
	advisor -collect survey -project-dir ./roofline_O1 -- ./matmul_O1.x 2000 2000 2000
	advisor -collect tripcounts -flop --enable-cache-simulation -project-dir ./roofline_O1 -- ./matmul_O1.x 2000 2000 2000

roofline_novec:
	advisor -collect survey -project-dir ./roofline_novec -- ./matmul_novec.x 2000 2000 2000
	advisor -collect tripcounts -flop --enable-cache-simulation -project-dir ./roofline_novec -- ./matmul_novec.x 2000 2000 2000

roofline_block:
	advisor -collect survey -project-dir ./roofline_block -- ./matmul_block.x 2000 2000 2000
	advisor -collect tripcounts -flop --enable-cache-simulation -project-dir ./roofline_block -- ./matmul_block.x 2000 2000 2000

gen9:
	advisor-python $(APM)/run_oa.py gen9 --config gen9_gt4 --collect basic --no-assume-dependencies -- ./matmul.x 4000 4000 4000

gen11: 
	advisor-python $(APM)/run_oa.py gen11 --config gen11_icl --collect basic --no-assume-dependencies -- ./matmul.x 4000 4000 4000

gen12: 
	advisor-python $(APM)/run_oa.py gen12 --config gen12_dg1 --collect basic --no-assume-dependencies -- ./matmul.x 4000 4000 4000

gen12_config:
	advisor-python $(APM)/run_oa.py gen12_config --config gen12_dg1 --config scalers.toml --collect basic --no-assume-dependencies -- ./matmul.x 4000 4000 4000

open-gui:
	advisor-gui advisor/advisor.advixeproj >/dev/null 2>&1 &

clean: 
	rm -f *.x *.optrpt

clean-all: 
	rm -f *.x *.optrpt 
	rm -rf ./vtune
