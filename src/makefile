# CCC = icc		# Using intel complier
CCC = g++
# CCC = /usr/local/opt/llvm/bin/clang
BOOST = /usr/local
# GSL = /share/apps/gsl-2.5
GSL = /usr/local
# omp = -fopenmp
omp =
FLAG = -O3
# -pg
# FLAG =

OUT_DIR = ../bin

MKDIR_P = mkdir -p

.PHONY: directories


# all: directories simcrc simgland
all: directories simgland simcell


directories: ${OUT_DIR}


${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}


# simcrc:
# 	$(CCC) $(FLAG) -std=gnu++11 sim_crc.cpp -o ../bin/simcrc -L$(BOOST)/lib/ -lboost_program_options -lboost_system -lboost_filesystem -L$(GSL)/lib/ -lgsl -lgslcblas -lm -lz -I$(BOOST)/include -I$(GSL)/include

simgland:
	$(CCC) $(FLAG) -std=gnu++11 sim_gland.cpp -o ${OUT_DIR}/simgland -L$(BOOST)/lib/ -lboost_program_options -lboost_system -lboost_filesystem -L$(GSL)/lib/ -lgsl -lgslcblas -lm -lz -I$(BOOST)/include -I$(GSL)/include

simcell:
	$(CCC) $(FLAG) -std=gnu++11 sim_single_gland.cpp -o ${OUT_DIR}/simcell -L$(BOOST)/lib/ -lboost_program_options -lboost_system -lboost_filesystem -L$(GSL)/lib/ -lgsl -lgslcblas -lm -lz -I$(BOOST)/include -I$(GSL)/include


clean:
	# rm ../bin/simcrc
	rm ../bin/simgland
