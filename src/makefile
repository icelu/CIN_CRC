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

# all: simcrc simgland
all: simgland

# sveta: sveta.cpp
# 	cd gzstream/ && make
# 	#  Add -lintlc when using intel complier
# 	$(CCC) $(FLAG) -std=gnu++11 sveta.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp -o sveta -L$(BOOST)/lib/ -lboost_program_options -lgsl -lgslcblas -L./gzstream -lgzstream -lz -I$(BOOST)/include -I./gzstream

# simcrc:
# 	$(CCC) $(FLAG) -std=gnu++11 sim_crc.cpp -o ../bin/simcrc -L$(BOOST)/lib/ -lboost_program_options -lboost_system -lboost_filesystem -L$(GSL)/lib/ -lgsl -lgslcblas -lm -lz -I$(BOOST)/include -I$(GSL)/include

simgland:
	$(CCC) $(FLAG) -std=gnu++11 sim_gland.cpp -o ../bin/simgland -L$(BOOST)/lib/ -lboost_program_options -lboost_system -lboost_filesystem -L$(GSL)/lib/ -lgsl -lgslcblas -lm -lz -I$(BOOST)/include -I$(GSL)/include

clean:
	# rm ../bin/simcrc
	rm ../bin/simgland
