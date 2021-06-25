EXEC=rydberg
CC = icpc
SRCDIR = src/
LOCALLIB = /usr/local/lib
LOCALINCLUDE = /usr/local/include
PROJINCLUDE = inc/
PROJLIB = lib/
MKLROOT = /opt/intel/oneapi/mkl/latest/
MKLINCLUDE = /opt/intel/oneapi/mkl/latest/include
INTELLIB = /opt/intel/oneapi/mkl/latest/lib


SRCEXT := cpp
SRC_FILES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
CFLAGS = -c -g -fvar-tracking -traceback -Wall -DMKL_ILP64 -qopenmp -fast -O3 -xhost -ip -fbuiltin -ipo -no-ftz -static-intel -std=c++11 -qmkl=parallel -Wl,-rpath,$(MKLROOT)/lib -I$(MKLINCLUDE) -I$(LOCALINCLUDE) -I$(PROJINCLUDE)
LFLAGS = -Wl,-rpath,$(MKLROOT)/lib -L$(MKLROOT)/Lib -L$(LOCALLIB) -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -L$(INTELLIB) -liomp5 -lpthread -lm

O_FILES = $(SRC_FILES:.cpp=.o) $(INI_SRC:.cpp=.o) 

print-%  : ; @echo $* = $($*)

$(EXEC): $(O_FILES)
	$(CC) -o $@ $(O_FILES) $(LFLAGS)
	
%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
