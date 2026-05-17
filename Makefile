CC = gcc
CXX = g++ # requires g++-9 or newer / icpc (with gcc compatibility 7.5 or newer) / clang++ with llvm-10 or newer
CXXFLAGS = -std=c++17 -fopenmp -Wall -Wfloat-conversion # need C++17 and OpenMP

SCTL_DIR=extern/SCTL
BIEST_DIR=extern/BIEST
CXXFLAGS += -I${SCTL_DIR}/include -I${BIEST_DIR}/include

#Optional flags
DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CXXFLAGS += -O0 -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract,float-divide-by-zero,float-cast-overflow -fno-sanitize-recover=all -fstack-protector # debug build
	CXXFLAGS += -DSCTL_MEMDEBUG # Enable memory checks
else
	CXXFLAGS += -O3 -march=native -DNDEBUG # release build
endif

OS = $(shell uname -s)
ifeq "$(OS)" "Darwin"
	CXXFLAGS += -g -rdynamic -Wl,-no_pie # for stack trace (on Mac)
else
	CXXFLAGS += -gdwarf-4 -g -rdynamic # for stack trace
endif

CXXFLAGS += -DSCTL_PROFILE=5 -DSCTL_VERBOSE # Enable profiling
CXXFLAGS += -DSCTL_SIG_HANDLER

#CXXFLAGS += -DSCTL_QUAD_T=_Float128 # Enable quadruple precision (GCC)
#CXXFLAGS += -Qoption,cpp,--extended_float_types -DSCTL_QUAD_T=__float128 # Enable quadruple precision (Intel)

#CXXFLAGS += -lblas -DSCTL_HAVE_BLAS # use BLAS
#CXXFLAGS += -llapack -DSCTL_HAVE_LAPACK # use LAPACK
#CXXFLAGS += -lopenblas -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK
#CXXFLAGS += -qmkl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK -DSCTL_HAVE_FFTW3_MKL # use MKL BLAS, LAPACK and FFTW (Intel compiler)
CXXFLAGS += -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK # use MKL BLAS and LAPACK (non-Intel compiler)
#CXXFLAGS += -DSCTL_HAVE_SVML

CXXFLAGS += -lfftw3_omp -DSCTL_FFTW_THREADS
CXXFLAGS += -lfftw3 -DSCTL_HAVE_FFTW
CXXFLAGS += -lfftw3f -DSCTL_HAVE_FFTWF
CXXFLAGS += -lfftw3l -DSCTL_HAVE_FFTWL

CXXFLAGS += -DSCTL_HAVE_LIBMVEC
#CXXFLAGS += -DSCTL_HAVE_SVML

RM = rm -f
MKDIRS = mkdir -p

BINDIR = ./bin
LIBDIR = ./lib
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include

TARGET_LIB = $(LIBDIR)/libvirtualcasing.a
TARGET_BIN = $(BINDIR)/virtual-casing \
						 $(BINDIR)/virtual-casing-c \
						 $(BINDIR)/axisymm

.PHONY: all clean

all : $(TARGET_BIN) $(TARGET_LIB)

$(BINDIR)/%: ./test/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $^ -o $@
ifeq "$(OS)" "Darwin"
	/usr/bin/dsymutil $@ -o $@.dSYM
endif

$(BINDIR)/%: ./test/%.c $(TARGET_LIB)
	-@$(MKDIRS) $(dir $@)
	$(CC) $(CXXFLAGS) -I$(INCDIR) $(TARGET_LIB) -lm -ldl -lstdc++ $^ -o $@
ifeq "$(OS)" "Darwin"
	/usr/bin/dsymutil $@ -o $@.dSYM
endif

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $^ -o $@

$(TARGET_LIB) : obj/virtual-casing.o
	-@$(MKDIRS) $(dir $@)
	ar rcs $@ $^

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/* $(LIBDIR)/*

