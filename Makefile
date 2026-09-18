# Autodetec Compiler: Priority: g++ > clang++ > icpc. 
CXX_CANDIDATES := g++ clang++ icpc
DETECTED_CXX   := $(firstword $(foreach c,$(CXX_CANDIDATES),$(if $(shell command -v $(c) 2>/dev/null),$(c))))
ifeq ($(DETECTED_CXX),)
  $(warning No C++ compiler found among: $(CXX_CANDIDATES). Set CXX manually.)
endif
# CXX has a built-in default (origin "default"), so only override it there.
ifeq ($(origin CXX),default)
  CXX := $(DETECTED_CXX)
endif
$(info Using compiler: $(CXX))

CPPFLAGS := -Ofast -march=native -Wall -Wconversion -Wshadow -Wnon-virtual-dtor -std=c++20 -fPIC -fopenmp
LD       := $(CXX)
LDFLAGS  :=

# BLAS/LAPACK: MKL > OpenBLAS > reference BLAS/LAPACK. Detected by test-linking
# each candidate. Override with `make LDLIBS=...`.
MKL_LIBS      := -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
OPENBLAS_LIBS := -lopenblas
BLAS_LIBS     := -lblas -llapack
check_link = $(shell echo 'int main(){return 0;}' | $(CXX) -x c++ - -o /dev/null $(1) $(LDFLAGS) >/dev/null 2>&1 && echo yes)

ifeq ($(call check_link,$(MKL_LIBS)),yes)
  BACKEND_LIBS := $(MKL_LIBS)
  BACKEND_NAME := Intel MKL
else ifeq ($(call check_link,$(OPENBLAS_LIBS)),yes)
  BACKEND_LIBS := $(OPENBLAS_LIBS)
  BACKEND_NAME := OpenBLAS
else ifeq ($(call check_link,$(BLAS_LIBS)),yes)
  BACKEND_LIBS := $(BLAS_LIBS)
  BACKEND_NAME := reference BLAS/LAPACK
else
  BACKEND_LIBS :=
  BACKEND_NAME := none found (set LDLIBS manually)
endif
LDLIBS ?= $(BACKEND_LIBS)
$(info Using BLAS/LAPACK backend: $(BACKEND_NAME))

# OpenMP runtime: libiomp5 (Intel) > libgomp (GNU). Detected the same way.
IOMP_LIBS := -liomp5 -lpthread
GOMP_LIBS := -lgomp -lpthread
ifeq ($(call check_link,$(IOMP_LIBS)),yes)
  OMP_LIBS := $(IOMP_LIBS)
  OMP_NAME := Intel libiomp5
else ifeq ($(call check_link,$(GOMP_LIBS)),yes)
  OMP_LIBS := $(GOMP_LIBS)
  OMP_NAME := GNU libgomp
else
  OMP_LIBS :=
  OMP_NAME := none found (set OMP_LIBS manually)
endif
LDLIBS += $(OMP_LIBS)
$(info Using OpenMP runtime: $(OMP_NAME))

# Files with their own int main() each become a binary; the rest is common code.
MAIN_SRCS   := main.cpp 3c_2c_transform.cpp occ_ri_k_step.cpp
ALL_SRCS    := $(wildcard *.cpp)
COMMON_SRCS := $(filter-out $(MAIN_SRCS),$(ALL_SRCS))
COMMON_OBJS := $(COMMON_SRCS:.cpp=.o)
MAIN_OBJS   := $(MAIN_SRCS:.cpp=.o)
BINS        := $(MAIN_SRCS:.cpp=)
HEADERS     := $(wildcard *.h) $(wildcard *.hpp)

.PHONY: all clean test check

all: $(BINS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) -c $< -o $@

$(BINS): %: %.o $(COMMON_OBJS)
	$(LD) -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean:
	rm -f $(MAIN_OBJS) $(COMMON_OBJS) $(BINS)
	rm -f *.a *.so
	rm -f tests/*.o tests/run_tests

# ---- test suite ---------------------------------------------------------
# All tests/*.cpp are compiled into a single run_tests binary and linked
# against the library's own COMMON_OBJS, so tests exercise the exact same
# compiled code as the real binaries (not a reimplementation). doctest.h is
# only included by tests/*.cpp, so it never touches the main build above.
TEST_SRCS    := $(wildcard tests/*.cpp)
TEST_OBJS    := $(TEST_SRCS:.cpp=.o)
TEST_HEADERS := $(wildcard tests/*.hpp) tests/thirdparty/doctest.h

tests/%.o: tests/%.cpp $(HEADERS) $(TEST_HEADERS)
	$(CXX) $(CPPFLAGS) -Itests/thirdparty -c $< -o $@

tests/run_tests: $(TEST_OBJS) $(COMMON_OBJS)
	$(LD) -o $@ $^ $(LDLIBS) $(LDFLAGS)

test check: tests/run_tests
	./tests/run_tests
