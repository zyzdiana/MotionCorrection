VERSION = DEBUG
#VERSION = RELEASE

CXXFLAGS += -I./
#CXXFLAGS += -I/Users/zyzdiana/GitHub/MotionCorrection/

ifeq ($(VERSION), DEBUG)                                                        
CXXFLAGS += -g
LDFLAGS += -g
endif                                                                           
ifeq ($(VERSION), RELEASE)                                                      
CXXFLAGS += -g
CXXFLAGS += -O3
LDFLAGS += -O3
endif  

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)                                                         
CXX = clang++
LDFLAGS += -framework Accelerate
CXXFLAGS += -DOSX
endif

ifeq ($(UNAME), Linux)
CXX = g++
CXXFLAGS += -DLINUX
endif


all: test test_read_file test_interp test_time_profiler test_rotate_coords

test_read_file:

test_interp:

test_time_profiler:

test_inverse:

test_rotate_coords:

test_tricubic:

test_gn:

#test: BinaryFile_tests.o interp3D_tests.o
TESTOBJECTS += BinaryFile_tests.o
TESTOBJECTS += Volume_tests.o
TESTOBJECTS += CentralDifferenceDifferentiator_tests.o
TESTOBJECTS += TrilinearInterpolator_tests.o

test: $(TESTOBJECTS)

clean:
	rm -f *.o
	rm -f test
	rm -f test_read_file
	rm -f test_interp
	rm -f test_time_profiler
	rm -f test_tricubic
	rm -f test_gn

PHONY: .clean
