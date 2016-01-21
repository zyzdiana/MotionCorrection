VERSION = DEBUG

ifeq ($(VERSION), DEBUG)                                                        
CXXFLAGS += -g
endif                                                                           
ifeq ($(VERSION), RELEASE)                                                      
CXXFLAGS += -O4
endif  

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)                                                         
CXX = clang++
LDFLAGS += -framework Accelerate
endif

all: 3d_interpolation test

3d_interpolation:

test: BinaryFile_tests.o

clean:
	rm -f *.o
	rm -f 3d_interpolation

PHONY: .clean
