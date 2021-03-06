# detect OS
ifeq ($(OS), Windows_NT)
	detected_OS := Windows
else
	detected_OS := $(shell uname -s)
endif

MPICXX = mpic++
MPIRUN = mpirun

# define any compile-time flags
CFLAGS = -std=c++11 -O2

# Red-Hat 
DEFAULT_INCLUDE = /usr/include

ifeq ($(detected_OS), Linux)
	REDHAT_LIB = /usr/lib64
	EXT_REDHAT_LIB = -L$(REDHAT_LIB)
else
	REDHAT_LIB =
	EXT_REDHAT_LIB =
endif

# HDF5
ifeq ($(detected_OS), Darwin) 
	HDF_INSTALL = /usr/local/Cellar/hdf5/1.10.1_1
endif
ifeq ($(detected_OS), Linux)
	HDF_INSTALL = /data/yiyang.zhang/lib/hdf5-1.10.2/hdf5
endif

HDF_LIB = $(HDF_INSTALL)/lib
EXTLIB = -L$(HDF_LIB)
HDF_LIB_FILES = $(HDF_LIB)/libhdf5.a $(HDF_LIB)/libhdf5_hl.a

LIB = -lsz -lz -lm -ldl -lfftw3
# LIB = -lz -lm -ldl

# OPEN-MPI
MPI_INSTALL = /usr/local/Cellar/open-mpi/3.1.0
MPI_LIB = $(MPI_INSTALL)/lib

# ParaSite path
PARASITE_PATH = ./ParaSite

# Test path
TEST_PATH = ./tests

# Eigen path
EIGEN_PATH = ./eigen

# EW_Model path
EWMODEL_PATH = ./EW_Model

# define any directories containing header files other than /usr/include
#
INCLUDES = -I$(PARASITE_PATH) \
			-I$(TEST_PATH) \
			-I$(HDF_INSTALL)/include \
			-I$(MPI_INSTALL)/include \
			-I$(DEFAULT_INCLUDE)

LIBSHDF = $(EXTLIB) \
		$(EXT_REDHAT_LIB) \
		$(HDF_LIB_FILES)
		

# define the C source files
MPI_SRCS = $(PARASITE_PATH)/MPIWrapper2D.cpp
HDF5_SRCS = $(PARASITE_PATH)/HDF5Wrapper.cpp
FFT_SRCS = $(PARASITE_PATH)/FFTWrapper.cpp
TEST_SRCS = $(TEST_PATH)/tests.cpp


SRCS = main.cpp \
		$(MPI_SRCS) $(HDF5_SRCS) $(FFT_SRCS) \
		$(EWMODEL_PATH)/EW_helper.cpp \
		$(EWMODEL_PATH)/EW_examples.cpp
		# $(TEST_SRCS)

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.cpp=.o)

# define the executable file 
TARGET = emlat.out
OUTPUT = result.txt

all:$(TARGET)
	@echo Compile success.

$(TARGET): $(OBJS)
	$(MPICXX) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LIBSHDF) $(LIB)
.cpp.o:
	$(MPICXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

#$(PARASITE_PATH)/HDF5Wrapper.o: $(PARASITE_PATH)/HDF5Wrapper.cpp $(PARASITE_PATH)/HDF5Wrapper.h
#	$(MPICXX) $(CFLAGS) $(INCLUDES) -c $(PARASITE_PATH)/HDF5Wrapper.cpp

#$(PARASITE_PATH)/

run:$(TARGET)
	$(MPIRUN) -np 1 ./$(TARGET)
run1:$(TARGET)
	$(MPIRUN) -np 1 ./$(TARGET) -r 1 -c 1 &> $(OUTPUT) &
run2:$(TARGET)
	$(MPIRUN) -np 2 ./$(TARGET) -r 1 -c 2
run4:$(TARGET)
	$(MPIRUN) -np 4 ./$(TARGET) -r 2 -c 2 &> $(OUTPUT) &
run8:$(TARGET)
	$(MPIRUN) -np 8 ./$(TARGET) -r 2 -c 4 &> $(OUTPUT) &
run32:$(TARGET)
	$(MPIRUN) -np 32 ./$(TARGET) -r 4 -c 8 &> $(OUTPUT) &
clean:
	$(RM) $(OBJS) $(TARGET)
cleandata:
	$(RM) *.h5 result.txt *_dtable.txt *_param.txt