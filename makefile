MPICXX = mpic++
MPIRUN = mpirun

# define any compile-time flags
CFLAGS = -std=c++14 -O2

# HDF5 
HDF_INSTALL = /usr/local/Cellar/hdf5/1.10.1_1
HDF_LIB = $(HDF_INSTALL)/lib
EXTLIB = -L$(HDF_LIB)
LIB = -lsz -lz -lm

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
INCLUDES = -I.$(PARASITE_PATH) \
			-I.$(TEST_PATH) \
			-I$(HDF_INSTALL)/include

LIBSHDF = $(EXTLIB) $(HDF_LIB)/libhdf5.a $(HDF_LIB)/libhdf5_hl.a

# define the C source files
MPI_SRCS = $(PARASITE_PATH)/MPIWrapper2D.cpp
HDF5_SRCS = $(PARASITE_PATH)/HDF5Wrapper.cpp
TEST_SRCS = $(TEST_PATH)/tests.cpp

SRCS = main.cpp $(EWMODEL_PATH)/EW_helper.cpp $(MPI_SRCS) $(HDF5_SRCS)

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
	$(MPIRUN) -np 1 ./$(TARGET) -r 1 -c 1 > $(OUTPUT) &
run4:$(TARGET)
	$(MPIRUN) -np 4 ./$(TARGET) -r 2 -c 2 > $(OUTPUT) &
run8:$(TARGET)
	$(MPIRUN) -np 8 ./$(TARGET) -r 2 -c 4 > $(OUTPUT) &
run24:$(TARGET)
	$(MPIRUN) -np 24 ./$(TARGET) -r 4 -c 6 > $(OUTPUT) &
clean:
	$(RM) $(OBJS) $(TARGET)