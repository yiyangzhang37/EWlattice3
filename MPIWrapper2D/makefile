MPICXX = mpic++
MPIRUN = mpirun

# define any compile-time flags
CFLAGS = -std=c++14

# HDF5 
HDF_INSTALL = /usr/local/Cellar/hdf5/1.10.1_1
HDF_LIB = $(HDF_INSTALL)/lib
EXTLIB = -L$(HDF_LIB)
LIB = -lsz -lz -lm

# define any directories containing header files other than /usr/include
#
INCLUDES = -I./ParaSite -I./tests \
			-I$(HDF_INSTALL)/include

LIBSHDF = $(EXTLIB) $(HDF_LIB)/libhdf5.a $(HDF_LIB)/libhdf5_hl.a

# define the C source files
MPI_SRCS = ./ParaSite/MPIWrapper2D.cpp
HDF5_SRCS = ./ParaSite/HDF5Wrapper.cpp
TEST_SRCS = ./tests/Test_MPIWrapper2D.cpp
SRCS = main.cpp EW_helper.cpp $(MPI_SRCS) $(TEST_SRCS) $(HDF5_SRCS)

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
TARGET = lat_test
OUTPUT = lat_test.txt

all:$(TARGET)
	@echo Compile success.

$(TARGET): $(OBJS)
	$(MPICXX) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LIBSHDF) $(LIB)
.cpp.o:
	$(MPICXX) $(CFLAGS) $(INCLUDES) -c $< -o $@


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